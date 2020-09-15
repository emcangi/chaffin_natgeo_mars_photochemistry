__precompile__()

module Analysis

using PyPlot
using PyCall
using HDF5, JLD
using LaTeXStrings
using Distributed
using DelimitedFiles
using SparseArrays
using LinearAlgebra
using PlotUtils
using Photochemistry: get_ncurrent, write_ncurrent, n_tot, effusion_velocity, 
                      speciesbcs, plot_bg, plotatm, Tpiecewise, Psat, Psat_HDO, 
                      searchdir, search_subfolders, create_folder, input

include("../../../PARAMETERS.jl")

export get_ncurrent, write_ncurrent, n_tot, 
       effusion_velocity, speciesbcs, 
       areadensity_to_micron_atm, molec_to_GEL, GEL_to_molecule, 
       get_flux, calculate_f,
       search_subfolders, create_folder, searchdir, input,
       get_colors, get_grad_colors, plot_bg, plotatm,
       Tpiecewise, Psat, Psat_HDO


# Conversions ==================================================================
function areadensity_to_micron_atm(numpercm2)
    #=
    Conversion: num/cm^2 * (cm^2/m^2) * (10μm-am/2.687e20 #/m^2)
    =#
    return numpercm2 * (1e4/1) * (10/2.687e20)
end

function molec_to_GEL(molecules, HorH2O)
    #=
    Converts molecules of H2O per cm^2 to a global equivalent layer in meters. 
    molecules: number of molecules of species HorH2O
    HorH2O: string, either "H" or "H2O" to specify which type the molecule is
    =#
    if HorH2O == "H"
        molecules = molecules / 2
    end
        
    return molecules * 1e4 * 18 * 1.67e-27 / 999.89
end

function GEL_to_molecule(GEL, HorH2O)
    #=
    Converts a global equivalent layer in meters to molecules per cm^2. 
    GEL: layer of water in meters
    HorH2O: string, either "H" or "H2O" to specify which molecule to return
    =#
    
    molec_H2O = (GEL * 999.89) / (1e4 * 18 * 1.67e-27)
    
    if HorH2O == "H"
        return 2 * molec_H2O
    elseif HorH2O == "H2O"
        return molec_H2O
    end
end

# Flux and f ===================================================================
function get_flux(species, readfile, temps; oflux=1.2e8, repro=false, therm_only=false)
    #=
    Retrieves the flux for either H or D at the top of the equilibrated 
    atmosphere.

    NOT TO BE CONFUSED WITH FUNCTION getflux DEFINED IN PHOTOCHEMISTRY MODULE!!

    species: species in question, :H or :D. no error control right now
    readfile: the file with simulation results
    oflux: flux of O in /cm^2s. 
    temps: array of [Ts, Tt, Te]
    repro: whether flux is being calculated for reproduction of a past study
    therm_only: whether to return flux_t only, false by default
    =#
    n_current = get_ncurrent(readfile)

    # find the experiment type for the purpose of setting the boundary 
    # conditions correcftly 
    exptype = match(r"[a-z]{0,5}(?=_.+)",readfile).match

    # the species which can have ither D or H in them for loss.
    bearer = Dict(:D=>[:D, :HD], :H=>[:H, :HD, :H2])
    num_D_or_H = Dict(:D=>[1, 1], :H=>[1, 1, 2])

    # this dict keeps track of loss due to each species. order: H, D, H2, HD
    contrib_t = Dict(:H=>0., :D=>0., :H2=>0., :HD=>0.)
    contrib_nt = Dict(:H=>0., :D=>0., :H2=>0., :HD=>0.)

    # total thermal and non-thermal fluxes, separated.
    flux_t = 0
    flux_nt = 0

    # set things up for accurate boundary conditions
    Temp(z::Float64) = Tpiecewise(z, temps[1], temps[2], temps[3])
    Temp_keepSVP(z::Float64) = Tpiecewise(z, meanTs, meanTt, meanTe)
    if exptype=="temp"
        H2Osat = map(x->Psat(x), map(Temp_keepSVP, alt))
        HDOsat = map(x->Psat_HDO(x), map(Temp_keepSVP, alt))
    else
        H2Osat = map(x->Psat(x), map(Temp, alt))
        HDOsat = map(x->Psat_HDO(x), map(Temp, alt))
    end
    surface_watersat = Dict("H2O"=>H2Osat[1], "HDO"=>HDOsat[1])
    H_effusion_velocity = effusion_velocity(Temp(zmax), 1.0, zmax)
    H2_effusion_velocity = effusion_velocity(Temp(zmax), 2.0, zmax)
    D_effusion_velocity = effusion_velocity(Temp(zmax), 2.0, zmax)
    HD_effusion_velocity = effusion_velocity(Temp(zmax), 3.0, zmax)

    # Used for passing a variable speciesbcs function
    v_eff = Dict("H"=>H_effusion_velocity, "D"=>D_effusion_velocity, 
                 "H2"=>H2_effusion_velocity, "HD"=>HD_effusion_velocity)

    # Calculate the thermal escape
    for (s, m) in zip(bearer[species], num_D_or_H[species])
        this_species_t_flux = m * n_current[s][end]*speciesbcs(s, surface_watersat, v_eff, oflux)[2,2]
        flux_t += this_species_t_flux
        contrib_t[s] += this_species_t_flux
    end
    
        # If inclusion of non-thermal escape is requested (therm_only==false), this section
    # will add the effect of non-thermal escape. 
    if therm_only==false
        if repro==false
            # Nonthermal ecsape velocities for temperatures: T_exo = 150K, 205K, 250K,
            # using ratios of thermal/nonthermal from Kras 2002, in cm/s.
            inds = Dict(150=>1, 205=>2, 250=>3) # indices for different exobase temps
            i = inds[Int(temps[3])]             # convert exobase temp to an index 
            v_nt = Dict(:H => [26.7, 34.4, 45.1], :H2 => [8.05, 10.2, 13.26], # assuming 2nd order polynomial
                        :D => [10.65, 13.47, 17.47], :HD => [5.65, 6.00, 7.42])  # in cm/s.
        else
            # If reproducing past studies, need the exact nonthermal escape velocities from Kras 2002.
            inds = Dict(200 => 1, 270 => 2, 350 => 3)
            i = inds[temps[3]]
            # Nonthermal: cm/s. Each species has a value recorded at T = 200K, 
            # 270K, and 350K.
            v_nt = Dict(:H => [38, 56, 89], :H2 => [12.9, 18.2, 28], :D => [17, 24, 37],
                        :HD => [8.2, 11.5, 17.7])  
        end

        for (s, m) in zip(bearer[species], num_D_or_H[species])
            this_species_nt_flux = m * n_current[s][end] * v_nt[s][i]
            flux_nt += this_species_nt_flux
            contrib_nt[s] += this_species_nt_flux
        end
    end

    if therm_only==true
        return flux_t, contrib_t
    else
        return flux_t, flux_nt, contrib_t, contrib_nt
    end
end

function calculate_f(thefile, flux_type, temps; oflux=1.2e8, reprod=false)
    #=
    A function to calculate f or a single simulation.

    thefile: an equilibrated atmosphere simulation for which to calculate f.
    flux_type: "thermal", "nonthermal", or "both"
    temps: list of the 3 parameters [Tsurf, Ttropo, Texo]
    oflux: escape flux of atomic O if varying. 
    reprod: whether we are calculating for a reproduction of an earlier study.
    =#
    ncur = get_ncurrent(thefile)

    # contrib dictionaries (of how each bearer species contributes to escape)
    # are not used in this function.

    t_flux_H, nt_flux_H, contrib_t_H, contrib_nt_H = get_flux(:H, thefile, temps, repro=reprod)
    t_flux_D, nt_flux_D, contrib_t_D, contrib_nt_D = get_flux(:D, thefile, temps, repro=reprod)

    if flux_type=="thermal"
        Hf = t_flux_H
        Df = t_flux_D
    elseif flux_type=="both"
        Hf = t_flux_H + nt_flux_H
        Df = t_flux_D + nt_flux_D
    elseif flux_type=="nonthermal"
        Hf = nt_flux_H
        Df = nt_flux_D
    else
        println("Invalid escape type: $(flux_type)")
    end
  
    return 2*(Df/Hf) / (ncur[:HDO][1]/ncur[:H2O][1])
end

# Plotting Functions ============================================================
function get_colors(L, cmap)
    #=
    Generates some colors based on a non-gradient color map for use in plotting a 
    bunch of lines all at once.
    L: number of colors to generate.
    cmap: color map name

    NOTE: This function refuses to work unless "using PlotUtils" is in the master
          module file. It's not enough to have "using PlotUtils" appear in the 
          same submodule file in which this function is defined for some reason.
          Thus I have this function here in the main file.
    =#

    Cmap = get_cmap(Symbol(cmap))
    colors_dumb = [Cmap(x) for x in range(0, stop=1, length=L)]
    c = Array{Float64}(undef, L, 3)

    for i in range(1, length=length(colors_dumb))
        c[i, 1] = colors_dumb[i][1]
        c[i, 2] = colors_dumb[i][2]
        c[i, 3] = colors_dumb[i][3]
    end
    return c
end

function get_grad_colors(L, cmap; strt=0, stp=1)
    #=
    Generates some colors based on a GRADIENT color map for use in plotting a 
    bunch of lines all at once.
    L: number of colors to generate.
    cmap: color map name

    AVAILABLE MAPS: blues, viridis, pu_or, magma, plasma, inferno
    =#

    colors_dumb = [cgrad(Symbol(cmap))[x] for x in range(strt, stop=stp, length=L)]
    c = Array{Float64}(undef, L, 3)

    for i in range(1, length=length(colors_dumb))
        c[i, 1] = red(colors_dumb[i])
        c[i, 2] = green(colors_dumb[i])
        c[i, 3] = blue(colors_dumb[i])
    end
    return c
end

end