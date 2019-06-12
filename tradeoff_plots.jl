################################################################################
# make_tradeoff_plots.jl - make trade off plots showing H and D flux and other
# variables as a function of random things such as water column. 

# Eryn Cangi
# 5 April 2019
# Currently tested for Julia: 0.7
################################################################################
using PyPlot
using HDF5
using LaTeXStrings
using PyCall
using PlotUtils
using JLD

# fundamental constants ========================================================
const boltzmannK = 1.38e-23;    # J/K
const bigG = 6.67e-11;          # N m^2/kg^2
const mH = 1.67e-27;            # kg
const marsM = 0.1075*5.972e24;  # kg
const radiusM = 3396e5;         # cm
DH = 5.5 * 1.6e-4               # SMOW value from Yung 1988

# Basic altitude stuff
alt = (0:2e5:200e5)
n_alt_index=Dict([z=>clamp((i-1),1, length(alt)-2) for (i, z) in enumerate(alt)])
const zmax = alt[end];

# Folder Location ==============================================================
# lead = "/data/GDrive-CU/Research/Results/TradeoffPlots/"
lead = "/home/emc/GDrive-CU/Research/Results/TradeoffPlots/"

# Functions ====================================================================
function Tpiecewise(z::Float64, Tsurf, Ttropo, Tinf, lapserate=-1.4e-5, ztropo=90e5)
    #=
    DO NOT MODIFY! If you want to change the temperature, define a
    new function or select different arguments and pass to Temp(z)

    a piecewise function for temperature as a function of altitude,
    using Krasnopolsky's 2010 temperatures for altitudes
    >htropo=90km, fixed at Ttropo=125K between htropo and
    htropo-htropowidth=60km, and rising at a constant lapse rate
    (1.4K/km) below.

    z: an altitude in cm.
    returns: temperature in K at the requested altitude.
    =#
    # allow varying troposphere width
    ztropowidth = ztropo - (1/lapserate)*(Ttropo-Tsurf)
    if ztropowidth < 0   # in case a profile is nonphysical, let lapse rate vary
        ztropowidth = 30e5
        m = (ztropo/1e5 - ztropowidth/1e5) / (Ttropo - Tsurf)
        lapserate = (1/m)*1e-5
    end

    if z >= ztropo
        return Tinf - (Tinf - Ttropo)*exp(-((z-ztropo)^2)/(8e10*Tinf))
    end
    if ztropo > z >= ztropo - ztropowidth
        return Ttropo
    end
    if ztropo-ztropowidth > z
        return Ttropo-lapserate*(ztropo-ztropowidth-z)
    end
end

function get_ncurrent(readfile)
    n_current_tag_list = map(Symbol, h5read(readfile,"n_current/species"))
    n_current_mat = h5read(readfile,"n_current/n_current_mat");
    n_current = Dict{Symbol, Array{Float64, 1}}()
    for ispecies in [1:length(n_current_tag_list);]
        n_current[n_current_tag_list[ispecies]] = reshape(n_current_mat[:,ispecies],length(alt)-2)
    end
    n_current
end

function get_H_fluxes(readfile, oflux, temps)
    #=
    Produces an array of H and D fluxes at the top of the atmosphere (200 km) due to
    introduction of water parcels of various ppm and various altitudes.
    =#

    n_current = get_ncurrent(readfile)

    # Calculate the D flux: sum of (H population @ 200 km) * (D flux rate)
    # and (H2 population @ 200 km) * 2*(H2 flux rate) and 
    # (HD population @ 200 km) * (HD flux rate)
    Hfluxes = (n_current[:H][end]*speciesbcs_VARIABLE(:H, oflux, temps)[2,2]
                  + 2*n_current[:H2][end]*speciesbcs_VARIABLE(:H2, oflux, temps)[2,2]
                  + n_current[:HD][end]*speciesbcs_VARIABLE(:HD, oflux, temps)[2,2])
    return Hfluxes
end

function get_D_fluxes(readfile, oflux, temps)
    #=
    Gets value of D flux out of the atmosphere in a converged file.
    =#
    n_current = get_ncurrent(readfile)

    # Calculate the D flux: sum of (D population @ 200 km) * (D flux rate)
    # and (HD population @ 200 km) * (HD flux rate)
    Dfluxes = (n_current[:D][end]*speciesbcs_VARIABLE(:D, oflux, temps)[2,2]
                  + n_current[:HD][end]*speciesbcs_VARIABLE(:HD, oflux, temps)[2,2])
    return Dfluxes
end

function effusion_velocity(Texo::Float64, m::Float64)
    #=
    Returns effusion velocity for a species.
    Texo: temperature of the exobase (upper boundary) in K
    m: mass of one molecule of species in amu
    =#
    lambda = (m*mH*bigG*marsM)/(boltzmannK*Texo*1e-2*(radiusM+zmax))
    vth=sqrt(2*boltzmannK*Texo/(m*mH))
    v = 1e2*exp(-lambda)*vth*(lambda+1)/(2*pi^0.5)
    return v
end

function n_tot(n_current, z)
    specieslist = [:CO2, :O2, :O3, :H2, :OH, :HO2, :H2O, :H2O2, :O, :CO,
                   :O1D, :H, :N2, :Ar, :CO2pl, :HOCO,
                   # species for deuterium chemistry:
                   :HDO, :OD, :HDO2, :D, :DO2, :HD, :DOCO];
    # get the total number density at a given altitude
    thisaltindex = n_alt_index[z]
    return sum( [n_current[s][thisaltindex] for s in specieslist] )
end

function LHDO(T)
   #=
   Latent heat of vaporization of HDO as a function of temperature in K.
   This analytical function was determined by fitting data from
   https://pubs.acs.org/doi/pdf/10.1021/cr60292a004 (Jancso 1974, page 734)
   and extrapolating. It is probably not accurate outside the range of the data,
   which was 0-100, but it shouldn't be too far off.

   The data was in cal/mol versus Celsius. We convert to kJ/mol below.
   Fit was done in Python in a separate script. parameters below are the output
   from said fit.
   =#
   a = -0.02806171415983739
   b = 89.51209910268079
   c = 11918.608639939 # this is cal/mol
   return (a*(T-b)^2 + c) / 239 # returns in kJ/mol
end

function Psat_HDO(T)
    #=
    Analytical expression for saturation vapor pressure of HDO, using analytical
    latent heat for HDO determined from fitting to data from Jancso 1974.
    The boiling temperature for HDO (100.7°C, 373.85K) and standard pressure are
    used as reference temp/pressure. The gas constant for HDO was calculated as
    R_HDO = kB/m, which probably came from Pierrehumbert's text.

    returns the pressure in #/cm^3. The conversion from N/m^2 to #/cm^3 is 
    divide by N*m (that is, kT) and multiply by the conversion to cm from m 
    (1e-6). 

    Input
        T: a single temperature in K
    Output:
        Pressure in #/cm^3
   =#
   R_HDO = 434.8 # J/kgK
   L_HDO = LHDO(T) * 1000 / 0.019 # kJ/mol * J/kJ * mol / kg = J / kg
   T0 = 373.85 # boiling temp for liquid HDO
   P0 = 101325 # standard pressure
   Psat_pascals = P0*exp(-(L_HDO/R_HDO)*(1/T - 1/T0))
   # To convert to #/cm^3 from Pa:
   # Pa = Nm^-2 (1/(Nm)) * (1e-6 m^3 / 1 cm^3)
   return Psat_pascals*(1e-6)/(boltzmannK*T) # return in #/cm^3
end

Psat(T::Float64) = (133.3/(10^6 * boltzmannK * T))*(10^(-2445.5646/T 
                    + 8.2312*log10(T) - 0.01677006*T + 1.20514e-5*T^2 
                    - 6.757169))

function speciesbcs_VARIABLE(species, oflux, temps)

    Temp(z::Float64) = Tpiecewise(z, temps[1], temps[2], temps[3])
    H2Osat = map(x->Psat(x), map(Temp, alt))
    HDOsat = map(x->Psat_HDO(x), map(Temp, alt))
    H_effusion_velocity = effusion_velocity(Temp(zmax),1.0)
    H2_effusion_velocity = effusion_velocity(Temp(zmax),2.0)
    D_effusion_velocity = effusion_velocity(Temp(zmax),2.0)
    HD_effusion_velocity = effusion_velocity(Temp(zmax),3.0)

    speciesbclist=Dict(
                    :CO2=>["n" 2.1e17; "f" 0.],
                    :Ar=>["n" 2.0e-2*2.1e17; "f" 0.],
                    :N2=>["n" 1.9e-2*2.1e17; "f" 0.],
                    :H2O=>["n" H2Osat[1]; "f" 0.], # bc doesnt matter if H2O fixed
                    :HDO=>["n" HDOsat[1]; "f" 0.],
                    :O=>["f" 0.; "f" oflux],
                    :H2=>["f" 0.; "v" H2_effusion_velocity],
                    :HD=>["f" 0.; "v" HD_effusion_velocity],
                    :H=>["f" 0.; "v" H_effusion_velocity],
                    :D=>["f" 0.; "v" D_effusion_velocity],
                   );
    get(speciesbclist,
        species,
        ["f" 0.; "f" 0.])
end

# Special functions for this file ==============================================

searchdir(path, key) = filter(x->occursin(key,x), readdir(path))

function normalize(arr, base_i)
    normed_arr = arr ./ arr[base_i]
    return normed_arr
end

function get_colors(L)
    #=
    Generates some colors based on a color map for use in plotting a bunch of
    lines all at once.
    L: number of colors to generate.
    =#

    colors_dumb = [cgrad(:viridis)[x] for x in range(0, stop=1, length=L)]
    c = Array{Float64}(undef, L, 3)

    for i in range(1, length=length(colors_dumb))
        c[i, 1] = red(colors_dumb[i])
        c[i, 2] = green(colors_dumb[i])
        c[i, 3] = blue(colors_dumb[i])
    end
    return c
end

function areadensity_to_micron_atm(numpercm2)
    # #/cm^2 * (cm^2/m^2) * (10μm-am/2.687e20 #/m^2)
    return numpercm2 * (1e4/1) * (10/2.687e20)
end

# Special plotting functions (apply to all main plotting funcs) ----------------

function DH_alt_prof_plot(DHproflist, exps, v, s, optext="", optlegend="")
    #=
    DHproflist: An array with D/H profiles in each row
    exps: a list of experiment identifiers as strings
    v: water, Oflux, or temp (for putting in correct folder)
    s: a subfolder "abs/" or "mr/" for either absolute abundances or mixing ratio
    optext: an optional extension that will be appended to v, for 
            the temperature case where it goes in the temp_tradeoff_plots 
            folder but the files also specify exobase, tropopause, etc.
    optlegend: an optional string to put into the legend
    =#

    # do the DH altitudinal profile plot
    # set up plot
    fig, ax = subplots(figsize=(6,4))
    ax.set_facecolor("#ededed")
    grid(zorder=0, color="white", which="both")
    for side in ["top", "bottom", "left", "right"]
        ax.spines[side].set_visible(false)
    end
    subplots_adjust(wspace=0, bottom=0.15)
    ax.set_xlabel("D/H ratio (in atomic D, H)")
    ax.set_ylabel("Altitude (km)")
    ax.set_yticks(ticks=collect(0:50:200))
    
    # generate colors
    c = get_colors(length(exps))

    # do actual plotting
    if typeof(exps[1]) != String
        exps = [string(x) for x in exps]
    end
    for i in range(1, length=length(exps))
        ax.plot(DHproflist[i, :], alt[2:end-1]./1e5, zorder=10, color=c[i, :], 
                linewidth=3, label=optlegend*"="*exps[i])
    end

    # set savepath
    plotpath = "../Results/TradeoffPlots/"*v*"_tradeoff_plots/"*s
    savepath = plotpath*v*optext*"_DH_prof.png"
    #legend(fontsize=12, bbox_to_anchor=[1.01,1], loc=2, borderaxespad=0)
    savefig(savepath, bbox_inches="tight")

    # save it again but with log x-axis 
    xscale("log")
    xticks(rotation=45)
    savepath = plotpath*"/"*v*optext*"_DH_prof_LOG.png"
    savefig(savepath, bbox_inches="tight")
    close(fig)
end

function CO_O2_plot(xvals, ydict, xlab, pathkey, meanX, s, tempkey="")
    #=
    xvals: the variations within the experiment. list of strings
    ydict: the dictionary containing the data
    xlab: label for the x axis
    pathkey: "water", "Oflux" or "temp"
    meanX: the nominal value to plot
    s: a subfolder "abs/" for absolute abundances or "mr/" for mixing ratio
    tempkey: "_surface", "_tropopause", "_exobase"
    =#

    # Pretty ugly, but data structs are always better than if/else, right?
    paststudies = Dict("water"=>Dict("yung"=>15, "nair"=>[3,8.8]), # these are in ppm
                            "Oflux"=>Dict("yung"=>5, "nair"=>9),
                            "temp_surface"=>Dict("yung"=>220, "nair"=>214), 
                            "temp_tropopause"=>Dict("yung"=>140, "nair"=>140), 
                            "temp_exobase"=>Dict("yung"=>364, "nair"=>288))

    lookupkey = pathkey*tempkey
    ystr = s == "abs/" ? L"Abundance, CO and O$_2$" : L"Mixing ratio, CO and O$_2$"
    
    # CO, O2, CO/O2 plot
    # set up plot
    fig, ax = subplots(figsize=(6,4))
    xlabel(xlab)
    if pathkey=="water"
        xscale("log")
    end
    xticks(rotation=45)

    # CO and O2 mixin ratio axis
    ax.set_facecolor("#ededed")
    ax.plot(xvals, ydict["CO"], marker="o", zorder=10, color="cornflowerblue") 
    ax.plot(xvals, ydict["O2"], marker="o", zorder=10, color="navy")
    ax.tick_params(axis="y", labelcolor="navy", which="both")
    ax.axvline(meanX, color="black", label="Nominal value", zorder=10)
    ax.set_ylabel(ystr, color="navy")
    ax.set_yscale("log")
    for side in ["top", "bottom", "left", "right"]
        ax.spines[side].set_visible(false)
    end
    
    # CO/O2 axis
    ax2 = ax.twinx() # set up the other axis
    
    ax2.grid(zorder=ax.get_zorder()-1, color="white", axis="both")
    ax2.plot(xvals, ydict["CO/O2"], marker="o", zorder=10, color="red")
    ax2.tick_params(axis="y", labelcolor="red")
    ax2.set_ylabel("CO/O"*L"_2"*" ratio", color="red")
    ytix = collect(0:0.2:1.2)
    if maximum(ydict["CO/O2"]) / minimum(ydict["CO/O2"]) >= 10
        ax2.set_yscale("log")
    end
    for side in ["top", "bottom", "left", "right"]
        ax2.spines[side].set_visible(false)
    end

    # past studies
    if pathkey=="water"
        ax2.errorbar(paststudies[lookupkey]["nair"][1], 0.14, capsize=2, zorder=10,
                     yerr=reshape([0.05; 0.43], 2, 1), color="xkcd:brick red",
                     ecolor="xkcd:brick red", marker="s")  # Nair 1994 low water
        ax2.errorbar(paststudies[lookupkey]["nair"][2], 0.1, capsize=2, zorder=10,
                     yerr=reshape([0.07; 0.4], 2, 1), color="xkcd:fire engine red",
                     ecolor="xkcd:fire engine red", marker="s")  # Nair 1994 high water
    else
        ax2.errorbar(paststudies[lookupkey]["nair"], 0.1, capsize=2, zorder=10,
                     yerr=reshape([0.07; 0.4], 2, 1), color="xkcd:brick red",
                     ecolor="xkcd:brick red", marker="s")  # Nair 1994 high water
    end
    ax2.scatter(paststudies[lookupkey]["yung"], 0.53, marker="*", s=100,
                color="xkcd:tangerine", zorder=10)  # Yung 1988

    # more ticks for exobase surface temp
    if lookupkey == "temp_exobase"
        ax.set_xticks(ticks=vcat(xvals, [325, 350]))
        [l.set_visible(false) for (i,l) in enumerate(ax.xaxis.get_ticklabels()) if i % 2 != 0]
    end

    plotpath = "../Results/TradeoffPlots/"*pathkey*"_tradeoff_plots/"*s
    savepath = plotpath*lookupkey*"_CO_and_O2.png"
    savefig(savepath, bbox_inches="tight")
    close(fig)
end

# Main plotting functions ------------------------------------------------------
function make_water_plots(water_x, d, DHdata, q, nom_i, s)
    #=
    Makes the plots for water variation experiment
    water_x: a list of the water mixing ratios, in string format
    d: a dictionary of values of different measurables as function of
       water mixing ratio in lower atmosphere
    DH data: vertical profiles of D/H by experiment
    q: " absolute" or " mixing ratio", just used for labeling
    nom_i: index in water_x of nominal value of water mixing ratio
    s: either "abs" or "mr" 
    =#
    normed_dict = Dict()

    # make plots pretty
    rcParams = PyDict(matplotlib."rcParams")
    rcParams["font.family"] = "sans-serif"
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 22
    rcParams["axes.labelsize"]= 24
    rcParams["xtick.labelsize"] = 22
    rcParams["ytick.labelsize"] = 22

    xlab = "Water"*q*", lower atmo"

    # Do the individual plots
    plots = filter!(e->e∉["CO", "O2", "CO/O2"], [k for k in keys(d)])
    for i in plots
        # set up plot
        fig, ax = subplots(figsize=(6,4))
        ax.set_facecolor("#ededed")
        grid(zorder=0, color="white")
        xscale("log")
        for side in ["top", "bottom", "left", "right"]
            ax.spines[side].set_visible(false)
        end
        plot(water_x, d[i], marker="o", zorder=10) 
        ax.axvline(nom_i, color="black", label="Nominal value")

        # past studies
        if i=="f"
            scatter(15, 0.32, marker="*", s=100, color="xkcd:tangerine", zorder=10)  # Yung 1988
        end
        
        # set axes labels
        xlabel(L"Total atmospheric water content (pr $\mu$m)")
        ylabdict = Dict("DH"=>"D/H ratio (/1.6e-4) @ 150 km",
                        "D"=>"D"*q*" at exobase",
                        "H"=>"H"*q*" at exobase",
                        "Dflux"=>L"$\phi_D$ (cm$^{-2}$s$^{-1}$)",
                        "Hflux"=>L"$\phi_H$ (cm$^{-2}$s$^{-1}$)",
                        "f"=>"fractionation factor (f)",
                        "O2"=>"O"*L"_2"*q,
                        "HD"=>"HD"*q,
                        "H2"=>"H"*L"_2"*q,
                        "O3"=>L"O$_3$ (#/cm$^{-2}$)")
        ylabel(ylabdict[i])

        # Certain measureables have a range not suited to log scale
        nologplease = ["Dflux", "CO/O2", "DH"]  # don't logscale everything
        if ~(i in nologplease)#maximum(odict[i])/minimum(odict[i]) > 10#
            yscale("log")
            if i in ["Hflux", "HD", "H", "H2", "O3", "O2"]
                ylim(minimum(d[i])/2, maximum(d[i])*2)
            end
        end

        # set savepath
        plotpath = "../Results/TradeoffPlots/water_tradeoff_plots/" * s
        savepath = "water_"*i*".png"
        
        savefig(plotpath*savepath, bbox_inches="tight")
        close(fig)
    end
    # CO, O2, and CO/O2 plot
    CO_O2_plot(water_x, d, xlab, "water", nom_i, s)

    # DH altitude plot
    DH_alt_prof_plot(DHdata, water_x, "water", s)
end

function make_Oflux_plots(phiO, phiO_str, d, DHdata, q, nom_i, s)
    #=
    Makes the plots for O flux variation experiment
    phiO: a list of the O flux values
    phiO_str: same thing but strings, used for plotting
    d: a dictionary of values of different measurables as function of
       O flux at exobase
    DHdata: altitude profiles of D/H by experiment.
    q: " absolute" or " mixing ratio", just used for labeling
    nom_i: index in phiO of nominal value of O flux
    s: either "abs" or "mr" 
    =#
    normed_dict = Dict()

    # make plots pretty
    rcParams = PyDict(matplotlib."rcParams")
    rcParams["font.family"] = "sans-serif"
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 22
    rcParams["axes.labelsize"]= 24
    rcParams["xtick.labelsize"] = 22
    rcParams["ytick.labelsize"] = 22

    xlab = L"$\phi_O$ (cm$^{-2}$s$^{-1}$)"

    # Make individual plots showing change 
    plots = filter!(e->e∉["CO", "O2", "CO/O2"], [k for k in keys(d)])
    for i in plots
        # basic plot stuff
        fig, ax = subplots(figsize=(6,4))
        ax.set_facecolor("#ededed")
        grid(zorder=0, color="white")
        for side in ["top", "bottom", "left", "right"]
            ax.spines[side].set_visible(false)
        end
        plot(phiO_str, d[i], marker="o", zorder=10) 
        ax.axvline(nom_i, color="black", label="Nominal value")

        # past studies
        if i=="f"
            scatter(findfirst(isequal(8e7), phiO)-1, 0.32, marker="*", s=100,
                    color="green", zorder=10)  # Yung 1988
        end
        
        # set axes labels
        xlabel(L"$\phi_O$ (cm$^{-2}$s$^{-1}$)",)
        xticks(rotation=45)
        ylabdict = Dict("DH"=>"D/H ratio (/1.6e-4) @ 150 km",
                        "D"=>"D"*q*"at exobase",
                        "H"=>"H"*q*"at exobase",
                        "Dflux"=>L"$\phi_D$ (cm$^{-2}$s$^{-1}$)",
                        "Hflux"=>L"$\phi_H$ (cm$^{-2}$s$^{-1}$)",
                        "f"=>"fractionation factor (f)",
                        "O2"=>"O"*L"_2"*q,
                        "HD"=>"HD"*q,
                        "H2"=>"H"*L"_2"*q,
                        "O3"=>"O"*L"_3"*q)
        ylabel(ylabdict[i])

        # Certain measureables have a range not suited to log scale
        nologplease = ["Dflux", "CO/O2", "DH"]  # don't logscale everything
        if ~(i in nologplease)#maximum(d[i])/minimum(d[i]) > 10#
            yscale("log")
            if i in ["Hflux", "HD", "O3"]
                ylim(minimum(d[i])/2, maximum(d[i])*2)
            end
        end

        # set savepath
        plotpath = "../Results/TradeoffPlots/Oflux_tradeoff_plots/"*s
        savepath = plotpath*"O_flux_"*i*".png"
        
        savefig(savepath, bbox_inches="tight")
        close(fig)
    end
    # CO, O2 and CO/O2 Plot
    CO_O2_plot(phiO_str, d, xlab, "Oflux", nom_i, s)

    DH_alt_prof_plot(DHdata, phiO_str, "Oflux", s)
end

function make_T_plots(T, T_str, d, DHdata, exp, q, nomT, s)
    #=
    Makes the plots for temperature variation experiment
    d: a dictionary of values of different measurables as function of
       temperatures at 3 points in atmosphere: surface, tropopause, exobase
    DHdata: altitude profiles of D/H by experiment.
    T: a list of the temperature values on the x axis
    T_str: same thing but strings, used for plotting
    q: " absolute" or " mixing ratio", just used for labeling
    nomT: a dictionary of nominal values of temperature at the 3 points
    s: either "abs" or "mr" 
    =#

    # make plots pretty
    rcParams = PyDict(matplotlib."rcParams")
    rcParams["font.family"] = "sans-serif"
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 22
    rcParams["axes.labelsize"]= 24
    rcParams["xtick.labelsize"] = 22
    rcParams["ytick.labelsize"] = 22
    
    xlab = exp*" temperature"

    # loop through the parameters of interest and plot them
    plots = filter!(e->e∉["CO", "O2", "CO/O2"], [k for k in keys(d)])
    for i in plots
        # basic plot stuff
        fig, ax = subplots(figsize=(6,4))
        ax.set_facecolor("#ededed")
        grid(zorder=0, color="white")
        for side in ["top", "bottom", "left", "right"]
            ax.spines[side].set_visible(false)
        end
        plot(T[exp], d[i], marker="o", zorder=10) 
        ax.axvline(nomT[exp], color="black", label="Nominal value")
        xlabel(xlab)
        xticks(ticks=T[exp], labels=T_str[exp], rotation=45)
        ylabdict = Dict("DH"=>"D/H ratio (/1.6e-4) @ 150 km",
                        "D"=>"D"*q*"at exobase",
                        "H"=>"H"*q*"at exobase",
                        "Dflux"=>L"$\phi_D$ (cm$^{-2}$s$^{-1}$)",
                        "Hflux"=>L"$\phi_H$ (cm$^{-2}$s$^{-1}$)",
                        "f"=>"fractionation factor (f)",
                        "O2"=>"O"*L"_2"*q,
                        "HD"=>"HD"*q,
                        "H2"=>"H"*L"_2"*q,
                        "O3"=>"O"*L"_3"*q)
        ylabel(ylabdict[i])

        # past studies
        if exp=="surface"
            if i == "f"
                scatter(220, 0.32, marker="d", color="xkcd:tangerine", 
                        zorder=10)  # Yung 1988
            end
        elseif exp=="tropopause"
            if i == "f"
                scatter(125, 0.055, marker="v", color="xkcd:golden", 
                        zorder=10)  #Kras 2002 solar minimum
                scatter(125, 0.082, marker="s", color="xkcd:blood orange", 
                        zorder=10)  #Kras 2002 solar mean
                scatter(125, 0.167, marker="^", color="xkcd:scarlet", 
                        zorder=10)  #Kras 2002 solar max
                errorbar(140, 0.14, yerr=reshape([0.05; 0.43], 2, 1), 
                         fmt="s", color="purple", ecolor="gray",
                         zorder=10, capsize=2)  # Nair 1994 low water
            end
        elseif exp=="exobase"
            if i == "f"
                scatter(200, 0.055, marker="v", color="xkcd:golden", 
                        zorder=10)  #Kras 2002 solar minimum
                scatter(270, 0.082, marker="s", color="xkcd:blood orange", 
                        zorder=10)  #Kras 2002 solar mean
                scatter(350, 0.167, marker="^", color="xkcd:scarlet", 
                        zorder=10)  #Kras 2002 solar max
                errorbar(288, 0.14, yerr=reshape([0.05; 0.43], 2, 1), 
                         fmt="s", color="xkcd:hunter green", 
                         ecolor="xkcd:dark sage", zorder=10, capsize=2)  # Nair 1994 low water
                #xlim(100,400)
                xticks(ticks=vcat(T[exp], [325, 350]), 
                       labels=vcat(T_str[exp], 
                       ["325", "350"]), rotation=45)
            end
        end 

        # Certain measureables have a range not suited to log scale
        nologplease = ["Dflux", "CO/O2", "DH"]  # don't logscale everything
        if ~(i in nologplease)#maximum(d[i])/minimum(d[i]) > 10#
            yscale("log")
            if i in ["Hflux", "HD"]
                ylim(minimum(d[i])/2, maximum(d[i])*2)
            end
        end

        # set savepath
        plotpath = "../Results/TradeoffPlots/temp_tradeoff_plots/"*s
        savepath = plotpath*join(["temp", exp], "_")*"_"*i*".png"
        savefig(savepath, bbox_inches="tight")
        close(fig)
    end

    # CO/O2 plot
    CO_O2_plot(T[exp], d, exp*" temperature", "temp", nomT[exp], s, 
               "_"*exp)

    # D/H altitude plot
    DH_alt_prof_plot(DHdata, T_str[exp], "temp", s, "_"*exp, 
                     latexstring("T_{$(exp[1:3])}"))
end

function make_rel_change_plots(output_MR, output_abs, ex)
    #=
    output_MR: 
    output_abs:
    ex: experiment type, the usual key of "water", "O flux", "surface" etc.
    =#
    # make plots pretty
    rcParams = PyDict(matplotlib."rcParams")
    rcParams["font.family"] = "sans-serif"
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 22
    rcParams["axes.labelsize"]= 24
    rcParams["xtick.labelsize"] = 22
    rcParams["ytick.labelsize"] = 22
    
    # the official values of each experiment that I ran 
    xvals = Dict(#"water"=>["1e-6", "5e-6", "1e-5", "5e-5", "1e-4", "5e-4", 
    #                        "1e-3", "5e-3", "1e-2"],
                 #"water"=>["5e-6", "5e-5", "8e-4", "2.2e-3", "4.4e-3", "9e-3", "1.35e-2"]
                 "water"=>[0.3, 2, 10, 25, 50, 100, 150],  # in pr μm, total column
                 "O flux"=>["3e7", "4e7", "5e7", "6e7", "7e7", "8e7", "9e7", 
                            "1e8", "1.1e8", "1.2e8", "1.3e8", "1.4e8", "1.5e8", 
                            "1.6e8"],
                 "surface"=>[180.0, 190.0, 200.0, 210.0, 220.0, 230.0, 240.0, 
                             250.0, 260.0, 270.0],
                 "tropopause"=>[70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 
                                140.0, 150.0, 160.0],
                 "exobase"=>[125.0, 150.0, 175.0, 200.0, 225.0, 250.0, 275.0, 
                             300.0, 325.0, 350.0])
    # something to add onto the x label for each plot
    xlab_addn = Dict("water"=>" mixing ratio (ppm)", "O flux"=>L" (cm$^{-2}s^{-1}$)", 
                     "surface"=>" temperature", "tropopause"=>" temperature", 
                     "exobase"=>" temperature")


    # initialize figure and stuff for both axes
    fig, ax = subplots(3, 1, sharex=true, sharey=false, figsize=(8, 24))
    subplots_adjust(hspace=0.05)

    # nominal value plot location or index, 1-indexed (Julia)
    nom_i_julia = Dict("exobase"=>200, "tropopause"=>110, "surface"=>190, 
                       "water"=>10,
                 #"water"=>findfirst(isequal("12.6"), xvals["water"]), 
                        "O flux"=>findfirst(isequal("1.2e8"), xvals["O flux"]))
    # You have to change to 0-indexing to pass to Python, so here they are in 
    # 0-indexed format. Water has 2 subtracted because we remove one of its 
    # values later (1e-6)
    nom_i_py = Dict("exobase"=>200, "tropopause"=>110, "surface"=>190, 
                 "water"=>10,#findfirst(isequal("2.7"), xvals["water"])-2, 
                 "O flux"=>findfirst(isequal("1.2e8"), xvals["O flux"])-1)
    
    # Data from MAVEN or other missions
    shade = "gainsboro"
    for axi in ax
        axi.set_facecolor("#ededed")
        axi.grid(color="white")
        for side in ["top", "bottom", "left", "right"]
            axi.spines[side].set_visible(false)
        end
        # specific formatting by experiment type
        if ex=="exobase"
            axi.set_xticks(150:25:350)
            axi.axvspan(150, 300, color=shade)
        elseif ex=="tropopause"
            axi.axvspan(75, 160, color=shade)
        elseif ex=="surface"
            axi.axvspan(180, 270, color=shade)
        elseif ex=="O flux"
            axi.axvspan(5, 6, color=shade) # doubled to 8.6e7 from 4.3e7 per Mike
        elseif ex=="water"
            axi.axvspan(80,200,color="xkcd:brick red", alpha=0.3) #dust storm
            axi.axvspan(8, 12, color=shade) #nominal
            axi.axvspan(30, 80, color="#67a9cf", alpha=0.3) # summer pole
        end
        axi.tick_params(rotation=45, which="x")
    end

    # first ax - compare with data ---------------------------------------------
    ax[1].set_ylabel(L"(Model-Obs)/$\sigma$")
    data = [6e-4, 1.2e-3, 15, 1]  # CO MR, O2 MR, H2 (ppm), O3 (#/10^10)
    s = [1.5e-4, 0.2e-3, 5, 0.7]    # sigmas (uncertainty) on each
    
    c1 = ["#10007A", "#2F7965", "#e16262", "#D59A07"]
    # calculate the relativeness of each point wrt data
    if ex in ["water", "O flux"]
        COdiff = (output_MR["CO"] .- data[1])/s[1]
        O2diff = (output_MR["O2"] .- data[2])/s[2]
        H2diff = ((output_MR["H2"]/1e-6) .- data[3])/s[3]  # in ppm
        O3diff = (areadensity_to_micron_atm(output_abs["O3"]) .- data[4])/s[4] 
    else
        COdiff = (output_MR[ex]["CO"] .- data[1])/s[1]
        O2diff = (output_MR[ex]["O2"] .- data[2])/s[2]
        H2diff = ((output_MR[ex]["H2"]/1e-6) .- data[3])/s[3]
        O3diff = (areadensity_to_micron_atm(output_abs[ex]["O3"]) .- data[4])/s[4] 
    end

    # plot limits
    if ex in ["surface", "exobase"]#, "water"]
        toplot = 2:1:length(xvals[ex])  # remove values that make plot unreadable
    else
        toplot = 1:1:length(xvals[ex])
    end

    # plot the actual stuff            
    ax[1].plot(xvals[ex][toplot], COdiff[toplot], marker="o", color=c1[1], 
         zorder=10, label="CO")
    ax[1].plot(xvals[ex][toplot], O2diff[toplot], marker="o", color=c1[2], 
         zorder=10, label=L"O$_2$")
    ax[1].plot(xvals[ex][toplot], H2diff[toplot], marker="o", color=c1[3], 
         zorder=10, label=L"H$_2$")
    ax[1].plot(xvals[ex][toplot], O3diff[toplot], marker="o", color=c1[4], 
         zorder=10, label=L"O$_3$")
    ax[1].axvline(nom_i_py[ex], color="#444444", zorder=5)
    ax[1].axhline(0, color="black")

    # text labels
    pts1 = Dict("exobase"=>Dict("CO"=>[155, 1], "O2"=>[155, -1.4], 
                                "H2"=>[157, 9], "O3"=>[165, -3.2]),
                 "tropopause"=>Dict("CO"=>[70, 1], "O2"=>[70, -1.5], 
                                    "H2"=>[74, 10], "O3"=>[85,-3.2]), 
                 "surface"=>Dict("CO"=>[192.5, 1], "O2"=>[202, -1.1], 
                                 "H2"=>[192, 2.4], "O3"=>[210, -2.5]), 
                 "water"=>Dict("CO"=>[0.5, 7], "O2"=>[0.5, -0.5],  
                               "H2"=>[0.5, 2.7], "O3"=>[0.5, -2.4]), 
                 "O flux"=>Dict("CO"=>[3.5, -1.2], "O2"=>[0.5, 1.2], 
                                "H2"=>[1.5, -0.5], "O3"=>[1.3, -2]))
    ax[1].text(pts1[ex]["CO"][1], pts1[ex]["CO"][2], L"$f_{\mathrm{CO}}$", 
               color=c1[1])
    ax[1].text(pts1[ex]["O2"][1], pts1[ex]["O2"][2], L"$f_{\mathrm{O}_2}$", 
               color=c1[2])
    ax[1].text(pts1[ex]["H2"][1], pts1[ex]["H2"][2], L"$f_{\mathrm{H}_2}$", 
               color=c1[3])
    ax[1].text(pts1[ex]["O3"][1], pts1[ex]["O3"][2], L"[O$_3]$", color=c1[4])
    if ex=="water"
        ax[1].text(81, 4, " Dust \nstorm", color="xkcd:brick red")
        ax[1].text(25, 6.5, "Summer\n   pole", color="navy")
    end
    # Text and coordinates for the nominal value label
    nomext = Dict("exobase"=>L"T$_{exobase}$", "tropopause"=>L"T$_{tropopause}$", 
                     "surface"=>L"T$_{surface}$", "water"=>"mixing ratio", 
                     "O flux"=>L"\phi_O")
    nomcoords = Dict("exobase"=>[205, 10], "tropopause"=>[112, 10], 
                     "surface"=>[192, 3], "water"=>[2.3, 8], "O flux"=>[9.2,4])
    nomlbl = "Nominal"
    ax[1].text(nomcoords[ex][1], nomcoords[ex][2], nomlbl, color="#444444")

    # second axis: H2, HD, H, D, Hflux, Dflux ----------------------------------
    ax[2].set_ylabel(L"X/X$_{nominal}$")
    c2 = ["#6270d9","#95b034","#d14f58","#5ca85c","#ce6d30"]#["#1175c1", "#004231", "#5e0b0c", "#158913", "#d83431"]
    # here, set a dummy var to handle fact that temperature has one more index 
    # for lookup. 
    if ex in ["water", "O flux"]
        dictvar = output_MR
    else
        dictvar = output_MR[ex]
    end

    # to do the division/normalization we need to re-find the index of the value
    # against which to normalize. To do this, look in the xvals by experiment, 
    # then index that according to whether we cut it or not. There doesn't seem 
    # to be a cleaner way to do the indexing of cut.
    normidx = Dict("exobase"=>findfirst(isequal(200), xvals["exobase"][2:1:length(xvals["exobase"])]), 
                   "tropopause"=>findfirst(isequal(110), xvals["tropopause"][1:1:length(xvals["tropopause"])]), 
                   "surface"=>findfirst(isequal(190), xvals["surface"][2:1:length(xvals["surface"])]), 
                   "water"=>findfirst(isequal(10), xvals["water"][1:1:length(xvals["water"])]),
                   "O flux"=>findfirst(isequal("1.2e8"), xvals["O flux"][1:1:length(xvals["O flux"])]))

    HDdiff = normalize(dictvar["HD"][toplot], normidx[ex])
    Hdiff = normalize(dictvar["H"][toplot], normidx[ex])
    Ddiff = normalize(dictvar["D"][toplot], normidx[ex])
    Hfdiff = normalize(dictvar["Hflux"][toplot], normidx[ex])
    Dfdiff = normalize(dictvar["Dflux"][toplot], normidx[ex])

    # plot the actual stuff
    ax[2].plot(xvals[ex][toplot], HDdiff, marker="o", color=c2[1], zorder=10, 
               label="HD")
    ax[2].plot(xvals[ex][toplot], Hdiff, marker="o", color=c2[2], zorder=10, 
               label="H")
    ax[2].plot(xvals[ex][toplot], Ddiff, marker="o", color=c2[3], zorder=10,
               label="D")
    ax[2].plot(xvals[ex][toplot], Hfdiff, marker="o", color=c2[4], zorder=10, 
               label=L"\phi_H")
    ax[2].plot(xvals[ex][toplot], Dfdiff, marker="o", color=c2[5], zorder=10, 
               label=L"\phi_D")
    ax[2].axvline(nom_i_py[ex], color="#444444", zorder=5)
    
    if ex=="exobase"
        ax[2].set_yscale("log")
    end

    # text labels
    pts2 = Dict("exobase"=>Dict("HD"=>[150, 3], "H"=>[150, 1.6], "D"=>[165, 0.6], 
                                "Hflux"=>[150, 0.5], "Dflux"=>[280, 50]),
                "tropopause"=>Dict("HD"=>[70, 2.55], "H"=>[70, 1.4], "D"=>[140, 1.2],
                                   "Hflux"=>[150, 1.2], "Dflux"=>[150, 6]), 
                "surface"=>Dict("HD"=>[210, 1.4], "H"=>[220, 0.4], "D"=>[193, 0.8],
                                "Hflux"=>[215, 1.1], "Dflux"=>[200, 1.7]), 
                "water"=>Dict("HD"=>[0.7, 1.15], "H"=>[0.3, 1.01], "D"=>[0.5, 1.165], 
                              "Hflux"=>[1, 1.01], "Dflux"=>[0.4, 1.18]), 
                "O flux"=>Dict("HD"=>[7, 0.75], "H"=>[0.2, 0.55], "D"=>[0, 0.82],
                               "Hflux"=>[2, 0.5], "Dflux"=>[3, 0.4])
                )
    ax[2].text(pts2[ex]["HD"][1], pts2[ex]["HD"][2], "HD", color=c2[1])
    ax[2].text(pts2[ex]["H"][1], pts2[ex]["H"][2], "H", color=c2[2])
    ax[2].text(pts2[ex]["D"][1], pts2[ex]["D"][2], "D", color=c2[3])
    ax[2].text(pts2[ex]["Hflux"][1], pts2[ex]["Hflux"][2], L"$\phi_H$", color=c2[4])
    ax[2].text(pts2[ex]["Dflux"][1], pts2[ex]["Dflux"][2], L"$\phi_D$", color=c2[5])    


    # third axis: f ------------------------------------------------------------
    C = "black"

    ax[3].set_xlabel(ex*xlab_addn[ex])
    ax[3].plot(xvals[ex][toplot], dictvar["f"][toplot], marker="o", color="xkcd:royal purple", 
               zorder=10)
    ax[3].axvline(nom_i_py[ex], color="#444444", zorder=5)

    ylims = Dict("tropopause"=>[10.0^(-4), 10.0^(-2.0)], "exobase"=>[10.0^(-5.0), 10.0^(-1.0)])

    if ex in ["exobase", "tropopause"]  # in these cases, we'll have 2 axes.
        # some adjustments to the main axes, which shows raw value
        ax[3].set_yscale("log") 
        ax[3].tick_params(axis="y", labelcolor="xkcd:royal purple", which="both")
        ax[3].set_ylim(ylims[ex][1], ylims[ex][2])

        # plot the secondary axis showing increase
        fdiff = normalize(dictvar["f"][toplot], normidx[ex])
        ax3_2 = ax[3].twinx()
        ax3_2.plot(xvals[ex][toplot], fdiff, marker="o", color="xkcd:forest")
        ax3_2.set_ylabel(L"$f$/$f_{nominal}$", color="xkcd:forest")
        ax3_2.tick_params(axis="y", labelcolor="xkcd:forest", which="both")
        for side in ["top", "bottom", "left", "right"]
            ax3_2.spines[side].set_visible(false)
        end
        C = "xkcd:royal purple"
    end
    ax[3].set_ylabel(L"$f$ (raw value)", color=C) # this has to be here because the color changes
    if ex=="water"
        xscale("log")
    end
    pathbase = "/home/emc/GDrive-CU/Research/Results/TradeoffPlots/"
    savefig(pathbase*"metrics_tradeoff_"*ex*".png", bbox_inches="tight")
end


# Analyzation functions (main routines) ----------------------------------------

function analyze_water(abs_or_mr, make_plots=false, path=lead)
    # Establish parameters, filenames, etc
    #watervals_str = ["1e-6", "5e-6", "1e-5", "5e-5", "1e-4", "5e-4", "1e-3", "5e-3", "1e-2"]
    watervals = [0.3, 2, 10, 25, 50, 100, 150]
    watervals_str = ["5e-6", "5e-5", "8e-4", "2.2e-3", "4.4e-3", "9e-3", "1.35e-2"]
    wfilelist = [path*"water_"*w*"/converged_standardwater_D_water_"*w*".h5" for w in watervals_str]
    temps = [192.0, 110.0, 199.0]
    oflux = 1.2e8
    q = abs_or_mr == "abs" ? " abundance" : " mixing ratio" # for labels
    mean_idx = findfirst(isequal(10), watervals) - 1 
    subfolder = abs_or_mr == "abs" ? "abs/" : "mr/"
    
    # Establish variables to store data on simulations
    # Easier to deal with D/H profiles separately due to different array size
    DHprofs = Array{Any}(undef, length(watervals_str), length(alt)-2)  
    wdict = Dict{String, Array}("O2"=>[], "HD"=>[], "H2"=>[], "H"=>[], "D"=>[], "CO"=>[], 
                 "Hflux"=>[], "Dflux"=>[], "f"=>[], "CO/O2"=>[], "O3"=>[], 
                 "DH"=>[])

    # loop water files and collect data 
    i = 1
    for wfile in wfilelist
        # get the current array
        ncur = get_ncurrent(wfile)

        N0 = abs_or_mr == "abs" ? 1 : n_tot(ncur, 0)
        Ntop = abs_or_mr == "abs" ? 1 : n_tot(ncur, 200e5)
        LA = collect(2e5:2e5:80e5)

        # Calculate the things we care about
        # O2 Mixing ratio at surface
        append!(wdict["O2"], ncur[:O2][1]/N0)
        # HD mixing ratio
        append!(wdict["HD"], ncur[:HD][1]/N0)
        # H2 mixing ratio
        append!(wdict["H2"], sum(ncur[:H2][1:length(LA)])/sum([n_tot(ncur, h) for h in LA]))
        # D Mixing ratio
        append!(wdict["D"], ncur[:D][end]/Ntop)
        # H mixing ratio
        append!(wdict["H"], ncur[:H][end]/Ntop)
        # D/H at 150 km
        append!(wdict["DH"], (ncur[:D][75]/ncur[:H][75])/(1.6e-4))
        # CO mixing ratio
        append!(wdict["CO"], ncur[:CO][1]/N0)
        # CO/O2 ratio
        append!(wdict["CO/O2"], ncur[:CO][1]/ncur[:O2][1])
        # O3 mixing ratio
        append!(wdict["O3"], sum(ncur[:O3])*2e5) # gets O3 in #/cm^2
        # H and D fluxes
        Hf = get_H_fluxes(wfile, oflux, temps)
        Df = get_D_fluxes(wfile, oflux, temps)
        append!(wdict["Hflux"], Hf)
        append!(wdict["Dflux"], Df)
        # fractionation factor 
        append!(wdict["f"], 2*(Df/Hf) / (ncur[:HDO][1]/ncur[:H2O][1]))
        # D/H profile
        DHprofs[i, :] = ncur[:D] ./ ncur[:H]  # altitude profile
        i += 1
    end

    if make_plots == true
        make_water_plots(watervals, wdict, DHprofs, q, mean_idx, subfolder)
    end

    println("Finished water plots")
    return wdict
end

function analyze_Oflux(abs_or_mr, make_plots=false, path=lead)
    # Establish important parameters, files, etc
    Ofluxvals = [3e7, 4e7, 5e7, 6e7, 7e7, 8e7, 9e7, 1e8, 1.1e8, 1.2e8, 1.3e8, 
                 1.4e8, 1.5e8, 1.6e8]
    Ofluxvals_str = ["3e7", "4e7", "5e7", "6e7", "7e7", "8e7", "9e7", "1e8", 
                     "1.1e8", "1.2e8", "1.3e8", "1.4e8", "1.5e8", "1.6e8"]
    Ofilelist = [path*"Oflux_"*o*"/converged_standardwater_D_Oflux_"*o*".h5" 
                 for o in Ofluxvals_str]
    temps = [192.0, 110.0, 199.0]
    mean_idx = findfirst(isequal(1.2e8), Ofluxvals) - 1
    q = abs_or_mr == "abs" ? " abundance " : " mixing ratio " # set label
    subfolder = abs_or_mr == "abs" ? "abs/" : "mr/"

    # Establish variables to store data on simulations
    odict = Dict("O2"=>[], "HD"=>[], "H2"=>[], "H"=>[], "D"=>[], "CO"=>[], 
                 "Hflux"=>[], "Dflux"=>[], "f"=>[], "CO/O2"=>[], "O3"=>[], 
                 "DH"=>[])

    # Easier to deal with D/H profiles separately due to different array size
    DHprofs = Array{Any}(undef, length(Ofluxvals_str), length(alt)-2)
    #HDprofs = Array{Any}(undef, length(Ofluxvals_str), length(alt)-2)

    # loop through O flux files 
    i = 1
    for (oflux, ofile) in zip(Ofluxvals, Ofilelist)
        # calculate the things we care about
        # get the current array
        ncur = get_ncurrent(ofile)

        N0 = abs_or_mr == "abs" ? 1 : n_tot(ncur, 0)
        Ntop = abs_or_mr == "abs" ? 1 : n_tot(ncur, 200e5)
        LA = collect(2e5:2e5:80e5)

        # Calculate the things we care about
        # O2 Mixing ratio at surface
        append!(odict["O2"], ncur[:O2][1]/N0)
        # HD mixing ratio
        append!(odict["HD"], ncur[:HD][1]/N0)
        # H2 mixing ratio
        append!(odict["H2"], sum(ncur[:H2][1:length(LA)])/sum([n_tot(ncur, h) for h in LA]))
        # D Mixing ratio
        append!(odict["D"], ncur[:D][end]/Ntop)
        # H mixing ratio
        append!(odict["H"], ncur[:H][end]/Ntop)
        # D/H at 150 km
        append!(odict["DH"], (ncur[:D][75]/ncur[:H][75])/(1.6e-4))
        # CO mixing ratio
        append!(odict["CO"], ncur[:CO][1]/N0)
        # CO/O2 ratio
        append!(odict["CO/O2"], ncur[:CO][1]/ncur[:O2][1])
        # O3 mixing ratio
        append!(odict["O3"], sum(ncur[:O3])*2e5) # gets O3 in #/cm^2
        # H and D fluxes
        Hf = get_H_fluxes(ofile, oflux, temps)
        Df = get_D_fluxes(ofile, oflux, temps)
        append!(odict["Hflux"], Hf)
        append!(odict["Dflux"], Df)
        # fractionation factor 
        append!(odict["f"], 2*(Df/Hf) / (ncur[:HDO][1]/ncur[:H2O][1]))
        # D/H profile
        DHprofs[i, :] = ncur[:D] ./ ncur[:H]  # altitude profiles
        #HDprofs[i, :] = ncur[:HD] ./ [n_tot(ncur, i) for i in alt[2:end-1]]
        i += 1
    end

    if make_plots == true
        make_Oflux_plots(Ofluxvals, Ofluxvals_str, odict, DHprofs, q, mean_idx, 
                         subfolder)
    end

    println("Finished O flux")
    return odict
end
    
function analyze_T(abs_or_mr, make_plots=false, path=lead)
    # Set up parameters, filenames, etc
    tvals = Dict("surface"=>[180.0, 190.0, 200.0, 210.0, 
                             220.0, 230.0, 240.0, 250.0, 260.0, 270.0],
                 "tropopause"=>[70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 
                                140.0, 150.0, 160.0],
                 "exobase"=>[125.0, 150.0, 175.0, 200.0, 225.0, 250.0, 
                             275.0, 300.0, 325.0, 350.0])
    # Dictionaries to store each experiment's data so it can be returned
    all_tdicts = Dict()

    tvals_str = Dict()
    tempfilelist = Dict()
    for k in keys(tvals)
        tvals_str[k] = [string(trunc(Int, x)) for x in tvals[k]]
        if k == "surface"
            tempfilelist[k] = [path*"temp_"*t*"_110_199"*"/converged_standardwater_D_temp_"*t*"_110_199.h5" for t in tvals_str[k]]
        elseif k == "tropopause"
            tempfilelist[k] = [path*"temp_192_"*t*"_199"*"/converged_standardwater_D_temp_192_"*t*"_199.h5" for t in tvals_str[k]]
        elseif k == "exobase"
            tempfilelist[k] = [path*"temp_192_110_"*t*"/converged_standardwater_D_temp_192_110_"*t*".h5" for t in tvals_str[k]]
        end
    end
    meanT = Dict("surface"=>190, "tropopause"=>110, "exobase"=>200)  # nominal 
    oflux = 1.2e8
    q = abs_or_mr == "abs" ? " abundance " : " mixing ratio " # set label
    subfolder = abs_or_mr == "abs" ? "abs/" : "mr/"

    # loop through which temp is varied and construct a list of datapoints
    for experiment in keys(tvals) # loop across the dictionary
        println("Now doing temperature ", experiment)
        tdict = Dict("O2"=>[], "HD"=>[], "H2"=>[], "H"=>[], "D"=>[], "CO"=>[], 
                     "Hflux"=>[], "Dflux"=>[], "f"=>[], "CO/O2"=>[], "O3"=>[], 
                     "DH"=>[])
        
        DHprofs = Array{Any}(undef, length(tvals_str[experiment]), length(alt)-2)
        #HDprofs = Array{Any}(undef, length(tvals_str[experiment]), length(alt)-2)
        

        # now loop through the values for each varied temp
        i = 1
        for (tv, tfile) in zip(tvals[experiment], tempfilelist[experiment])
            # set the temperature profile
            if experiment == "surface" 
                temps = [tv, 110.0, 199.0]
            elseif experiment == "tropopause"
                temps = [192.0, tv, 199.0]
            elseif experiment == "exobase"
                temps = [192.0, 110.0, tv]
            end

            # get the current array
            ncur = get_ncurrent(tfile)
            N0 = abs_or_mr == "abs" ? 1 : n_tot(ncur, 0)
            Ntop = abs_or_mr == "abs" ? 1 : n_tot(ncur, 200e5)
            LA = collect(2e5:2e5:80e5)

            # Calculate the things we care about
            # O2 Mixing ratio at surface
            append!(tdict["O2"], ncur[:O2][1]/N0)
            # HD mixing ratio
            append!(tdict["HD"], ncur[:HD][1]/N0)
            # H2 mixing ratio
            append!(tdict["H2"], sum(ncur[:H2][1:length(LA)])/sum([n_tot(ncur, h) for h in LA]))
            # D Mixing ratio
            append!(tdict["D"], ncur[:D][end]/Ntop)
            # H mixing ratio
            append!(tdict["H"], ncur[:H][end]/Ntop)
            # D/H at 150 km
            append!(tdict["DH"], (ncur[:D][75]/ncur[:H][75])/(1.6e-4))
            # CO mixing ratio
            append!(tdict["CO"], ncur[:CO][1]/N0)
            # CO/O2 ratio
            append!(tdict["CO/O2"], ncur[:CO][1]/ncur[:O2][1])
            # O3 mixing ratio
            append!(tdict["O3"], sum(ncur[:O3])*2e5) # gets O3 in #/cm^2
            # H and D fluxes
            Hf = get_H_fluxes(tfile, oflux, temps)
            Df = get_D_fluxes(tfile, oflux, temps)
            append!(tdict["Hflux"], Hf)
            append!(tdict["Dflux"], Df)
            # fractionation factor 
            append!(tdict["f"], 2*(Df/Hf) / (ncur[:HDO][1]/ncur[:H2O][1]))
            # D/H profile
            DHprofs[i, :] = ncur[:D] ./ ncur[:H]  # altitude profile
            #HDprofs[i, :] = ncur[:HD] ./ [n_tot(ncur, i) for i in alt[2:end-1]]
            i += 1
        end

        if make_plots == true
            make_T_plots(tvals, tvals_str, tdict, DHprofs, experiment, q, meanT, 
                         subfolder)
        end
        
        println("Finished temperature ", experiment)

        all_tdicts[experiment] = tdict
    end
    println("Finished temps")
    return all_tdicts
end

# Function calls ===============================================================

makeplots = false
write_new_files = false

water_data_mr = analyze_water("mr", makeplots)
water_data_abs = analyze_water("abs", makeplots)

println()

o_data_abs = analyze_Oflux("abs", makeplots)
o_data_mr = analyze_Oflux("mr", makeplots)

println()

T_data_mr = analyze_T("mr", makeplots)
T_data_abs = analyze_T("abs", makeplots)

if write_new_files
    wd_mr = jldopen("water_MR_data.jld", "w")
    @write wd_mr water_data_mr
    close(wd_mr)
    wd_abs = jldopen("water_abs_data.jld", "w")
    @write wd_abs water_data_abs
    close(wd_abs)

    O_mr = jldopen("O_MR_data.jld", "w")
    @write O_mr o_data_mr
    close(O_mr)
    O_abs = jldopen("O_abs_data.jld", "w")
    @write O_abs o_data_abs
    close(O_abs)

    T_mr = jldopen("T_MR_data.jld", "w")
    @write T_mr T_data_mr
    close(T_mr)
    T_abs = jldopen("T_abs_data.jld", "w")
    @write T_abs T_data_abs
    close(T_abs)
end


# Relative changes plot with two panels ========================================
# panel 1: metrics CO, O2, O3, and H2
# panel 2: all the other measureablesbles except fractionation factor
# panel 3: frationatoin factor

make_rel_change_plots(water_data_mr, water_data_abs, "water")
make_rel_change_plots(o_data_mr, o_data_abs, "O flux")
make_rel_change_plots(T_data_mr, T_data_abs, "surface")
make_rel_change_plots(T_data_mr, T_data_abs, "tropopause")
make_rel_change_plots(T_data_mr, T_data_abs, "exobase")