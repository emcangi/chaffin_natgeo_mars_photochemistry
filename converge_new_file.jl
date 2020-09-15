###############################################################################
# converge_new_file.jl
# TYPE: (1) Model files - required
# DESCRIPTION: does initial convergence for a photochemistry experiment of the
# Martian atmosphere.
# 
# Eryn Cangi
# Created 2018
# Last edited: 21 July 2020
# Currently tested for Julia: 1.4.1
###############################################################################

using PyPlot
using PyCall
using HDF5, JLD
using LaTeXStrings
using Distributed
using DelimitedFiles
using SparseArrays
using LinearAlgebra
using ProgressMeter
using Photochemistry  # custom module for this project
# The following functions are imported directly in order to overload them, making
# multiple methods available. 
import Photochemistry.fluxcoefs, Photochemistry.Dcoef, Photochemistry.scaleH


include("PARAMETERS.jl")

################################################################################
#                                   FUNCTIONS                                  #
################################################################################

#=
These functions are required to be in this file for one of three reasons:
1) because they, or one of their overrides, call on the dynamically-defined 
   Temp(z) function,
2) Because they call a function that is couched in an @eval statement, which 
   cannot be relocated,
3) They are the main routine functions and will not be shared by any other scripts.

The Temp(z) function is defined during each run as a shortcut to Tpiecewise(), 
which takes as arguments the temperatures at the surface, tropopause, and exobase.
Temp(z) is defined so that those arguments need not be passed constantly. This 
may be fixed in the future.

For additional functions, see the Photochemistry module.
=#

# transport/scale height =======================================================

function scaleH(z, species::Symbol)
    #=
    Same as first scaleH, but for a particular atomic/molecular species.
    =#
    T = Temp(z)
    mm = speciesmolmasslist[species]
    return scaleH(z, T, mm)
end

# transport functions ==========================================================

#=
    at each level of the atmosphere, density can be transferred up or
    down with a specified rate coefficient.

                             | n_i+1
          ^  tspecies_i_up   v tspecies_i+1_down
      n_i |
          v tspecies_i_down  ^ tspecies_i-1_up
                             | n_i-1

    the flux at each cell boundary is the sum of the upward flux from
    the cell below and the downward flux of the cell above. These fluxes
    are determined using flux coefficients that come from the diffusion
    equation. Care must be taken at the upper and lower boundary so that
    tspecies_top_up and tspecies_bottom_down properly reflect the
    boundary conditions of the atmosphere.

    This is handled in the code with the population of the appropriate
    reactions, with variable rate coeffecients that are populated
    between each timestep (similar to the way photolysis rates are
    included). We need to make reactions at each interior altitude
    level:
             n_i -> n_i+1  tspecies_i_up
             n_i -> n_i-1  tspecies_i_down

    At the upper and lower boundary we omit the species on the RHS, so
    that these reactions are potentially non-conservative:

            n_top    -> NULL  tspecies_top_up
            n_bottom -> NULL  tspecies_bottom_down

    These coefficients then describe the diffusion velocity at the top
    and bottom of the atmosphere.
=#

function fluxcoefs(z, dz, species, n_current)
    #= 
    overload to generate the coefficients K, D, T, Hs if they are not supplied (most common)

    z: altitude in cm
    dz: altitudinal layer thickness in cm
    species: species symbol for which flux coefficients are calculated
    n_current: current atmospheric state array

    p: upper layer ("plus")
    0: this layer
    m: lower layer ("minus")
    =#

    ntp = n_tot(n_current, z+dz)
    nt0 = n_tot(n_current, z)
    ntm = n_tot(n_current, z-dz)
    Kp = Keddy(z+dz, ntp)
    K0 = Keddy(z, nt0)
    Km = Keddy(z-dz, ntm)
    Tp = Temp(z+dz)
    T0 = Temp(z)
    Tm = Temp(z-dz)
    Dp = Dcoef(Tp, ntp, species)
    D0 = Dcoef(T0, nt0, species)
    Dm = Dcoef(Tm, ntm, species)
    Hsp = scaleH(z+dz, species)
    Hs0 = scaleH(z, species)
    Hsm = scaleH(z-dz, species)
    H0p = scaleH(z+dz, Tp, n_current)
    H00 = scaleH(z, T0, n_current)
    H0m = scaleH(z-dz, Tm, n_current)

    # return the coefficients
    return fluxcoefs(z, dz, [Km , K0, Kp], [Dm , D0, Dp], [Tm , T0, Tp],
                     [Hsm, Hs0, Hsp], [H0m, H00, H0p], species)
end

function lower_up(z, dz, species, n_current)
    #= 
    define transport coefficients for a given atmospheric layer for
    transport from that layer to the one above. 
    p: layer above ("plus"), 0: layer at altitude z, m: layer below ("minus") 

    z: altitude in cm
    dz: altitude layer thickness ("resolution"), in cm
    species: Symbol; species for which this coefficients are calculated
    n_current: Array; species number density by altitude

    returns: return of fluxcoefs
    =#
    ntp = n_tot(n_current, z+dz)
    nt0 = n_tot(n_current, z)
    ntm = 1
    Kp = Keddy(z+dz, ntp)
    K0 = Keddy(z,nt0)
    Km = 1
    Tp = Temp(z+dz)
    T0 = Temp(z)
    Tm = 1
    Dp = Dcoef(Tp, ntp, species)
    D0 = Dcoef(T0, nt0, species)
    Dm = 1
    Hsp = scaleH(z+dz, species)
    Hs0 = scaleH(z,species)
    Hsm = 1
    H0p = scaleH(z+dz, Tp, n_current)
    H00 = scaleH(z,T0, n_current)
    H0m = 1

    # return the coefficients
    return fluxcoefs(z, dz,
              [Km , K0, Kp],
              [Dm , D0, Dp],
              [Tm , T0, Tp],
              [Hsm, Hs0, Hsp],
              [H0m, H00, H0p],
              species)[2]
end

function upper_down(z, dz, species, n_current)
    #= 
    define transport coefficients for a given atmospheric layer for
    transport from that layer to the one below. 
    p: layer above ("plus"), 0: layer at altitude z, m: layer below ("minus") 

    z: altitude in cm
    dz: altitude layer thickness ("resolution"), in cm
    species: Symbol; species for which this coefficients are calculated
    n_current: Array; species number density by altitude

    returns: return of fluxcoefs
    =#
    ntp = 1
    nt0 = n_tot(n_current, z)
    ntm = n_tot(n_current, z-dz)
    Kp = 1
    K0 = Keddy(z, nt0)
    Km = Keddy(z-dz, ntm)
    Tp = 1
    T0 = Temp(z)
    Tm = Temp(z-dz)
    Dp = 1
    D0 = Dcoef(T0, nt0, species)
    Dm = Dcoef(Tm, ntm, species)
    Hsp = 1
    Hs0 = scaleH(z, species)
    Hsm = scaleH(z-dz, species)
    H0p = 1
    H00 = scaleH(z, T0, n_current)
    H0m = scaleH(z-dz, Tm, n_current)

    # return the coefficients
    return fluxcoefs(z, dz,
              [Km , K0, Kp],
              [Dm , D0, Dp],
              [Tm , T0, Tp],
              [Hsm, Hs0, Hsp],
              [H0m, H00, H0p],
              species)[1]
end

function boundaryconditions(species, speciesbclist, dz, n_current)
    #= 
    returns the symbolic transport coefficients that encode the
    boundary conditions for the null-pointing equations

    n_1->NULL t_lower_bc
    n_f->NULL t_upper_bc

    this defines two additional symbols for each species that need
    to be resolved in the function call macro:
                tspecies_lower_up and
                tspecies_upper_down
    these are found by passing the appropriate values to fluxcoefs
    and selecting the correct output.

    species: Symbol
    speciesbclist: Has to be passed in to pass on to speciesbcs. Because water
                   saturation at the surface can vary, but is a bc.
    dz: Float64; layer thickness in cm
    n_current: Array; species number density by altitude

    returns: 2x2 boundary condition array where the first row is for the surface
             layer and second row is for the top of the atmosphere. 
    =#

    bcs = speciesbcs(species, speciesbclist)
    if issubset([species],notransportspecies)
        bcs = ["f" 0.; "f" 0.]
    end

    # first element returned corresponds to lower BC, second to upper
    # BC transport rate. Within each element, the two rates correspond
    # to the two equations
    # n_b  -> NULL (first rate, depends on species concentration)
    # NULL -> n_b  (second rate, independent of species concentration 
    bcvec = Float64[0 0;0 0]

    # LOWER
    if bcs[1, 1] == "n"
        bcvec[1,:]=[fluxcoefs(alt[2], dz, species, n_current)[1],
                    lower_up(alt[1], dz, species, n_current)*bcs[1, 2]]
    elseif bcs[1, 1] == "f"
        bcvec[1,:] = [0.0, bcs[1, 2]/dz]
    elseif bcs[1, 1] == "v"
        bcvec[1,:] = [bcs[1, 2]/dz, 0.0]
    else
        throw("Improper lower boundary condition!")
    end

    # UPPER
    if bcs[2, 1] == "n"
        bcvec[2,:] = [fluxcoefs(alt[end-1],dz, species, n_current)[2],
                    upper_down(alt[end],dz, species, n_current)*bcs[1, 2]]
    elseif bcs[2, 1] == "f"
            bcvec[2,:] = [0.0,-bcs[2, 2]/dz]
    elseif bcs[2, 1] == "v"
        bcvec[2,:] = [bcs[2, 2]/dz, 0.0]
    else
        throw("Improper upper boundary condition!")
    end

    return bcvec
end

# chemistry functions ==========================================================

function ratefn(nthis, inactive, inactivespecies, activespecies, Jrates, T, M, tup, tdown, tlower, tupper)
    # at each altitude, get the appropriate group of concentrations,
    # coefficients, and rates to pass to ratefn_local
    nthismat = reshape(nthis, (length(activespecies), length(intaltgrid)))
    inactivemat = reshape(inactive,(length(inactivespecies),length(intaltgrid)))

    returnrates = zero(nthismat)

    # fill the first altitude entry with information for all species
    returnrates[:,1] = ratefn_local([nthismat[:,1]; nthismat[:,2];
                                    fill(1.0, length(activespecies));
                                    inactivemat[:,1]; Jrates[:,1]; T[1]; M[1];
                                    tup[:,1]; tlower[:,1]; tdown[:,2];
                                    tlower[:,2]]...)

    # iterate through other altitudes except the last level, filling the info in
    for ialt in 2:(length(intaltgrid)-1)
        returnrates[:,ialt] = ratefn_local([nthismat[:,ialt];
                                          nthismat[:,ialt+1];
                                          nthismat[:,ialt-1];
                                          inactivemat[:,ialt];
                                          Jrates[:,ialt];
                                          T[ialt]; M[ialt];
                                          tup[:,ialt]; tdown[:,ialt];
                                          tdown[:,ialt+1]; tup[:,ialt-1]]...)
    end

    # fill in the last level of altitude (200 km)
    returnrates[:,end] = ratefn_local([nthismat[:,end];
                                       fill(1.0, length(activespecies));
                                       nthismat[:,end-1];
                                       inactivemat[:,end];
                                       Jrates[:,end];
                                       T[end]; M[end];
                                       tupper[:,1]; tdown[:,end];
                                       tupper[:,2]; tup[:,end-1]]...)
    return [returnrates...;]
end

function chemJmat(nthis, inactive, activespecies, inactivespecies, Jrates, T, M, tup, tdown, tlower, tupper, dt)
    nthismat = reshape(nthis, (length(activespecies), length(intaltgrid)))
    inactivemat = reshape(inactive, (length(inactivespecies), length(intaltgrid)))
    chemJi = Int64[]
    chemJj = Int64[]
    chemJval = Float64[]

    (tclocal, tcupper, tclower) = chemJmat_local([nthismat[:,1]; nthismat[:,2];
                                                  fill(1.0, length(activespecies));
                                                  inactivemat[:,1]; Jrates[:,1];
                                                  T[1]; M[1]; tup[:,1]; tlower[:,1];
                                                  tdown[:,2]; tlower[:,2];dt]...)
    #add the influence of the local densities
    append!(chemJi, tclocal[1])
    append!(chemJj, tclocal[2])
    append!(chemJval, tclocal[3])
    #and the upper densities
    append!(chemJi, tcupper[1])
    append!(chemJj, tcupper[2] .+ length(activespecies))
    append!(chemJval, tcupper[3])

    for ialt in 2:(length(intaltgrid)-1)
        (tclocal, tcupper, tclower) = chemJmat_local([nthismat[:,ialt];
                                                      nthismat[:,ialt+1];
                                                      nthismat[:,ialt-1];
                                                      inactivemat[:,ialt];
                                                      Jrates[:,ialt]; T[ialt]; 
                                                      M[ialt]; tup[:,ialt]; 
                                                      tdown[:,ialt];
                                                      tdown[:,ialt+1];
                                                      tup[:,ialt-1]; dt]...)
        # add the influence of the local densities
        append!(chemJi, tclocal[1].+(ialt-1)*length(activespecies))
        append!(chemJj, tclocal[2].+(ialt-1)*length(activespecies))
        append!(chemJval, tclocal[3])
        # and the upper densities
        append!(chemJi, tcupper[1].+(ialt-1)*length(activespecies))
        append!(chemJj, tcupper[2].+(ialt  )*length(activespecies))
        append!(chemJval, tcupper[3])
        # and the lower densities
        append!(chemJi, tclower[1].+(ialt-1)*length(activespecies))
        append!(chemJj, tclower[2].+(ialt-2)*length(activespecies))
        append!(chemJval, tclower[3])
    end

    (tclocal, tcupper, tclower) = chemJmat_local([nthismat[:,end];
                                              fill(1.0, length(activespecies));
                                              nthismat[:,end-1];
                                              inactivemat[:,end];
                                              Jrates[:,end];
                                              T[end]; M[end];
                                              tupper[:,1]; tdown[:,end];
                                              tupper[:,2]; tup[:,end-1];dt]...)
    # add the influence of the local densities
    append!(chemJi, tclocal[1].+(length(intaltgrid)-1)*length(activespecies))
    append!(chemJj, tclocal[2].+(length(intaltgrid)-1)*length(activespecies))
    append!(chemJval, tclocal[3])
    # and the lower densities
    append!(chemJi, tclower[1].+(length(intaltgrid)-1)*length(activespecies))
    append!(chemJj, tclower[2].+(length(intaltgrid)-2)*length(activespecies))
    append!(chemJval, tclower[3])

    # make sure to add 1's along the diagonal
    append!(chemJi,[1:length(nthis);])
    append!(chemJj,[1:length(nthis);])
    append!(chemJval, fill(1.0, length(nthis)))

    sparse(chemJi, chemJj, chemJval, length(nthis), length(nthis), +);
end

# diffusion functions ==========================================================

# Overload
Dcoef(z, species::Symbol, n_current) = Dcoef(Temp(z), n_tot(n_current, z), species)

# main routine functions =======================================================
function update_Jrates!(n_current::Dict{Symbol, Array{Float64, 1}})
    #=
    this function updates the photolysis rates stored in n_current to
    reflect the altitude distribution of absorbing species
    =#
    #    global solarabs::Array{Array{Float64, 1},1}

    # Initialize an array, length=num of altitude levels - 2.
    # Each sub-array is an array of length 2000, corresponding to 2000 wavelengths.
    solarabs = Array{Array{Float64}}(undef, length(alt)-2)
    for i in range(1, length=length(alt)-2)
        solarabs[i] = zeros(Float64, 2000)
    end

    nalt = size(solarabs, 1)
    nlambda = size(solarabs[1],1)

    for jspecies in Jratelist
        species = absorber[jspecies]

        jcolumn = 0.
        for ialt in [nalt:-1:1;]
            #get the vertical column of the absorbing constituient
            jcolumn += n_current[species][ialt]*dz
            # if jspecies==:JO2toOpO
            #     println(string("At alt = ",alt[ialt+1],
            #                    ", n_",species," = ",n_current[species][ialt],
            #                    ", jcolumn = ",jcolumn))
            #     println("and solarabs[ialt] is $(solarabs[ialt]) before we do axpy")
            #     readline(STDIN)
            # end

            # add the total extinction to solarabs:
            # multiplies air column density at all wavelengths by crosssection
            # to get optical depth. This is an override of axpy! to use the 
            # full arguments. For the equation Y' = alpha*X + Y:
            # ARG 1: n (length of arrays in ARGS 3, 5)
            # ARG 2: alpha, a scalar.
            # ARG 3: X, an array of length n.
            # ARG 4: the increment of the index values of X, maybe?
            # ARG 5: Y, an array of length n
            # ARG 6: increment of index values of Y, maybe?
            BLAS.axpy!(nlambda, jcolumn, crosssection[jspecies][ialt+1], 1,
                       solarabs[ialt],1)
        end
    end

    #solarabs now records the total optical depth of the atmosphere at
    #each wavelength and altitude

    # actinic flux at each wavelength is solar flux diminished by total
    # optical depth
    for ialt in [1:nalt;]
        solarabs[ialt] = solarflux[:,2].*exp.(-solarabs[ialt])
    end

    # each species absorbs according to its cross section at each
    # altitude times the actinic flux.
    # BLAS.dot includes an integration (sum) across wavelengths, i.e:
    # (a·b) = aa + ab + ab + bb etc that kind of thing
    for j in Jratelist
        for ialt in [1:nalt;]
            n_current[j][ialt] = BLAS.dot(nlambda, solarabs[ialt], 1,
                                          crosssection[j][ialt+1], 1)
        end
       # this section for testing sensitivity to J rates
        # if contains(string(j), "H2O") | contains(string(j), "HDO")
        # if contains(string(j), "CO2toCOpO")
           # n_current[j] = n_current[j] ./ 10
        # end
    end
end

function timeupdate(mytime)
    for i = 1:15
        plotatm(n_current, t=mytime, iter=i)
        # println("dt: $(mytime)")
        update!(n_current, mytime)
    end
    # show()
end

function next_timestep(nstart::Array{Float64, 1}, nthis::Array{Float64, 1},
                       inactive::Array{Float64, 1}, activespecies, inactivespecies, 
                       Jrates::Array{Float64, 2},
                       T::Array{Float64, 1}, M::Array{Float64, 1},
                       tup::Array{Float64, 2}, tdown::Array{Float64, 2},
                       tlower::Array{Float64, 2}, tupper::Array{Float64, 2},
                       dt::Float64)
    #=
    moves to the next timestep using Newton's method on the linearized
    coupled transport and chemical reaction network.
    =#
    eps = 1.0 # ensure at least one iteration
    iter = 0
    while eps>1e-8
        nold = deepcopy(nthis)

        # stuff concentrations into update function and jacobian
        fval = nthis - nstart - dt*ratefn(nthis, inactive, inactivespecies, activespecies, Jrates, T, M, tup,
                                          tdown, tlower, tupper)
        updatemat = chemJmat(nthis, inactive, activespecies, inactivespecies, Jrates, T, M, tup, tdown, tlower,
                             tupper, dt)

        # update
        nthis = nthis - (updatemat \ fval)
        # check relative size of update
        eps = maximum(abs.(nthis-nold)./nold)
        iter += 1
        if iter>1e3; throw("too many iterations in next_timestep!"); end;
    end
    return nthis
end

function update!(n_current::Dict{Symbol, Array{Float64, 1}},dt)
    # update n_current using the coupled reaction network, moving to
    # the next timestep

    #set auxiliary (not solved for in chemistry) species values, photolysis rates
    inactive = deepcopy(Float64[[n_current[sp][ialt] for sp in inactivespecies, ialt in 1:length(intaltgrid)]...])
    Jrates = deepcopy(Float64[n_current[sp][ialt] for sp in Jratelist, ialt in 1:length(intaltgrid)])

    # extract concentrations
    nstart = deepcopy([[n_current[sp][ialt] for sp in activespecies, ialt in 1:length(intaltgrid)]...])
    M = sum([n_current[sp] for sp in fullspecieslist])

    # set temperatures
    T = Float64[Temp(a) for a in alt[2:end-1]]  

    # take initial guess
    nthis = deepcopy(nstart)

    # these are the sum of the transport flux coefficients D+K, divided by Δz², units 1/s
    tup = Float64[issubset([sp],notransportspecies) ? 0.0 : fluxcoefs(a, dz, sp, n_current)[2] for sp in specieslist, a in alt[2:end-1]]
    tdown = Float64[issubset([sp],notransportspecies) ? 0.0 : fluxcoefs(a, dz, sp, n_current)[1] for sp in specieslist, a in alt[2:end-1]]

    # put the lower layer and upper layer boundary conditions in separate arrays; but they are not the
    # right shape! 
    tlower_temporary = [boundaryconditions(sp, speciesbclist, dz, n_current)[1,:] for sp in specieslist]
    tupper_temporary = [boundaryconditions(sp, speciesbclist, dz, n_current)[2,:] for sp in specieslist]

    # reshape tlower and tupper into 2x2 arrays
    tlower = zeros(Float64, length(tlower_temporary), 2)
    tupper = zeros(Float64, length(tupper_temporary), 2)

    # tlower_temporary & tupper_temporary have same length; OK to use lower for the range
    for r in range(1, length=length(tlower_temporary))
        tlower[r, :] = tlower_temporary[r]
        tupper[r, :] = tupper_temporary[r]
    end

    # update to next timestep
    nthis = next_timestep(nstart, nthis, inactive, activespecies, inactivespecies, Jrates, T, M, tup, tdown,
                          tlower, tupper, dt)
    nthismat = reshape(nthis,(length(activespecies),length(intaltgrid)))

    # write found values out to n_current
    for s in 1:length(activespecies)
        for ia in 1:length(intaltgrid)
            tn = nthismat[s, ia]
            n_current[activespecies[s]][ia] = tn > 0. ? tn : 0.
        end
    end

    update_Jrates!(n_current)
end

################################################################################
#                                  MAIN SETUP                                  #
################################################################################

# Set up simulation files and experiment type ==================================

# Note: directory paths are in PARAMETERS.jl
# get command line arguments for experiment type. format:
# <experiment type> <parameters> <solar cycle type>
# examples: 
# temp Tsurf Ttropo Texo mean; water <mixing ratio> mean; dh <multiplier> mean; 
# Oflux <cm^-2s^-1> mean
# examples: temp 190 110 200 mean; water 1e-3 mean; dh 8 mean; Oflux 1.2e8 mean.
# last argument is solar cycle: min, mean, or max.
args = Any[ARGS[i] for i in 1:1:length(ARGS)]

# Establish a pattern for filenames. FNext = filename extension
if args[1]=="temp"
    FNext = "temp_$(args[2])_$(args[3])_$(args[4])"
elseif args[1]=="water"
    FNext = "water_$(args[2])"
elseif args[1]=="dh"
    FNext = "dh_$(args[2])"
elseif args[1]=="Oflux"
    FNext = "Oflux_$(args[2])"
else
    throw("Error! Bad experiment type")
end

# Set up the folder if it doesn't exist
create_folder(FNext, results_dir)

# Case where this file will be used to converge an atmosphere of a different extent
make_new_alt_grid = input("Would you like to use the script to converge a new atmosphere of a different extent? (y/n): ")
if make_new_alt_grid=="y"
    readfile = research_dir*"converged_200km_atmosphere.h5"
    const alt = convert(Array, (0:2e5:200e5))
    n_current = get_ncurrent(readfile)

    new_zmax = parse(Int64, input("Enter the new top of the atmosphere in km: "))
    extra_entries = Int64((new_zmax - 200)/(dz/1e5))

    # Extend the grid 
    for (k,v) in zip(keys(n_current), values(n_current))
       append!(v, fill(v[end], extra_entries))  # repeats the last value in the array for the upper atmo as an initial value.
    end

    const alt = convert(Array, (0:2e5:new_zmax*10^5))

    if new_zmax != 250
        println("Warning: You entered $(new_zmax) for the new max altitude but 
                 250 is hard-coded in the PARAMETERS.jl file. I haven't made this
                 general yet. So the code is probably about to break")
    end
elseif make_new_alt_grid=="n"
    # Set up the converged file to read from and load the simulation state at init.
    file_to_use = input("Enter the name of a file containing a converged, 250 km atmosphere to use (press enter to use default): ")
    readfile = file_to_use == "" ? "converged_250km_atmosphere.h5" : file_to_use
    n_current = get_ncurrent(readfile)
else
    throw("Didn't understand response")
end


# Set solar cycle file and alert user 
cycle = args[end]
solar_data_file = Dict("max"=>"marssolarphotonflux_solarmax.dat", 
                       "mean"=>"marssolarphotonflux_solarmean.dat", 
                       "min"=>"marssolarphotonflux_solarmin.dat")
solarfile = solar_data_file[cycle]

# Convert the arguments to numbers so we can use them to do maths
for i in 2:1:length(args)-1
    args[i] = parse(Float64, args[i])
end

# Let the user know what is being done 
println("ALERT: running sim for $(FNext)")
println("ALERT: Using file: ", readfile)
if cycle != "mean"
    println("ALERT: Solar $(cycle) data being used")
end

# Plot styles ==================================================================
rcParams = PyCall.PyDict(matplotlib."rcParams")
rcParams["font.sans-serif"] = ["Louis George Caf?"]
rcParams["font.monospace"] = ["FreeMono"]
rcParams["font.size"] = 22
rcParams["axes.labelsize"]= 24
rcParams["xtick.labelsize"] = 22
rcParams["ytick.labelsize"] = 22

################################################################################
#                             ADD DEUTERATED SPECIES                           # 
################################################################################

# D/H ratio ====================================================================
# General D/H ratio for mars, 5.5*SMOW, Atmosphere & Climate of Mars 2017 
DH = 5.5 * 1.6e-4
if args[1] == "dh"
    DH = parse(Float64, args[2]) * 1.6e-4
end

# modify n_current with deuterated species profiles ============================
n_current[:HDO] = n_current[:H2O] * DH
n_current[:OD] = n_current[:OH] * DH
n_current[:HDO2] = n_current[:H2O2] * DH
n_current[:D] = n_current[:H] * DH
n_current[:DO2] = n_current[:HO2] * DH
n_current[:HD] = n_current[:H2] * DH
n_current[:DOCO] = n_current[:HOCO] * DH

# add the new Jrates --the values will get populated automatically =============
n_current[:JHDOtoHpOD] = zeros(length(alt))
n_current[:JHDOtoDpOH] = zeros(length(alt))
n_current[:JHDO2toOHpOD] = zeros(length(alt))
n_current[:JHDOtoHDpO1D] = zeros(length(alt)) 
n_current[:JHDOtoHpDpO] = zeros(length(alt))
n_current[:JODtoOpD] = zeros(length(alt))
n_current[:JHDtoHpD] = zeros(length(alt))
n_current[:JDO2toODpO] = zeros(length(alt))
n_current[:JHDO2toDO2pH] = zeros(length(alt))
n_current[:JHDO2toHO2pD] = zeros(length(alt))
n_current[:JHDO2toHDOpO1D] = zeros(length(alt))
n_current[:JODtoO1DpD] = zeros(length(alt))

################################################################################
#                       TEMPERATURE/PRESSURE PROFILES                          #
################################################################################

println("Mean temperatures as entered in code:")
println("Surface: $(meanTs), Tropopause: $(meanTt), Exobase: $(meanTe)")

#If changes to the temperature are needed, they should be made here 
if args[1]=="temp"
    Temp(z::Float64) = Tpiecewise(z, args[2], args[3], args[4])
    Temp_keepSVP(z::Float64) = Tpiecewise(z, meanTs, meanTt, meanTe) # for testing temp without changing SVP.
else 
    Temp(z::Float64) = Tpiecewise(z, meanTs, meanTt, meanTe)
end

plot_temp_prof([Temp(a) for a in alt], results_dir*FNext, alt)


################################################################################
#                               WATER PROFILES                                 #
################################################################################

# set SVP to be fixed or variable with temperature
if args[1] == "temp"
    fix_SVP = true
    H2Osat = map(x->Psat(x), map(Temp_keepSVP, alt)) # for holding SVP fixed
    HDOsat = map(x->Psat_HDO(x), map(Temp_keepSVP, alt))  # for holding SVP fixed
    println("ALERT: SVP will be fixed to that of the mean temperature profile ")
else
    fix_SVP = false
    H2Osat = map(x->Psat(x), map(Temp, alt)) # array in #/cm^3 by altitude
    HDOsat = map(x->Psat_HDO(x), map(Temp, alt))
end

# H2O Water Profile ============================================================
surface_watersat = Dict("H2O"=>H2Osat[1], "HDO"=>HDOsat[1])
H2Osatfrac = H2Osat./map(z->n_tot(n_current, z), alt)  # get SVP as fraction of total atmo
# set H2O SVP fraction to minimum for all alts above first time min is reached
H2Oinitfrac = H2Osatfrac[1:something(findfirst(isequal(minimum(H2Osatfrac)), H2Osatfrac), 0)]
H2Oinitfrac = [H2Oinitfrac;   # ensures no supersaturation
               fill(minimum(H2Osatfrac), length(alt)-2-length(H2Oinitfrac))]

# make profile constant in the lower atmosphere (well-mixed).
# when doing water experiments, the temperature profile is such that the minimum
# in the SVP curve occurs at about 55 km alt. This was found by manual tweaking.
# thus we need to only set the lower atmo mixing ratio below that point, or 
# there will be a little spike in the water profile.
if args[1] == "water"
    H2Oinitfrac[findall(x->x<hygropause_alt, alt)] .= args[2]
    MR = args[2] # mixing ratio 
else  # i believe 30 km is supposed to be approximately the hygropause.
    MR = MR_mean_water
    H2Oinitfrac[findall(x->x<hygropause_alt, alt)] .= MR # 10 pr μm
end

for i in [1:length(H2Oinitfrac);]
    H2Oinitfrac[i] = H2Oinitfrac[i] < H2Osatfrac[i+1] ? H2Oinitfrac[i] : H2Osatfrac[i+1]
end

# HDO water profile ============================================================
HDOsatfrac = HDOsat./map(z->n_tot(n_current, z), alt)
# use D/H ratio to set population of HDO
HDOinitfrac = H2Oinitfrac * DH  # initial profile for HDO

# Compute total water column in pr μm starting with mixing ratio array =========
# H2O #/cm^3 = sum(MR * n_tot) for each alt
H2O_per_cc = sum([MR; H2Oinitfrac] .* map(z->n_tot(n_current, z), alt[1:end-1]))
HDO_per_cc = sum([MR*DH; HDOinitfrac] .* map(z->n_tot(n_current, z), alt[1:end-1]))

# pr μm = (H2O #/cm^3) * cm * (mol/#) * (H2O g/mol) * (1 cm^3/g) * (10^4 μm/cm)
# where the lone cm is a slice of the atmosphere of thickness dz, #/mol=6.023e23, 
# H2O or HDO g/mol = 18 or 19, cm^3/g = 1 or 19/18 for H2O or HDO.
# written as conversion factors for clarity.
H2Oprum = (H2O_per_cc * dz) * (18/1) * (1/6.02e23) * (1/1) * (1e4/1)
HDOprum = (HDO_per_cc * dz) * (19/1) * (1/6.02e23) * (19/18) * (1e4/1)

# Write out water content to a file ============================================
f = open(results_dir*FNext*"/water_column_"*FNext*".txt", "w")
write(f, "Total H2O col: $(H2O_per_cc*2e5)\n")
write(f, "Total HDO col: $(HDO_per_cc*2e5)\n")
write(f, "Total water col: $((H2O_per_cc + HDO_per_cc)*2e5)\n")
write(f, "H2O+HDO at surface: $((H2O_per_cc[1] + HDO_per_cc[1])*2e5)\n")
write(f, "Total H2O (pr μm): $(H2Oprum)\n")
write(f, "Total HDO (pr μm): $(HDOprum)\n")
write(f, "Total H2O+HDO, no enhancement: $(H2Oprum + HDOprum)")
close(f)

# Plot the water profiles (as raw mixing ratio) ================================
fig, ax = subplots(figsize=(6,9))
plot_bg(ax)
semilogx(H2Oinitfrac, alt[2:end-1]/1e5, color="cornflowerblue", linewidth=3, 
         label=L"H$_2$O")
semilogx(HDOinitfrac, alt[2:end-1]/1e5, color="cornflowerblue", linestyle="--", 
         linewidth=3, label="HDO")
semilogx(H2Osatfrac, alt[1:end]/1e5, color="black", alpha=0.5, linewidth=3, 
         label=L"H$_2$O saturation")
semilogx(HDOsatfrac, alt[1:end]/1e5, color="black", alpha=0.5, linestyle="--", 
         linewidth=3, label="HDO saturation")
xlabel("Mixing ratio", fontsize=18)
ylabel("Altitude [km]", fontsize=18)
title(L"H$_2$O and HDO model profiles", fontsize=20)
ax.tick_params("both",labelsize=16)
legend()
savefig(results_dir*FNext*"/water_profile_rawMR.png")

################################################################################
#                             BOUNDARY CONDITIONS                              #
################################################################################

H_effusion_velocity = effusion_velocity(Temp(zmax), 1.0, zmax)
H2_effusion_velocity = effusion_velocity(Temp(zmax), 2.0, zmax)
D_effusion_velocity = effusion_velocity(Temp(zmax), 2.0, zmax)
HD_effusion_velocity = effusion_velocity(Temp(zmax), 3.0, zmax)

# Used for passing a variable speciesbcs function
v_eff = Dict("H"=>H_effusion_velocity, "D"=>D_effusion_velocity, 
             "H2"=>H2_effusion_velocity, "HD"=>HD_effusion_velocity)

#=
    boundary conditions for each species (mostly from Nair 1994, Yung 1988). For 
    most species, default boundary condition is zero (f)lux at top and bottom. 
    Atomic/molecular hydrogen and deuterated analogues have a nonzero effusion 
    (v)elocity at the upper layer of the atmosphere. Some species have a (n)umber 
    density boundary condition.
=#

if args[1] == "Oflux"
    global const Of = args[2]
else
    global const Of = 1.2e8
end

# This has to be defined here, because it uses as a boundary condition the H2O 
# and HDO saturation at the surface.
global const speciesbclist=Dict(
                :CO2=>["n" 2.1e17; "f" 0.],
                :Ar=>["n" 2.0e-2*2.1e17; "f" 0.],
                :N2=>["n" 1.9e-2*2.1e17; "f" 0.],
                :H2O=>["n" H2Osat[1]; "f" 0.], # bc doesnt matter if H2O fixed
                :HDO=>["n" HDOsat[1]; "f" 0.],
                :O=>["f" 0.; "f" Of],
                :H2=>["f" 0.; "v" H2_effusion_velocity],
                :HD=>["f" 0.; "v" HD_effusion_velocity],
                :H=>["f" 0.; "v" H_effusion_velocity],
                :D=>["f" 0.; "v" D_effusion_velocity],
               );

################################################################################
#                       COMBINED CHEMISTRY AND TRANSPORT                       #
################################################################################

#=
    We now have objects that return the list of indices and coefficients
    for transport, assuming no other species in the atmosphere
    (transportmat), and for chemistry, assuming no other altitudes
    (chemical_jacobian). We need to perform a kind of outer product on
    these operators, to determine a fully coupled set of equations for
    all species at all altitudes.
=#

# need to get a list of all species at all altitudes to iterate over
const intaltgrid = round.(Int64, alt/1e5)[2:end-1]
const replacespecies = [fullspecieslist, Jratelist, [:T,:M]]

#=
    the rates at each altitude can be computed using the reaction network
    already in place, plus additional equations describing the transport
    to and from the cells above and below:
=#
upeqns = [Any[Any[[s], [Symbol(string(s)*"_above")],Symbol("t"*string(s)*"_up")],
          Any[[Symbol(string(s)*"_above")],[s],Symbol("t"*string(s)*"_above_down")]]
          for s in specieslist]

downeqns = [Any[Any[[s], [Symbol(string(s)*"_below")],Symbol("t"*string(s)*"_down")],
            Any[[Symbol(string(s)*"_below")],[s],Symbol("t"*string(s)*"_below_up")]]
            for s in specieslist]

local_transport_rates = [[[Symbol("t"*string(s)*"_up") for s in specieslist]
                          [Symbol("t"*string(s)*"_down") for s in specieslist]
                          [Symbol("t"*string(s)*"_above_down") for s in specieslist]
                          [Symbol("t"*string(s)*"_below_up") for s in specieslist]]...;]

transportnet = [[upeqns...;]; [downeqns...;]]

# define names for all the species active in the coupled rates:
activespecies = union(chemspecies, transportspecies)
active_above = [Symbol(string(s)*"_above") for s in activespecies]
active_below = [Symbol(string(s)*"_below") for s in activespecies]
inactivespecies = intersect(nochemspecies, notransportspecies)


# obtain the rates and jacobian for each altitude
const rates_local = Expr(:vcat, map(x->getrate(reactionnet, transportnet, x),activespecies)...);
const chemJ_local = chemical_jacobian(reactionnet, transportnet, activespecies, activespecies);
const chemJ_above = chemical_jacobian(reactionnet, transportnet, activespecies, active_above);
const chemJ_below = chemical_jacobian(reactionnet, transportnet, activespecies, active_below);

arglist_local = [activespecies; active_above; active_below; inactivespecies;
                 Jratelist; :T; :M; local_transport_rates; :dt]

arglist_local_typed=[:($s::Float64) for s in arglist_local]

@eval begin
    function ratefn_local($(arglist_local_typed[1:end-1]...))
        $rates_local # evaluates the rates_local expression
    end
end

@eval begin
    function chemJmat_local($(arglist_local_typed...))
        localchemJi = $(chemJ_local[1])
        localchemJj = $(chemJ_local[2])
        localchemJval = -dt*$(chemJ_local[3])

        abovechemJi = $(chemJ_above[1])
        abovechemJj = $(chemJ_above[2])
        abovechemJval = -dt*$(chemJ_above[3])

        belowchemJi = $(chemJ_below[1])
        belowchemJj = $(chemJ_below[2])
        belowchemJval = -dt*$(chemJ_below[3])

        ((localchemJi, localchemJj, localchemJval),
         (abovechemJi, abovechemJj, abovechemJval),
         (belowchemJi, belowchemJj, belowchemJval))
    end
end

@eval begin
    function reactionrates_local($(specieslist...), $(Jratelist...), T, M)
        #= a function to return chemical reaction rates for specified species
           concentrations =#
        $(Expr(:vcat, map(x->Expr(:call,:*,x[1]..., x[3]), reactionnet)...))
    end
end

################################################################################
#                         PHOTOCHEMICAL CROSS SECTIONS                         #
################################################################################

# Change following line as needed depending on local machine
xsecfolder = research_dir * "uvxsect/";

# Crosssection Files ===========================================================
co2file = "CO2.dat"
co2exfile = "binnedCO2e.csv" # added to shield short λ of sunlight in upper atmo
h2ofile = "h2oavgtbl.dat"
hdofile = "HDO.dat"#"HDO_250K.dat"# # TODO: change back
h2o2file = "H2O2.dat"
hdo2file = "H2O2.dat" #TODO: do HDO2 xsects exist?
o3file = "O3.dat"
o3chapfile = "O3Chap.dat"
o2file = "O2.dat"
o2_130_190 = "130-190.cf4"
o2_190_280 = "190-280.cf4"
o2_280_500 = "280-500.cf4"
h2file = "binnedH2.csv"
hdfile = "binnedH2.csv" # TODO: change this to HD file if xsects ever exist
ohfile = "binnedOH.csv"
oho1dfile = "binnedOHo1D.csv"
odfile = "OD.csv"

# Loading Data =================================================================
# CO2 --------------------------------------------------------------------------
# temperature-dependent between 195-295K
co2xdata = readandskip(xsecfolder*co2file,'\t',Float64, skipstart = 4)

# CO2 photoionization (used to screen high energy sunlight)
co2exdata = readandskip(xsecfolder*co2exfile,',',Float64, skipstart = 4)

# H2O & HDO --------------------------------------------------------------------
h2oxdata = readandskip(xsecfolder*h2ofile,'\t',Float64, skipstart = 4)

# These crosssections for HDO are for 298K.
hdoxdata = readandskip(xsecfolder*hdofile,'\t', Float64, skipstart=12)

# H2O2 + HDO2 ------------------------------------------------------------------
# the data in the following table cover the range 190-260nm
h2o2xdata = readandskip(xsecfolder*h2o2file,'\t',Float64, skipstart=3)
hdo2xdata = readandskip(xsecfolder*hdo2file,'\t',Float64, skipstart=3)

# O3 ---------------------------------------------------------------------------
# including IR bands which must be resampled from wavenumber
o3xdata = readandskip(xsecfolder*o3file,'\t',Float64, skipstart=3)
o3ls = o3xdata[:,1]
o3xs = o3xdata[:,2]
o3chapxdata = readandskip(xsecfolder*o3chapfile,'\t',Float64, skipstart=3)
o3chapxdata[:,1] = map(p->1e7/p, o3chapxdata[:,1])
for i in [round(Int, floor(minimum(o3chapxdata[:,1]))):round(Int, ceil(maximum(o3chapxdata))-1);]
    posss = getpos(o3chapxdata, x->i<x<i+1)
    dl = diff([map(x->o3chapxdata[x[1],1], posss); i])
    x = map(x->o3chapxdata[x[1],2],posss)
    ax = reduce(+,map(*,x, dl))/reduce(+,dl)
    global o3ls = [o3ls; i+0.5]
    global o3xs = [o3xs; ax]
end
o3xdata = reshape([o3ls; o3xs],length(o3ls),2)

# O2 ---------------------------------------------------------------------------
o2xdata = readandskip(xsecfolder*o2file,'\t',Float64, skipstart = 3)
o2schr130K = readandskip(xsecfolder*o2_130_190,'\t',Float64, skipstart = 3)
o2schr130K[:,1] = map(p->1e7/p, o2schr130K[:,1])
o2schr130K = binupO2(o2schr130K)
o2schr190K = readandskip(xsecfolder*o2_190_280,'\t',Float64, skipstart = 3)
o2schr190K[:,1] = map(p->1e7/p, o2schr190K[:,1])
o2schr190K = binupO2(o2schr190K)
o2schr280K = readandskip(xsecfolder*o2_280_500,'\t',Float64, skipstart = 3)
o2schr280K[:,1] = map(p->1e7/p, o2schr280K[:,1])
o2schr280K = binupO2(o2schr280K)

# HO2 & DO2 --------------------------------------------------------------------
ho2xsect = [190.5:249.5;]
ho2xsect = reshape([ho2xsect; map(ho2xsect_l, ho2xsect)],length(ho2xsect),2)
do2xsect = deepcopy(ho2xsect)

# H2 & HD ----------------------------------------------------------------------
h2xdata = readandskip(xsecfolder*h2file,',',Float64, skipstart=4)
hdxdata = readandskip(xsecfolder*hdfile,',',Float64, skipstart=4)

# OH & OD ----------------------------------------------------------------------
ohxdata = readandskip(xsecfolder*ohfile,',',Float64, skipstart=4)
ohO1Dxdata = readandskip(xsecfolder*oho1dfile,',',Float64, skipstart=4)
odxdata = readandskip(xsecfolder*odfile,',',Float64, skipstart=3)

# PHOTODISSOCIATION ============================================================

const solarflux=readandskip(research_dir*solarfile,'\t', Float64,skipstart=4)[1:2000,:]
solarflux[:,2] = solarflux[:,2]/2  # Account for only doing the dayside (I think?)

absorber = Dict(:JCO2ion =>:CO2,
                :JCO2toCOpO =>:CO2,
                :JCO2toCOpO1D =>:CO2,
                :JO2toOpO =>:O2,
                :JO2toOpO1D =>:O2,
                :JO3toO2pO =>:O3,
                :JO3toO2pO1D =>:O3,
                :JO3toOpOpO =>:O3,
                :JH2toHpH =>:H2,
                :JHDtoHpD => :HD,
                :JOHtoOpH =>:OH,
                :JOHtoO1DpH =>:OH,
                :JODtoOpD =>:OD,
                :JODtoO1DpD => :OD,
                :JHO2toOHpO =>:HO2,
                :JDO2toODpO => :DO2,
                :JH2OtoHpOH =>:H2O,
                :JH2OtoH2pO1D =>:H2O,
                :JH2OtoHpHpO =>:H2O,
                :JH2O2to2OH =>:H2O2,
                :JH2O2toHO2pH =>:H2O2,
                :JH2O2toH2OpO1D =>:H2O2,
                :JHDO2toHDOpO1D => :HDO2,
                :JHDOtoHpOD=>:HDO,
                :JHDO2toOHpOD=>:HDO2,
                :JHDO2toDO2pH => :HDO2,
                :JHDO2toHO2pD => :HDO2,
                :JHDOtoDpOH=>:HDO,
                :JHDOtoHpDpO=>:HDO,
                :JHDOtoHDpO1D=>:HDO
                );

#=
    this is a dictionary of the 1-nm photodissociation or photoionization
    cross-sections important in the atmosphere. keys are symbols found in
    jratelist. each entry is an array of arrays, yielding the wavelengths
    and cross-sections for each altitude in the atmosphere.

    NOTE: jspecies refers to the photodissociation or photoionization
    cross section for a particular species which produces a UNIQUE SET OF
    PRODUCTS. In this sense, crosssection has already folded in quantum
    efficiency considerations.
=#
crosssection = Dict{Symbol, Array{Array{Float64}}}()

# CO2 photodissociation --------------------------------------------------------
setindex!(crosssection, fill(co2exdata, length(alt)), :JCO2ion)
#CO2+hv->CO+O
setindex!(crosssection,
          map(xs->quantumyield(xs,((l->l>167, 1), (l->95>l, 0.5))),
          map(t->co2xsect(co2xdata, t),map(Temp, alt))), :JCO2toCOpO)
#CO2+hv->CO+O1D
setindex!(crosssection,
          map(xs->quantumyield(xs,((l->95<l<167, 1), (l->l<95, 0.5))),
          map(t->co2xsect(co2xdata, t),map(Temp, alt))), :JCO2toCOpO1D)

# O2 photodissociation ---------------------------------------------------------
#O2+hv->O+O
setindex!(crosssection,
          map(xs->quantumyield(xs,((x->x>175, 1),)), map(t->o2xsect(o2xdata, o2schr130K, o2schr190K, o2schr280K, t),map(Temp, alt))),
          :JO2toOpO)
#O2+hv->O+O1D
setindex!(crosssection,
          map(xs->quantumyield(xs,((x->x<175, 1),)), map(t->o2xsect(o2xdata, o2schr130K, o2schr190K, o2schr280K, t),map(Temp, alt))),
          :JO2toOpO1D)

# O3 photodissociation ---------------------------------------------------------
# O3+hv->O2+O
setindex!(crosssection,
          map(t->quantumyield(o3xdata,
                              (
                               (l->l<193, 1-(1.37e-2*193-2.16)),
                               (l->193<=l<225, l->(1 .- (1.37e-2*l-2.16))),
                               (l->225<=l<306, 0.1),
                               (l->306<=l<328, l->(1 .- O3O1Dquantumyield(l, t))),
                               (l->328<=l<340, 0.92),
                               (l->340<=l, 1.0)
                              )), map(Temp, alt)), :JO3toO2pO)
# O3+hv->O2+O1D
setindex!(crosssection,
          map(t->quantumyield(o3xdata,
                              (
                               (l->l<193, 1.37e-2*193-2.16),
                               (l->193<=l<225, l->(1.37e-2*l-2.16)),
                               (l->225<=l<306, 0.9),
                               (l->306<=l<328, l->O3O1Dquantumyield(l, t)),
                               (l->328<=l<340, 0.08),
                               (l->340<=l, 0.0)
                              )), map(Temp, alt)), :JO3toO2pO1D)
# O3+hv->O+O+O
setindex!(crosssection,
          fill(quantumyield(o3xdata,((x->true, 0.),)),length(alt)),
          :JO3toOpOpO)

# H2 and HD photodissociation --------------------------------------------------
# H2+hv->H+H
setindex!(crosssection, fill(h2xdata, length(alt)), :JH2toHpH)
# HD+hν -> H+D 
setindex!(crosssection, fill(hdxdata, length(alt)), :JHDtoHpD)

# OH and OD photodissociation --------------------------------------------------
# OH+hv->O+H
setindex!(crosssection, fill(ohxdata, length(alt)), :JOHtoOpH)
# OH+hv->O1D+H
setindex!(crosssection, fill(ohO1Dxdata, length(alt)), :JOHtoO1DpH)
# OD + hv -> O+D  
setindex!(crosssection, fill(odxdata, length(alt)), :JODtoOpD)
# OD + hν -> O(¹D) + D 
setindex!(crosssection, fill(ohO1Dxdata, length(alt)), :JODtoO1DpD)

# HO2 and DO2 photodissociation ------------------------------------------------
# HO2 + hν -> OH + O
setindex!(crosssection, fill(ho2xsect, length(alt)), :JHO2toOHpO)
# DO2 + hν -> OD + O
setindex!(crosssection, fill(do2xsect, length(alt)), :JDO2toODpO)

# H2O and HDO photodissociation ------------------------------------------------
# H2O+hv->H+OH
setindex!(crosssection,
          fill(quantumyield(h2oxdata,((x->x<145, 0.89),(x->x>145, 1))),length(alt)),
          :JH2OtoHpOH)

# H2O+hv->H2+O1D
setindex!(crosssection,
          fill(quantumyield(h2oxdata,((x->x<145, 0.11),(x->x>145, 0))),length(alt)),
          :JH2OtoH2pO1D)

# H2O+hv->H+H+O
setindex!(crosssection,
          fill(quantumyield(h2oxdata,((x->true, 0),)),length(alt)),
          :JH2OtoHpHpO)

# HDO + hν -> H + OD
setindex!(crosssection,
          fill(quantumyield(hdoxdata,((x->x<145, 0.5*0.89),(x->x>145, 0.5*1))),length(alt)),
          :JHDOtoHpOD)

# HDO + hν -> D + OH
setindex!(crosssection,
          fill(quantumyield(hdoxdata,((x->x<145, 0.5*0.89),(x->x>145, 0.5*1))),length(alt)),
          :JHDOtoDpOH)

# HDO + hν -> HD + O1D
setindex!(crosssection,
          fill(quantumyield(hdoxdata,((x->x<145, 0.11),(x->x>145, 0))),length(alt)),
          :JHDOtoHDpO1D)

# HDO + hν -> H + D + O
setindex!(crosssection,
          fill(quantumyield(hdoxdata,((x->true, 0),)),length(alt)),
          :JHDOtoHpDpO)


# H2O2 and HDO2 photodissociation ----------------------------------------------
setindex!(crosssection,
          map(xs->quantumyield(xs,((x->x<230, 0.85),(x->x>230, 1))),
          map(t->h2o2xsect(h2o2xdata, t), map(Temp, alt))), :JH2O2to2OH)

# H2O2+hv->HO2+H
setindex!(crosssection,
          map(xs->quantumyield(xs,((x->x<230, 0.15),(x->x>230, 0))),
          map(t->h2o2xsect(h2o2xdata, t), map(Temp, alt))), :JH2O2toHO2pH)

# H2O2+hv->H2O+O1D
setindex!(crosssection,
          map(xs->quantumyield(xs,((x->true, 0),)), map(t->h2o2xsect(h2o2xdata, t),
          map(Temp, alt))), :JH2O2toH2OpO1D)

# HDO2 + hν -> OH + OD
setindex!(crosssection,
          map(xs->quantumyield(xs,((x->x<230, 0.85),(x->x>230, 1))),
          map(t->hdo2xsect(hdo2xdata, t), map(Temp, alt))), :JHDO2toOHpOD)

# HDO2 + hν-> DO2 + H
setindex!(crosssection,
          map(xs->quantumyield(xs,((x->x<230, 0.5*0.15),(x->x>230, 0))),
          map(t->hdo2xsect(hdo2xdata, t), map(Temp, alt))), :JHDO2toDO2pH)

# HDO2 + hν-> HO2 + D
setindex!(crosssection,
          map(xs->quantumyield(xs,((x->x<230, 0.5*0.15),(x->x>230, 0))),
          map(t->hdo2xsect(hdo2xdata, t), map(Temp, alt))), :JHDO2toHO2pD)

# HDO2 + hν -> HDO + O1D
setindex!(crosssection,
          map(xs->quantumyield(xs,((x->true, 0),)), map(t->hdo2xsect(hdo2xdata, t),
          map(Temp, alt))), :JHDO2toHDOpO1D)

# Solar Input ------------------------------------------------------------------
lambdas = Float64[]
for j in Jratelist, ialt in 1:length(alt)
    global lambdas = union(lambdas, crosssection[j][ialt][:,1])
end

if !(setdiff(solarflux[:,1],lambdas)==[])
    throw("Need a broader range of solar flux values!")
end

# pad all cross-sections to solar
for j in Jratelist, ialt in 1:length(alt)
    crosssection[j][ialt] = padtosolar(solarflux, crosssection[j][ialt])
end

# we need some global objects for the Jrates calculation:
# intensity as a function of wavelength at each altitude

# this is the unitialized array for storing values
solarabs = fill(fill(0.,size(solarflux, 1)),length(alt)-2);

################################################################################
#                             CONVERGENCE CODE                                 #
################################################################################

# set the water profiles =======================================================
n_current[:H2O] = H2Oinitfrac.*map(z->n_tot(n_current, z), alt[2:end-1])
n_current[:HDO] = HDOinitfrac.*map(z->n_tot(n_current, z), alt[2:end-1])

# Plot initial water profile ===================================================
# ...in ppm --------------------------------------------------------------------
fig, ax = subplots(figsize=(6,8))
plot_bg(ax)
ax.tick_params(axis="x", which="minor", bottom=true, top=true)
semilogx(H2Oinitfrac/1e-6, alt[2:end-1]/1e5, color="cornflowerblue", linewidth=2,
         label=L"H$_2$O")
semilogx(HDOinitfrac/1e-6, alt[2:end-1]/1e5, color="cornflowerblue", 
         linestyle="--", linewidth=2, label="HDO")
xlabel("Volume Mixing Ratio [ppm]")
ylabel("Altitude [km]")
title(L"H$_2$O and HDO model profiles")
legend()
savefig(results_dir*FNext*"/initial_water_MR.png")

# ...and raw abundance ---------------------------------------------------------
fig, ax = subplots(figsize=(6,8))
plot_bg(ax)
ax.tick_params(axis="x", which="minor", bottom=true, top=true)
semilogx(n_current[:H2O], alt[2:end-1]/1e5, color="cornflowerblue", linewidth=2,
         label=L"H$_2$O")
semilogx(n_current[:HDO], alt[2:end-1]/1e5, color="cornflowerblue", 
         linestyle="--", linewidth=2, label="HDO")
xlabel(L"Species density [cm$^{-2}$")
ylabel("Altitude [km]")
title(L"H$_2$O and HDO model profiles")
legend()
savefig(results_dir*FNext*"/initial_water_number.png")
#show()

# do the convergence ===========================================================
# initialize whole atmosphere figure 
convfig, convax = subplots(figsize=(8,6))
# convfig.show()
# convfig.canvas.draw()

@showprogress 0.1 "Converging over 10 My..." [timeupdate(t) for t in [10.0^(1.0*i) for i in -3:14]]
@showprogress 0.1 "Last convergence steps..." for i in 1:100
    plotatm(n_current, t="1e14", iter=i)
    # println("dt: 1e14 iter $(i)")
    update!(n_current, 1e14)
end

# write out the new converged file to matching folder. 
towrite = results_dir*FNext*"/converged_"*FNext*".h5"
write_ncurrent(n_current, towrite)
println("Wrote $(towrite)")

# save the figure
savefig(results_dir*FNext*"/converged_"*FNext*".png", bbox_inches="tight")
println("Saved figure to same folder")

################################################################################
#                                 LOGGING                                      #
################################################################################

# crosssection dict for logging purposes =======================================
xsect_dict = Dict("CO2"=>[co2file, co2exfile], 
                  "H2O, HDO"=>[h2ofile, hdofile],
                  "H2O2, HDO2"=>[h2o2file, hdo2file],
                  "O3"=>[o3file, o3chapfile],
                  "O2"=>[o2file, o2_130_190, o2_190_280, o2_280_500],
                  "H2, HD"=>[h2file, hdfile],
                  "OH, OD"=>[ohfile, oho1dfile, odfile])

# Log temperature and water parameters =========================================
if args[1]=="temp"
    input_string = "T_0=$(args[2]), T_tropo=$(args[3]), T_exo=$(args[4])" * 
                   "\nwater init=$(MR)\nDH=5.5 \nOflux=1.2e8" #*
                   #"\nlapse rate=$(lapserate_logme)\n"
elseif args[1]=="water"
    input_string = "T_s=$(meanTs), T_tropo=$(meanTt), T_exo=$(meanTe)\n" *
                   "water init=$(args[2])\nDH=5.5\nOflux=1.2e8\n" #*
                   #"lapse rate=$(lapserate_logme)\n"
elseif args[1]=="dh"
    input_string = "T_s=$(meanTs), T_tropo=$(meanTt), T_exo=$(meanTe)\nwater=(MR)\n" *
                   "DH=$(args[2]) \nOflux=1.2e8\n"#lapse rate=$(lapserate_logme)\n"
elseif args[1]=="Oflux"
    input_string = "T_s=$(meanTs), T_tropo=$(meanTt), T_exo=$(meanTe)\nwater=(MR)" *
                   "\nDH=5.5\nOflux=$(Of)\n"#lapse rate=$(lapserate_logme)\n"
end

# Write the log ================================================================
f = open(results_dir*FNext*"/convergence_data_"*FNext*".txt", "w")
write(f, "Finished convergence for $(args[1]) experiment with control parameters: \n")
write(f, input_string)
write(f, "\nSVP fixed: $(fix_SVP)\n")
write(f, "Solar cycle status: solar $(cycle)\n")


write(f, "\nCROSS SECTIONS: \n")
for k in keys(xsect_dict)  # cross sections
    write(f, k*": "*join(xsect_dict[k], ", ")*"\n")
end
# boundary conditions
write(f, "\nBOUNDARY CONDITIONS: \n")
write(f, "n: number density at surface, f: flux at top, v: velocity at top\n")
for k2 in keys(speciesbclist)
    bcstring = join([join(speciesbclist[k2][1, :], "="), 
                     join(speciesbclist[k2][2, :], "=")], ", ")
    write(f, string(k2)*": "*bcstring*"\n")
end

close(f)

println("ALERT: Finished")
println()