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

t1 = time()

using Revise
using PyPlot
using PyCall
using HDF5, JLD
using LaTeXStrings
using Distributed
using DelimitedFiles
using SparseArrays
using LinearAlgebra
using ProgressMeter
using Photochemistry  # custom module

t2 = time()

println("Time to load modules: $(round(t2-t1, digits=1)) seconds")

user_input_paramfile = input("Enter a parameter file or press enter to use default (PARAMETERS.jl): ")
paramfile = user_input_paramfile == "" ? "PARAMETERS.jl" : user_input_paramfile
t3 = time()
include(paramfile)
t4 = time()
println("Time to load PARAMETERS: $(round(t4-t3, digits=1)) seconds")
t5 = time()

################################################################################
#                                   FUNCTIONS                                  #
################################################################################

function ratefn(nthis, inactive, inactivespecies, activespecies, Jrates, T, tup, tdown, tlower, tupper)
    #=
    at each altitude, get the appropriate group of concentrations,
    coefficients, and rates to pass to ratefn_local
    =#
    nthismat = reshape(nthis, (length(activespecies), num_layers))
    inactivemat = reshape(inactive, (length(inactivespecies), num_layers))
    returnrates = zero(nthismat)

    # fill the first altitude entry with information for all species
    returnrates[:,1] .= ratefn_local([nthismat[:,1]; nthismat[:,2];
                                     fill(1.0, length(activespecies));
                                     inactivemat[:,1]; Jrates[:,1]; 
                                     T[1]; tup[:,1]; tlower[:,1]; tdown[:,2]; tlower[:,2]]...)

    # iterate through other altitudes in the lower atmosphere
    for ialt in 2:(num_layers-1)
        returnrates[:,ialt] .= ratefn_local([nthismat[:,ialt];
                                            nthismat[:,ialt+1];
                                            nthismat[:,ialt-1];
                                            inactivemat[:,ialt];
                                            Jrates[:,ialt];
                                            T[ialt]; 
                                            tup[:,ialt]; 
                                            tdown[:,ialt];
                                            tdown[:,ialt+1]; 
                                            tup[:,ialt-1]]...)
    end

    # fill in the last level of altitude
    returnrates[:,end] .= ratefn_local([nthismat[:,end];
                                       fill(1.0, length(activespecies));
                                       nthismat[:,end-1];
                                       inactivemat[:,end];
                                       Jrates[:,end];
                                       T[end]; 
                                       tupper[:,1]; 
                                       tdown[:,end];
                                       tupper[:,2]; 
                                       tup[:,end-1]]...)

    # NEW: Overwrite the entries for water in the lower atmosphere with 0s so that it will behave as fixed.
    # Only runs when water is in the activespecies list. If neutrals are set to inactive, it will be taken care of already.
    if in(:H2O, activespecies) && in(:HDO, activespecies)
        returnrates[H2Oi, 1:upper_lower_bdy_i] .= 0
        returnrates[HDOi, 1:upper_lower_bdy_i] .= 0
    end
    
    return [returnrates...;]
end

function chemJmat(nthis, inactive, activespecies, inactivespecies, Jrates, T, tup, tdown, tlower, tupper, dt)
    #=
    Collects coordinate tuples of (I, J, V) [row index, column index, value] for a sparse matrix
    representing the chemical jacobian of the atmospheric system. 

    nthis: The atmospheric densities array, but flattened, in the form [n_CO(z=0), n_CO2(z=0)...n_N2Dpl(z=0), n_CO(z=2)...n_N2Dpl(z=250)]
    inactive: A flattened array of the atmospheric densities of any inactive species, same format as nthis. Functionally constant.
    activespecies: List of active species (const)
    inactivespecies: List of inactivespecies (const)
    Jrates: Flattened array of Jrates, same format as nthis.
    Tn, Ti, Te: Temperature-vs-altitude arrays for neutrals, ions, electrons. (const)
    tup, tdown: Transport coefficients
    tlower, tupper: Transport coefficients

    dt: Can be specified to manually control the timestep by which the derivatives are multiplied.
        If not supplied, then dt=1 so that the external solver can manage the multiplication by time.
    =#

    nthismat = reshape(nthis, (length(activespecies), num_layers))
    inactivemat = reshape(inactive, (length(inactivespecies), num_layers))
    chemJi = Int64[]
    chemJj = Int64[]
    chemJval = Float64[]

    # tc___ are the coordinate tuples containing (I, J, V) to be used to fill a sparse matrix.
    (tclocal, tcupper, tclower) = chemJmat_local([nthismat[:,1]; nthismat[:,2]; fill(1.0, length(activespecies));
                                                  inactivemat[:,1]; Jrates[:,1];
                                                  T[1]; 
                                                  tup[:,1]; tlower[:,1];
                                                  tdown[:,2]; tlower[:,2]; dt]...) 


    # add the influence of the local densities
    append!(chemJi, tclocal[1])
    append!(chemJj, tclocal[2])
    append!(chemJval, tclocal[3])

    # and the upper densities
    append!(chemJi, tcupper[1])
    append!(chemJj, tcupper[2] .+ length(activespecies))
    append!(chemJval, tcupper[3])

    for ialt in 2:(num_layers-1)
        (tclocal, tcupper, tclower) = chemJmat_local([nthismat[:,ialt];
                                                      nthismat[:,ialt+1];
                                                      nthismat[:,ialt-1];
                                                      inactivemat[:,ialt];
                                                      Jrates[:,ialt]; T[ialt]; 
                                                      tup[:,ialt];
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
                                              T[end];
                                              tupper[:,1]; tdown[:,end];
                                              tupper[:,2]; tup[:,end-1]; dt]...)

    # add the influence of the local densities
    append!(chemJi, tclocal[1].+(num_layers-1)*length(activespecies))
    append!(chemJj, tclocal[2].+(num_layers-1)*length(activespecies))
    append!(chemJval, tclocal[3])

    # and the lower densities
    append!(chemJi, tclower[1].+(num_layers-1)*length(activespecies))
    append!(chemJj, tclower[2].+(num_layers-2)*length(activespecies))
    append!(chemJval, tclower[3])

    # NEW: fix water below whatever we set as upper/lower atmosphere boundary.
    # This only runs if water is designated as an active species; if it's in inactivespecies, this doesn't need to run,
    # and in fact, CAN'T run. When it is active, this finds all the H2O and HDO indices for the lower atmosphere. 
    # It's like above where we add (ialt-1)*length(activespecies), but this way it's outside the loop.
    if in(:H2O, activespecies) && in(:HDO, activespecies)
        H2Opositions = H2Oi .+ length(activespecies)*collect(0:upper_lower_bdy_i-1) # if H2Oi is 5, this returns 5, 69, 133...
        HDOpositions = HDOi .+ length(activespecies)*collect(0:upper_lower_bdy_i-1)
        water_positions = sort(union(H2Opositions, HDOpositions))

        i_remove = findall(x->in(x, water_positions), chemJi)
        j_remove = findall(x->in(x, water_positions), chemJj)
        remove_these = sort(union(i_remove, j_remove)) # This makes a set, since it describes the locations where the H2O and HDO indices are.
                                                       # Kinda confusing since we're talking about indices of indices.s
        deleteat!(chemJi, remove_these)
        deleteat!(chemJj, remove_these)
        deleteat!(chemJval, remove_these)
    end
    
    # make sure to add 1's along the diagonal
    append!(chemJi,[1:length(nthis);])
    append!(chemJj,[1:length(nthis);])
    append!(chemJval, fill(1.0, length(nthis))) 

    return sparse(chemJi, chemJj, chemJval, length(nthis), length(nthis), +);
end


# main routine functions =======================================================
function update_Jrates!(n_cur_densities::Dict{Symbol, Array{Float64, 1}})
    #=
    this function updates the photolysis rates stored in n_current to
    reflect the altitude distribution of absorbing species
    =#

    # Initialize an array, length=number of active layers
    # Each sub-array is an array of length 2000, corresponding to 2000 wavelengths.
    solarabs = Array{Array{Float64}}(undef, num_layers)
    for i in range(1, length=num_layers)
        solarabs[i] = zeros(Float64, 2000)
    end

    nalt = size(solarabs, 1)
    nlambda = size(solarabs[1],1)

    for jspecies in Jratelist
        species = absorber[jspecies]

        jcolumn = 0.
        for ialt in [nalt:-1:1;]
            #get the vertical column of the absorbing constituient
            jcolumn += n_cur_densities[species][ialt]*dz

            # add the total extinction to solarabs:
            # multiplies air column density (N, #/cm^2) at all wavelengths by crosssection (σ)
            # to get optical depth (τ). This is an override of axpy! to use the
            # full arguments. For the equation Y' = alpha*X + Y:
            # ARG 1: n (length of arrays in ARGS 3, 5)
            # ARG 2: alpha, a scalar.
            # ARG 3: X, an array of length n.
            # ARG 4: the increment of the index values of X, maybe?
            # ARG 5: Y, an array of length n
            # ARG 6: increment of index values of Y, maybe?
            BLAS.axpy!(nlambda, jcolumn, crosssection[jspecies][ialt+1], 1, solarabs[ialt], 1)
        end
    end

    # solarabs now records the total optical depth of the atmosphere at
    # each wavelength and altitude

    # actinic flux at each wavelength is solar flux diminished by total
    # optical depth
    for ialt in [1:nalt;]
        solarabs[ialt] = solarflux[:,2] .* exp.(-solarabs[ialt])
    end

    # each species absorbs according to its cross section at each
    # altitude times the actinic flux.
    # BLAS.dot includes an integration (sum) across wavelengths, i.e:
    # (a·b) = aa + ab + ab + bb etc that kind of thing
    for j in Jratelist
        n_cur_densities[j] .= 0. .* similar(n_cur_densities[:CO2])
        for ialt in [1:nalt;]
            n_cur_densities[j][ialt] = BLAS.dot(nlambda, solarabs[ialt], 1, crosssection[j][ialt+1], 1)
        end
    end
end

function timeupdate(mytime, thefolder, D_arr)
    # plot rates to inspect them! Used while converging an atmosphere.

    numiters = 15
    for i = 1:numiters
        update!(n_current, mytime, D_arr)
    end

    # plot_atm(n_current, [fullspecieslist], thefolder*"/atm_dt=$(mytime).png", t=mytime)  # Turned off to improve speed. Can be run in separate script.

    # Write out the current atmospheric state for making plots later 
    write_ncurrent(n_current, thefolder*"/ncurrent_$(mytime).h5") 

end

function next_timestep(nstart::Array{Float64, 1}, nthis::Array{Float64, 1},
                       inactive::Array{Float64, 1}, activespecies, inactivespecies,
                       Jrates::Array{Float64, 2},
                       T::Array{Float64, 1},
                       tup::Array{Float64, 2}, tdown::Array{Float64, 2},
                       tlower::Array{Float64, 2}, tupper::Array{Float64, 2},
                       dt::Float64, abs_tol_vec, rel_tol)
    #=
    moves to the next timestep using Newton's method on the linearized
    coupled transport and chemical reaction network.
    =#
    eps = 1.0 # ensure at least one iteration
    iter = 0

    # Set up the array to track whether all elements have met criterion
    met_criterion = BitArray(undef, length(nthis))

    while !all(met_criterion) # will quit when all elements of met_criterion are true.
        nold = deepcopy(nthis)

        # stuff concentrations into update function and jacobian. This is the Eulerian (P*L) * dt
        fval = nthis - nstart - dt*ratefn(nthis, inactive, inactivespecies, activespecies, Jrates, T, tup, tdown, tlower, tupper) 
        
        updatemat = chemJmat(nthis, inactive, activespecies, inactivespecies, Jrates, T, tup, tdown, tlower, tupper, dt)

        nthis = nthis - (updatemat \ fval)

        # Check whether error tolerance is met         
        # abs_eps = abs.(nthis-nold) ./ abs_tol_vec  # calculated for every species.
        # rel_eps = (abs.(nthis-nold)./nold) ./ rel_tol

        eps = abs.((nthis-nold) ./ (abs_tol_vec .+ nold .* rel_tol))

        # Check each element to see if < 1. 
        # abs_bool = abs_eps .< 1
        # rel_bool = rel_eps .< 1

        # Assign values to the array that governs when the loop runs.
        met_criterion = eps .< 1#abs_bool .| rel_bool

        # if all(met_criterion)
        #     println("All elements of nthis have met the tolerance after $(iter) iterations")
        # end

        iter += 1
        if iter>1e3
            throw("too many iterations in next_timestep! number of elements that met the criteria: $(count(met_criterion))/$(length(met_criterion))")
        end
    end
    return nthis
end

function update!(n_current::Dict{Symbol, Array{Float64, 1}}, dt, D_arr; plotJratesflag=false, iter=nothing)
    # update n_current using the coupled reaction network, moving to
    # the next timestep

    #set auxiliary (not solved for in chemistry) species values, photolysis rates
    # inactive = deepcopy(Float64[[n_current[sp][ialt] for sp in inactivespecies, ialt in 1:length(intaltgrid)]...])
    Jrates = deepcopy(Float64[n_current[sp][ialt] for sp in Jratelist, ialt in 1:length(non_bdy_layers)])

    # extract concentrations
    nstart = flatten_atm(n_current, activespecies)
    # nstart = deepcopy([[n_current[sp][ialt] for sp in activespecies, ialt in 1:length(intaltgrid)]...])
    # Next line calculates 1 ppt of the total density at each altitude - used for absolute error tolerance.
    # This gets it into the same shape as nstart.
    ppt_vals = 1e-12 .* [[n_tot(n_current, a) for sp in activespecies, a in non_bdy_layers]...]
    rel_tol = 1e-6

    # set temperatures
    # T = Float64[Temp(a) for a in non_bdy_layers]

    # take initial guess
    nthis = deepcopy(nstart)

    # Make the transport coefficients, diffusion coefficients, and scale height dicts
    Keddy_arr = similar(alt)
    Keddy_arr .= map(z->Keddy(z, n_tot(n_current, z)), alt) # Eddy diffusion: K coefficient by altitude (array)
    Dcoef_dict = Dict{Symbol, Vector{Float64}}([s=>deepcopy(Dcoef!(D_arr, Tn_arr, s, n_current)) for s in fullspecieslist])
    H0_dict = Dict{String, Vector{Float64}}("neutral"=>map((z,t)->scaleH(z, t, n_current), alt, Tn_arr),
                                            "ion"=>map((z,t)->scaleH(z, t, n_current), alt, Tn_arr)) # NOTE: ions will be wrong

    # these are the sum of the transport flux coefficients D+K, divided by Δz², units 1/s
    # tup = Float64[issubset([sp],notransportspecies) ? 0.0 : fluxcoefs(a, dz, sp, n_current, [T_surf, T_tropo, T_exo])[2] for sp in specieslist, a in non_bdy_layers]
    # tdown = Float64[issubset([sp],notransportspecies) ? 0.0 : fluxcoefs(a, dz, sp, n_current, [T_surf, T_tropo, T_exo])[1] for sp in specieslist, a in non_bdy_layers]
    fluxcoefs_all = fluxcoefs(Tn_arr, Tn_arr, Keddy_arr, Dcoef_dict, H0_dict, Hs_dict)
    tup = fill(-999., length(specieslist), num_layers)
    tdown = fill(-999., length(specieslist), num_layers)
    for (s, i) in zip(fullspecieslist, collect(1:length(fullspecieslist)))
        tup[i, :] .= fluxcoefs_all[s][2:end-1, 2]
        tdown[i, :] .= fluxcoefs_all[s][2:end-1, 1]
    end

    bc_dict_new = boundaryconditions(fluxcoefs_all, speciesbclist)
    tlower = permutedims(reduce(hcat, [bc_dict_new[sp][1,:] for sp in specieslist]))
    tupper = permutedims(reduce(hcat, [bc_dict_new[sp][2,:] for sp in specieslist]))

    # update to next timestep
    nthis = next_timestep(nstart, nthis, inactive, activespecies, inactivespecies, Jrates, Tn_arr, tup, tdown,
                          tlower, tupper, dt, ppt_vals, rel_tol)
    nthismat = reshape(nthis, (length(activespecies), num_layers))

    # write found values out to n_current
    for s in 1:length(activespecies)
        for ia in 1:num_layers
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

user_input_folder_name = input("Enter a name for the results folder or press enter to use default: ")
const sim_folder_name = user_input_folder_name == "" ? FNext : user_input_folder_name
# Set up the folder if it doesn't exist
create_folder(sim_folder_name, results_dir)

# Case where this file will be used to converge an atmosphere of a different extent
make_new_alt_grid = input("Would you like to use the script to converge a new atmosphere of a different extent? (y/n): ")
while make_new_alt_grid != "y" && make_new_alt_grid != "n"
    println("Bad entry!")
    global make_new_alt_grid = input("Would you like to use the script to converge a new atmosphere of a different extent? (y/n): ")
end
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
end

# Various other needful things ===============================================================================

# Densities of inactive species, which won't change by definition
const inactive = flatten_atm(n_current, inactivespecies)
const Dcoef_arr_template = 0. .* similar(alt)

# Set up the timesteps
const mindt = dt_min_and_max["neutrals"][1]
const maxdt = dt_min_and_max["neutrals"][2]

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
#                       TEMPERATURE/PRESSURE PROFILES                          #
################################################################################

#If changes to the temperature are needed, they should be made here
if args[1]=="temp"
    const T_surf = args[2]
    const T_tropo = args[3]
    const T_exo = args[4]
    Temp(z::Float64) = Tpiecewise(z, args[2], args[3], args[4])
    Temp_keepSVP(z::Float64) = Tpiecewise(z, meanTs, meanTt, meanTe) # for testing temp without changing SVP.
else
    const T_surf = meanTs
    const T_tropo = meanTt
    const T_exo = meanTe
    Temp(z::Float64) = Tpiecewise(z, meanTs, meanTt, meanTe)
end

const Tn_arr = [Temp(a) for a in alt];
const controltemps = [T_surf, T_tropo, T_exo]

# Species-specific scale heights - has to be done here once the control temps are set
Hs_dict = Dict{Symbol, Vector{Float64}}([sp=>map(z->scaleH(z, sp, controltemps), alt) for sp in fullspecieslist])

plot_temp_prof([Temp(a) for a in alt], savepath=results_dir*sim_folder_name)


################################################################################
#                               WATER PROFILES                                 #
################################################################################
println("Setting up the water profile...")
# set SVP to be fixed or variable with temperature
if args[1] == "temp"
    fix_SVP = true
    const H2Osat = map(x->Psat(x), map(Temp_keepSVP, alt)) # for holding SVP fixed
    const HDOsat = map(x->Psat_HDO(x), map(Temp_keepSVP, alt))  # for holding SVP fixed
else
    fix_SVP = false
    const H2Osat = map(x->Psat(x), map(Temp_n, alt)) # array in #/cm^3 by altitude
    const HDOsat = map(x->Psat_HDO(x), map(Temp_n, alt))
end

# H2O Water Profile ============================================================
surface_watersat = Dict("H2O"=>H2Osat[1], "HDO"=>HDOsat[1])
H2Osatfrac = H2Osat./map(z->n_tot(n_current, z), alt)  # get SVP as fraction of total atmo
# set H2O SVP fraction to minimum for all alts above first time min is reached
H2Oinitfrac = H2Osatfrac[1:something(findfirst(isequal(minimum(H2Osatfrac)), H2Osatfrac), 0)]
H2Oinitfrac = [H2Oinitfrac;   # ensures no supersaturation
               fill(minimum(H2Osatfrac), num_layers-length(H2Oinitfrac))]

# make profile constant in the lower atmosphere (well-mixed).
# when doing water experiments, the temperature profile is such that the minimum
# in the SVP curve occurs at about 55 km alt. This was found by manual tweaking.
# thus we need to only set the lower atmo mixing ratio below that point, or 
# there will be a little spike in the water profile.
if args[1] == "water"
    H2Oinitfrac[findall(x->x<hygropause_alt, alt)] .= args[2]
    MR = args[2] # mixing ratio 
else
    MR = MR_mean_water
    H2Oinitfrac[findall(x->x<hygropause_alt, alt)] .= MR # 10 pr μm
end

for i in [1:length(H2Oinitfrac);]
    H2Oinitfrac[i] = H2Oinitfrac[i] < H2Osatfrac[i+1] ? H2Oinitfrac[i] : H2Osatfrac[i+1]
end

# set the water profiles =======================================================
n_current[:H2O] = H2Oinitfrac.*map(z->n_tot(n_current, z), non_bdy_layers)
n_current[:HDO] = 2 * DH * n_current[:H2O] # This should be the correct way to set the HDO profile.
# n_current[:HDO] = HDOinitfrac.*map(z->n_tot(n_current, z), non_bdy_layers) # OLD WAY that is wrong.

# ADD EXCESS WATER AS FOR DUST STORMS.
dust_storm_on = "no"#input("Add excess water to simulate a dust storm? (y/n)")
if dust_storm_on == "y" || dust_storm_on == "yes"
    H2Oppm = 1e-6*map(x->250 .* exp(-((x-42)/12.5)^2), non_bdy_layers/1e5) + H2Oinitfrac  # 250 ppm at 42 km (peak)
    HDOppm = 1e-6*map(x->0.350 .* exp(-((x-38)/12.5)^2), non_bdy_layers/1e5) + HDOinitfrac  # 350 ppb at 38 km (peak)
    n_current[:H2O] = H2Oppm .* map(z->n_tot(n_current, z), non_bdy_layers)
    n_current[:HDO] = HDOppm .* map(z->n_tot(n_current, z), non_bdy_layers)
end

# We still have to calculate the HDO initial fraction in order to calculate the pr um 
# and make water plots.
HDOinitfrac = n_current[:HDO] ./ map(z->n_tot(n_current, z), non_bdy_layers)  

# Compute total water column for logging and checking that we did things right.
# H2O #/cm^3 (whole atmosphere) = sum(MR * n_tot) for each alt. Then convert to 
# pr μm: (H2O #/cm^3) * cm * (mol/#) * (H2O g/mol) * (1 cm^3/g) * (10^4 μm/cm)
# where the lone cm is a slice of the atmosphere of thickness dz, #/mol=6.023e23, 
# H2O or HDO g/mol = 18 or 19, cm^3/g = 1 or 19/18 for H2O or HDO.
# written as conversion factors for clarity.
H2O_per_cc = sum([MR; H2Oinitfrac] .* map(z->n_tot(n_current, z), alt[1:end-1]))
HDO_per_cc = sum([MR*DH; HDOinitfrac] .* map(z->n_tot(n_current, z), alt[1:end-1]))
H2Oprum = (H2O_per_cc * dz) * (18/1) * (1/6.02e23) * (1/1) * (1e4/1)
HDOprum = (HDO_per_cc * dz) * (19/1) * (1/6.02e23) * (19/18) * (1e4/1)

################################################################################
#                             BOUNDARY CONDITIONS                              #
################################################################################

H_effusion_velocity = effusion_velocity(Temp(zmax), 1.0, zmax)
H2_effusion_velocity = effusion_velocity(Temp(zmax), 2.0, zmax)
D_effusion_velocity = effusion_velocity(Temp(zmax), 2.0, zmax)
HD_effusion_velocity = effusion_velocity(Temp(zmax), 3.0, zmax)

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

    the rates at each altitude can be computed using the reaction network
    already in place, plus additional equations describing the transport
    to and from the cells above and below:
=#
const upeqns = [Any[Any[[s], [Symbol(string(s)*"_above")],Symbol("t"*string(s)*"_up")],
                    Any[[Symbol(string(s)*"_above")],[s],Symbol("t"*string(s)*"_above_down")]]
                    for s in specieslist]

const downeqns = [Any[Any[[s], [Symbol(string(s)*"_below")],Symbol("t"*string(s)*"_down")],
                      Any[[Symbol(string(s)*"_below")],[s],Symbol("t"*string(s)*"_below_up")]]
                      for s in specieslist]

const local_transport_rates = [[[Symbol("t"*string(s)*"_up") for s in specieslist]
                                [Symbol("t"*string(s)*"_down") for s in specieslist]
                                [Symbol("t"*string(s)*"_above_down") for s in specieslist]
                                [Symbol("t"*string(s)*"_below_up") for s in specieslist]]...;]

const transportnet = [[upeqns...;]; [downeqns...;]]

# define names for all the species active in the coupled rates:
const active_above = [Symbol(string(s)*"_above") for s in activespecies]
const active_below = [Symbol(string(s)*"_below") for s in activespecies]

# obtain the rates and jacobian for each altitude
const rates_local = Expr(:vcat, map(x->getrate(reactionnet, transportnet, x),activespecies)...);
const chemJ_local = chemical_jacobian(reactionnet, transportnet, activespecies, activespecies);
const chemJ_above = chemical_jacobian(reactionnet, transportnet, activespecies, active_above);
const chemJ_below = chemical_jacobian(reactionnet, transportnet, activespecies, active_below);

# TODO: Turn these notes into a test that looks for "+()" in rates_local instead of debugging lines. 
# -----------------------------------------------------------------------------------------------
# These lines are useful for troubleshooting if you're getting weird errors with ratefn.
# "no method matching +()" means the chemical system is unbalanced and you've got a
# species that has production but no consumption, or vice versa.
# other errors may occur. Uncomment these 3 lines to inspect what goes into rates_local.
# println("The contents of rates_local: ")
# println(rates_local)
# println("That was rates_local")

const arglist_local = [activespecies; active_above; active_below; inactivespecies;
                       Jratelist; :T; local_transport_rates; :dt]

const arglist_local_typed=[:($s::Float64) for s in arglist_local]

# This expression which is evaluated below enables a more accurate assessment of M
# by calculating at the time of being called rather than only at each timestep.
const Mexpr = Expr(:call, :+, fullspecieslist...)

# NOTE: These functions within @eval cannot be moved. Do not move them.
@eval begin
    function ratefn_local($(arglist_local_typed[1:end-1]...))

        # M is calculated here to ensure that the right number of particles is used.
        # It is for only the altitude at which this function was called 
        # (i.e. all the arguments to the function, when it's called, are for only one altitude)
        M = $Mexpr
        $rates_local # evaluates the rates_local expression
    end
end


@eval begin
    function chemJmat_local($(arglist_local_typed...))

        M = $Mexpr
        
        localchemJi = $(chemJ_local[1])
        localchemJj = $(chemJ_local[2])
        localchemJval = -dt*$(Expr(:vcat, chemJ_local[3]...)) 

        abovechemJi = $(chemJ_above[1])
        abovechemJj = $(chemJ_above[2])
        abovechemJval = -dt*$(Expr(:vcat, chemJ_above[3]...))

        belowchemJi = $(chemJ_below[1])
        belowchemJj = $(chemJ_below[2])
        belowchemJval = -dt*$(Expr(:vcat, chemJ_below[3]...))

        ((localchemJi, localchemJj, localchemJval),
         (abovechemJi, abovechemJj, abovechemJval),
         (belowchemJi, belowchemJj, belowchemJval))
    end
end


################################################################################
#                         PHOTOCHEMICAL CROSS SECTIONS                         #
################################################################################
println("Populating cross section dictionary...")
const crosssection = populate_xsect_dict(controltemps, extended_network=false)


# Solar Input ------------------------------------------------------------------
solarflux=readdlm(research_dir*solarfile,'\t', Float64, comments=true, comment_char='#')[1:2000,:]
solarflux[:,2] = solarflux[:,2]/2  # To roughly put everything at an SZA=60° (from a Kras comment)

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
#                                 LOGGING                                      #
################################################################################

println("Creating the simulation log file...")

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
f = open(results_dir*sim_folder_name*"/simulation_params_"*FNext*".txt", "w")
write(f, "$(args[1]) experiment: \n")
write(f, input_string)
write(f, "\n")

# Mean temperatures
write(f, "Mean temperatures used:\n")
write(f, "Surface: $(meanTs) K, Tropopause: $(meanTt) K, Exobase: $(meanTe) K\n\n")

# which species are turned on
write(f, "All species: $(join([string(i) for i in fullspecieslist], ", "))\n")
write(f, "No-chem species: $(join([string(i) for i in nochemspecies], ", "))\n")
write(f, "No-transport species: $(join([string(i) for i in notransportspecies], ", "))\n")
write(f, "Active species: $(join(sort([string(i) for i in activespecies]), ", "))\n\n")

# Solar cycle and SVP stuff
write(f, "\nSVP fixed: $(fix_SVP)\n\n")
write(f, "Solar cycle status: solar $(cycle)\n\n")

# Water profile information
write(f, "Water profile information: \n")
write(f, "Total H2O col: $(H2O_per_cc*2e5)\n")
write(f, "Total HDO col: $(HDO_per_cc*2e5)\n")
write(f, "Total water col: $((H2O_per_cc + HDO_per_cc)*2e5)\n")
write(f, "H2O+HDO at surface: $((H2O_per_cc[1] + HDO_per_cc[1])*2e5)\n")
write(f, "Total H2O (pr μm): $(H2Oprum)\n")
write(f, "Total HDO (pr μm): $(HDOprum)\n")
write(f, "Total H2O+HDO, no enhancement: $(H2Oprum + HDOprum)\n\n")


# cross sections
write(f, "\nCROSS SECTIONS: \n")
for k in keys(xsect_dict)  # cross sections
    write(f, k*": "*join(xsect_dict[k], ", ")*"\n")
end
write(f, "\n")

# boundary conditions
write(f, "\nBOUNDARY CONDITIONS: \n")
write(f, "n: number density at surface, f: flux at top, v: velocity at top\n")
for k2 in keys(speciesbclist)
    bcstring = join([join(speciesbclist[k2][1, :], "="),
                     join(speciesbclist[k2][2, :], "=")], ", ")
    write(f, string(k2)*": "*bcstring*"\n")
end
write(f, "\n")

# timestep etc
write(f, "Initial atmosphere state: $(readfile)\n\n")

write(f, "Absolute error tolerance = 1 ppt\n")


################################################################################
#                             CONVERGENCE CODE                                 #
################################################################################

# Uncomment this line if you'd like to add an extra parcel to some species. You must specify the species.
# n_current[:D] = map(x->1e5*exp(-((x-184)/20)^2), non_bdy_layers/1e5) + n_current[:D]

# write initial atmospheric state ==============================================
write_ncurrent(n_current, results_dir*sim_folder_name*"/initial_state.h5")

# Plot initial water profile ===================================================
plot_water_profile(H2Oinitfrac, HDOinitfrac, n_current[:H2O], n_current[:HDO], results_dir*sim_folder_name, watersat=H2Osatfrac)

# Plot initial atmosphere condition  ===========================================
println("Plotting the initial condition...")
plot_atm(n_current, [fullspecieslist], results_dir*sim_folder_name*"/initial_atmosphere.png")


# Record setup time
t6 = time()
write(f, "Setup time $(format_sec_or_min(t6-t5))\n")

# do the convergence ===========================================================
img_savepath = results_dir*sim_folder_name*"/converged_atm_"*FNext*".png"
println("Beginning Convergence")
t7 = time()
@showprogress 0.1 "Converging over 10 My..." [timeupdate(t, results_dir*sim_folder_name, Dcoef_arr_template) for t in [10.0^(1.0*i) for i in -3:14]]
@showprogress 0.1 "Last convergence steps..." for i in 1:100
    # plot_atm(n_current, [fullspecieslist], img_savepath, t="1e14", iter=i) # Turned off to improve speed. Can be run in separate script.
    update!(n_current, 1e14, Dcoef_arr_template)
end
t8 = time()

write(f, "Simulation runtime $(format_sec_or_min(t8-t7))\n")

# Plot the final atmosphere
plot_atm(n_current, [fullspecieslist], img_savepath, t="final, converged")

# write out the new converged file to matching folder.
towrite = results_dir*FNext*"/converged_"*FNext*".h5"
write_ncurrent(n_current, towrite)
println("Wrote $(towrite)")

# save the figure
savefig(results_dir*FNext*"/converged_"*FNext*".png", bbox_inches="tight")
println("Saved figure to same folder")

println("ALERT: Finished")
println()

t9 = time()

write(f, "Total runtime $(format_sec_or_min(t9-t1))\n")
close(f)
