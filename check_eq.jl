################################################################################
# check_eq.jl
# TYPE: (2a) Analysis - required
# DESCRIPTION: Check if the atmosphere is in equilibrium by looking for
# Φ_O = 0.5Φ_{H+D}.
#
# Eryn Cangi
# Created: 11 September 2018
# Last edited: 20 July 2020
# Currently tested for Julia: 1.4.1
################################################################################
using HDF5
using Analysis

include("PARAMETERS.jl")

# Experiment type and loading the file =========================================
# get arguments to use for accessing files 
argarray = Any[ARGS[i] for i in 1:1:length(ARGS)] 

# the options are "temp", "water", "dh"
if argarray[1] == "temp"
    FNext = "temp_$(argarray[2])_$(argarray[3])_$(argarray[4])"
    temparr = [parse(Float64, a) for a in argarray[2:end]]
else
    FNext = "$(argarray[1])_$(argarray[2])"
    temparr = [meanTs, meanTt, meanTe]
end

# load readfile and altitude array
main_results_path = results_dir * main_cases_dir
detail_results_path = results_dir * det_cases_dir
sim_rel_path = FNext * "/converged_" * FNext * ".h5"

# Set up the file path, automatically look for the simulation without external 
# specification of where it might be. 
# there may be duplicates in both folders but that's fine, it will do this once.
if isfile(main_results_path * sim_rel_path)
    rf = main_results_path * sim_rel_path
elseif isfile(detail_results_path * sim_rel_path)
    rf = detail_results_path * sim_rel_path  
else
    throw("Invalid: that simulation does not exist")
end


# Calculate the flux ratio =====================================================
Of = 1.2e8
Hf = get_flux(:H, rf, temparr, therm_only=true)
HDf = Hf + get_flux(:D, rf, temparr, therm_only=true)


if round(HDf/Of) == 2
    println("Simulation: $(FNext)")
    println("Equilibrium: YES")
    println()
else
    println("Simulation: $(FNext)")
    println("Equilibrium: NO")
    println("Ratio: ", HDf/Of)
    println()
end


