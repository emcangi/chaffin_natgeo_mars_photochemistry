################################################################################
# reproductions_plot.jl
# TYPE: (2) Analysis - required
# DESCRIPTION: Makes a plot of reproductions of past studies. Is a separate 
# file for simplicity and because it allows more easy access to modify the way
# I calculate the reproduction of Krasnopolsky 2002 (which can be done either
# as the models were run, or by including his nonthermal escape velocities)
#
# Eryn Cangi
# Created 31 May 2019
# Last edited: 21 July 2020
# Currently tested for Julia: 1.4.1
################################################################################

using PyPlot
using PyCall
using HDF5
using Printf
using LaTeXStrings
using Analysis

include("PARAMETERS.jl")

global rcParams = PyCall.PyDict(matplotlib."rcParams")

# Main =========================================================================

function make_reproduction_plot()
    #=
    Plots our attempts at reproducing past studies against those studies.
    =#
    Oflux = 1.2e8 
    filemin = results_dir*"Replications/Kras2002-rep/temp_213_125_200/converged_temp_213_125_200.h5"
    filemean = results_dir*"Replications/Kras2002-rep/temp_213_125_270/converged_temp_213_125_270.h5"
    filemax = results_dir*"Replications/Kras2002-rep/temp_213_125_350/converged_temp_213_125_350.h5"


    # Calculate f values -------------------------------------------------------
    f_k02min_therm = calculate_f(filemin, "thermal", [213.0, 125.0, 200.0], reprod=true)
    f_k02mean_therm = calculate_f(filemean, "thermal", [213.0, 125.0, 270.0], reprod=true)
    f_k02max_therm = calculate_f(filemax, "thermal", [213.0, 125.0, 350.0], reprod=true)

    f_k02min_both = calculate_f(filemin, "both", [213.0, 125.0, 200.0], reprod=true)
    f_k02mean_both = calculate_f(filemean, "both", [213.0, 125.0, 270.0], reprod=true)
    f_k02max_both = calculate_f(filemax, "both", [213.0, 125.0, 350.0], reprod=true)

    println(f_k02min_therm)

    flist_repro_therm = [["past", "Replication", 0.26],
               ["past", "Yung+1988", 0.32],
               ["past", "Replication", round(f_k02min_therm, sigdigits=2)],
               ["past", L"Kras. 2002 $\odot$ min", 0.055],
               ["past", "Replication", round(f_k02mean_therm, sigdigits=2)],
               ["past", L"Kras. 2002 $\odot$ mean", 0.082],
               ["past", "Replication", round(f_k02max_therm, sigdigits=2)],
               ["past", L"Kras. 2002 $\odot$ max", 0.167]]

    flist_repro_both = [["past", "Replication", 0.26],
               ["past", "Yung+1988", 0.32],
               ["past", "Replication", round(f_k02min_both, sigdigits=2)],
               ["past", L"Kras. 2002 $\odot$ min", 0.055],
               ["past", "Replication", round(f_k02mean_both, sigdigits=2)],
               ["past", L"Kras. 2002 $\odot$ mean", 0.082],
               ["past", "Replication", round(f_k02max_both, sigdigits=2)],
               ["past", L"Kras. 2002 $\odot$ max", 0.167]]

    # these lines reshape the above arrays to the same shape, but a different type.
    # I don't know. Julia is really annoying with array shapes and types.
    # Basically these lines are witchcraft that I somehow figured out once.
    flist_repro_therm = permutedims(reshape(hcat(flist_repro_therm...), 
                                            (length(flist_repro_therm[1]), 
                                             length(flist_repro_therm))))

    flist_repro_both = permutedims(reshape(hcat(flist_repro_both...), 
                                        (length(flist_repro_both[1]), 
                                         length(flist_repro_both))))

    # REPRODUCTION PLOT ========================================================
    fig, ax = subplots(1, 2, sharey=true, sharex=true, figsize=(14,7))
    for a in ax
        plot_bg(a)
    end
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 16
    rcParams["axes.labelsize"]= 20
    rcParams["xtick.labelsize"] = 18
    rcParams["ytick.labelsize"] = 18
    carray_repros = ["#74c476", "#74c476", "#a3bfdc", "#a3bfdc", "#8c96c6", 
                 "#8c96c6", "#8c6bb1", "#8c6bb1"]

    ax[1].barh(collect(1:1:size(flist_repro_therm)[1]), flist_repro_therm[:, 3], 
               0.9, color=carray_repros, zorder=10)
    ax[1].set_xscale("log")
    ax[1].set_xlabel("Fractionation Factor", fontsize=20)
    ax[1].set_yticks(collect(1:1:size(flist_repro_therm)[1]))
    ax[1].set_yticklabels(flist_repro_therm[:, 2], fontsize=18)
    ax[1].tick_params(labelsize=18)
    ax[1].set_title("Thermal escape only", fontsize=22)
    ax[1].set_xticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1])
    ax[1].set_xlim(1e-4, 1)
    # next block just lists the value on the bar
    for (i, v) in enumerate(flist_repro_therm[:, 3])
        if v <= 2e-5
            m = 4
        else
            m = -0.05
        end

        printvalstyle = "$(v)"
        ax[1].text(v+m*v, i-0.1, printvalstyle, color="black", zorder=16,
                ha="right")
    end

    ax[2].barh(collect(1:1:size(flist_repro_both)[1]), flist_repro_both[:, 3], 
               0.9, color=carray_repros, zorder=10)
    ax[2].set_xscale("log")
    ax[2].set_xlabel("Fractionation Factor", fontsize=20)
    ax[2].set_yticks(collect(1:1:size(flist_repro_both)[1]))
    ax[2].set_yticklabels(flist_repro_both[:, 2], fontsize=18)
    ax[2].tick_params(labelsize=18)
    ax[2].set_title("Thermal + non-thermal escape", fontsize=22)
    ax[2].set_xticks([1e-4, 1e-3, 1e-2, 1e-1, 1])
    for (i, v) in enumerate(flist_repro_both[:, 3])
        if v <= 2e-5
            m = 4
        else
            m = -0.05
        end

        printvalstyle = "$(v)"
        ax[2].text(v+m*v, i-0.1, printvalstyle, color="black", zorder=16,
                ha="right")
    end

    savefig(results_dir*"ALL STUDY PLOTS/f-reproductions-plot.png", bbox_inches="tight")
end

# ==============================================================================

make_reproduction_plot()