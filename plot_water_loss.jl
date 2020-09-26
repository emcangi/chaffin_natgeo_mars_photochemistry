################################################################################
# plot_water_loss.jl
# TYPE: (2) Analysis - required
# DESCRIPTION: Estimates water lost from Mars for a given array of fractionation 
# factors and plots it, along with the Rayleigh distillation equation.
#
# Eryn Cangi
# Created 2019
# Last edited: 21 July 2020
# Currently tested for Julia: 1.4.1
################################################################################

using Analysis
using PyPlot
using PyCall
using JLD

include("PARAMETERS.jl")

function h2olost(cur_h2o, DHnow=5.5, DH0=1, f=0.001)
    #=
    Calculates water lost, given:

    cur_h2o: current water inventory. Can be in whatever unit you like, m GEL common
    DHnow: present day atmospheric D/H in water
    DH0: ancient Mars atmospheric D/H in water
    f: fractionation factor 
    =#
    return cur_h2o * ((DHnow/DH0)^(1/(1-f)) - 1)
end

function required_escape_rate(W_lost)
    #=
    Calculates the required escape rate of H atoms, assuming a constant escape
    rate and a linear change with time. NOTE: this assumes all water lost 
    (W_lost) is in the form H2O.

    1.419e17 is 4.5 billion years in seconds.
    =#
    H_lost = GEL_to_molecule(W_lost, "H")
    flux = H_lost / 1.419e17 # units #/cm^2s 
end

function plot_water_loss(f_therm, f_both, f_therm_range=nothing, f_both_range=nothing)
    #=
    Plot water loss vs. current water inventory estimate.
    Accounts for fractionation factor results and also plots a line
    according to the D/H ratio used. 

    f_therm: mean fractionation factor for only thermal escape
    f_both: mean fractionation factor for both thermal + nonthermal
    f_therm_range: min and max fractionation factor for thermal escape
    f_both_range: min and max fractionation factor for all types of escape
    =#
    # set up the domain 
    cur_h2o = range(20, stop=30, length=11)

    # make plot look nice
    fig = figure(figsize=(7,5))
    rcParams = PyCall.PyDict(matplotlib."rcParams")
    rcParams["font.family"] = "sans-serif"
    rcParams["font.sans-serif"] = "Louis George Caf?"
    rcParams["font.monospace"] = "FreeMono"
    rcParams["font.size"] = 18
    rcParams["axes.labelsize"] = 20
    rcParams["xtick.labelsize"] = 18
    rcParams["ytick.labelsize"] = 18
    ax = gca()
    plot_bg(ax)
    xlabel("Current exchangeable water (m GEL)")
    ylabel("Total water loss (m GEL)")

    # plot 
    blues = get_grad_colors(2, "Blues", strt=0.5)

    # calculate and plot the loss for mean f
    loss_t = h2olost(cur_h2o, 5.5, 1.275, f_therm)
    loss_b = h2olost(cur_h2o, 5.5, 1.275, f_both)

    ax.plot(cur_h2o, loss_t, zorder=2, color=blues[1, :], marker="D", label="f=0.002")
    ax.plot(cur_h2o, loss_b, zorder=2, color=blues[2, :], marker="D", label="f=0.06")

    if f_therm_range != nothing
        therm_min_loss = h2olost(cur_h2o, 5.5, 1.275, f_therm_range[1])
        therm_max_loss = h2olost(cur_h2o, 5.5, 1.275, f_therm_range[2])
        ax.fill_between(cur_h2o, therm_min_loss, therm_max_loss, zorder=2,
                        alpha=0.4, color=blues[1, :])
    end

    if f_both_range != nothing
        both_min_loss = h2olost(cur_h2o, 5.5, 1.275, f_both_range[1])
        both_max_loss = h2olost(cur_h2o, 5.5, 1.275, f_both_range[2])
        ax.fill_between(cur_h2o, both_min_loss, both_max_loss, zorder=2,
                        alpha=0.4, color=blues[2, :])
    end

    ax.text(21, 83, "thermal+non-thermal escape", rotation=24, color=blues[2, :])
    ax.text(22, 68, "thermal escape only", rotation=23, color=blues[1, :])

    ax.legend()

    savefig(results_dir*"AllResultPlots/h2oloss.png", bbox_inches="tight", dpi=300)
    savefig(results_dir*"MainCases/h2oloss.png", bbox_inches="tight", dpi=300)

    println("Thermal escape water loss (based on f for mean temperature profile):")
    println(Array(loss_t))
    println("All escape water loss (based on f for mean temperature profile): ")
    println(Array(loss_b))

    println()

    println("Thermal escape minimum water loss (f minima):")
    println(Array(therm_min_loss))
    println("Thermal escape maximum water loss (f maxima):")
    println(Array(therm_max_loss))

    println()

    println("All escape minimum water loss (f minima):")
    println(Array(both_min_loss))
    println("All escape maximum water loss (f maxima):")
    println(Array(both_max_loss))
end

function plot_water_loss_vs_temp()
    #=
    Plot water loss as a result of fractionation factor results and also a line
    according to the D/H ratio used. 

    In this plot, the x-axis is exobase temperature. It was created to convince
    a reviewer that for a specific current water inventory estimate, 
    that increased temps => increased f => increased water loss.

    f_therm: mean fractionation factor for only thermal escape
    f_both: mean fractionation factor for both thermal + nonthermal
    f_therm_range: min and max fractionation factor for thermal escape
    f_both_range: min and max fractionation factor for all types of escape
    =#
    # set up the domain 
    cur_h2o = 25 # a random value in the range

    # make plot look nice
    fig = figure(figsize=(7,5))
    rcParams["font.family"] = "sans-serif"
    rcParams["font.sans-serif"] = "Louis George Caf?"
    rcParams["font.monospace"] = "FreeMono"
    rcParams["font.size"] = 18
    rcParams["axes.labelsize"] = 20
    rcParams["xtick.labelsize"] = 18
    rcParams["ytick.labelsize"] = 18
    ax = gca()
    plot_bg(ax)
    xlabel("Temperature (K)")
    ylabel("Total water loss (m GEL)")

    # plot 
    blues = get_grad_colors(2, "Blues", strt=0.5)

    # calculate and plot the loss for mean f
    T_abs_dicts = load(results_dir*det_cases_dir*"T_abs_data.jld")
    f_by_temp = T_abs_dicts["T_data_abs"]["exobase"]["f"]
    println(f_by_temp)

    loss_t = [h2olost(cur_h2o, 5.5, 1.275, f) for f in f_by_temp]

    ax.plot(150:25:350, loss_t, zorder=2, color=blues[1, :], marker="D")
    ax.text(150, 120, "Water loss via thermal escape\nfor a current inventory of 25 m GEL", color=medgray, va="top")
    savefig(results_dir*"AllResultPlots/h2oloss-bytemp.png", bbox_inches="tight", dpi=300)
    savefig(results_dir*"MainCases/h2oloss-bytemp.png", bbox_inches="tight", dpi=300)
end

function plot_simple_water_loss(cur_h2o=25)
    #=
    Plot just the Rayleigh distillation equation, no fractionation factor results.
    For demonstration purposes. Bruce asked for this. Not in the paper.

    cur_h2o: assumed water inventory in m GEL.
    =#
    farray = range(0.000005, stop=1, length=100)

    l1 = [h2olost(cur_h2o, 5.5, 1.275, f) for f in farray]

    # make plot
    fig = figure(figsize=(10, 7))

    # style.use("default")
    rc("text", usetex=false)
    rcParams = PyCall.PyDict(matplotlib."rcParams")
    rcParams["font.family"] = "sans-serif"
    rcParams["font.sans-serif"] = "Louis George Caf?"
    rcParams["font.monospace"] = "FreeMono"
    rcParams["font.size"] = 26
    rcParams["axes.labelsize"] = 28
    rcParams["xtick.labelsize"] = 26
    rcParams["ytick.labelsize"] = 26

    ax = gca()
    plot_bg(ax)

    xlabel(r"Fractionation factor $f$")
    ylabel("Watr lost (m GEL) in 4.5 Ga")

    # plot things
    ax.semilogx(farray, l1, zorder=2, color="cornflowerblue")
    # legend(fontsize=24)
    text(0.000005, 50, "Water loss constant \n"+L" with $f$ when $f\ll0.05$")
    text(0.0006, 200, "Water loss depends \nstrongly on " + L"$f$ for $f\gtrapprox 0.05$")
    ylim(0, 240)
    savefig(results_dir*"AllResultPlots/h2oloss_simple.png", bbox_inches="tight")
end

println("ALERT: f values are hard-coded in this file. Make sure they're correct."*
        " They must be filled in manually, using output from plot_f_results.jl")
f_mean_thermal = 0.001899480829925377
f_therm_range = [3.26998e-5, 0.0171938]
f_mean_both = 0.05964736746762264
f_both_range = [0.0354146, 0.0971759]

println("Echoing f inputs for clarity:")
println("f_mean_thermal: $(f_mean_thermal)")
println("f_therm_range: $(f_therm_range)")
println("f_mean_both: $(f_mean_both)")
println("f_both_range: $(f_both_range)")

plot_water_loss(f_mean_thermal, f_mean_both, f_therm_range, f_both_range)

plot_water_loss_vs_temp()