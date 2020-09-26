################################################################################
# make_tradeoff_plots.jl
# TYPE: (2) Analysis - required
# DESCRIPTION: make results plots:
# Main plot of f results, and comparison with other studies
# H- and D-bearing species concentrations relative to stnadard atmosphere
# Model output comparison to data for O2, CO, H2, O3 
# etc.

# Eryn Cangi
# Created 5 April 2019
# Last edited: 20 July 2020
# Currently tested for Julia: 1.4.1
################################################################################
using PyPlot
using HDF5
using LaTeXStrings
using PyCall
using PlotUtils
using JLD
using Analysis
using DataFrames

include("PARAMETERS.jl")

# fundamental constants ========================================================
global tvals = Dict("surface"=>[150.0, 160.0, 170.0, 180.0, 190.0, 200.0, 210.0, 220.0,
                                230.0, 240.0, 250.0, 260.0, 270.0],
                    "tropopause"=>[100.0, 110.0, 120.0, 130.0, #70.0, 80.0, 90.0,
                                   140.0, 150.0, 160.0],
                    "exobase"=>[150.0, 175.0, 200.0, 225.0, 250.0, 275.0,
                                300.0, 325.0, 350.0])
global Oflux_vals = ["8e7", "9e7", "1.0e8", "1.1e8", "1.2e8", "1.3e8",
                            "1.4e8", "1.5e8", "1.6e8"]

global watervals = [1, 10, 25, 50, 100]#, 150]
global watervals_str = ["1.32e-5", "1.38e-4", "3.63e-4", "8.05e-4", "2.76e-3"]#, "1.358e-2"]

# nominal value plot location or index, 1-indexed for Julia
global nom_i_julia = Dict("exobase"=>meanTe, "tropopause"=>meanTt, "surface"=>meanTs,
                   "O flux"=>findfirst(isequal("1.2e8"), Oflux_vals),
                   "water"=>10)
# Passing the values to PyPlot requires 0-indexing so here it is:
global nom_i_py = Dict("exobase"=>meanTe, "tropopause"=>meanTt, "surface"=>meanTs,
                "O flux"=>findfirst(isequal("1.2e8"), Oflux_vals)-1,
                "water"=>10)

# set the data for comparison. Order:
# CO MR (Trainer+2019), O2 MR (Trainer+2019), H2 abundance (Kras&Feldman 2001),
# O3 μm-atm (Clancy 2016)
global data = [5.8e-4, 1.61e-3, 15, 1.18]  # CO MR, O2 MR, H2 ppm, O3 μm-atm
global s = [0.8e-4, 0.09e-3, 5, 0.7]     # sigmas (uncertainty) on each

# color settings to make plots have a nice design
Hcolor = "#5CA85C"
Hcolor_dark = "#152815"
Dcolor = "#1D4926"
Dcolor_dark = "#0C1D0F"
H2color = "#6888C5"
H2color_dark = "#223558"
HDcolor= "#384394"
HDcolor_dark = "#161B3B"
Hflux_color = "#DE9323"
Dflux_color = "#D74950"
COcolor = "#B684A1"
O2color = "#584268"
O3color = "#256E51"
f_main_color = "#574D9D"
f_sec_color = "#588C54"
sz = 10

# Functions ====================================================================

# Some utility functions -------------------------------------------------------

function normalize(arr, base_i)
    normed_arr = arr ./ arr[base_i]
    return normed_arr
end

function normalize_val(arr, normval)
    normed_arr = arr ./ normval
    return normed_arr
end

# The main results plotting function -------------------------------------------
function plot_results_caltech_together(base)
    #=
    Plots the thermal and thermal+nonthermal results together in the same plot.

    base: the folder of results
    =#
    n_alt_index=Dict([z=>clamp((i-1),1, length(alt)-2) for (i, z) in enumerate(alt)])
    Oflux = 1.2e8  # whatever Oflux is, in format "1.2e8" cm^-2 s^-1
    

    # get high and low f in surfaces
    f_surf = Array{Float64}(undef, 3, 4) #array to store f
    f_surf[:, 1] = [lowTs, meanTs, hiTs]
    i = 1
    for Ts in f_surf[:, 1]
        # construct file name
        F = "temp_$(Int(Ts))_$(Int(meanTt))_$(Int(meanTe))"
        filetouse = base*F*"/converged_"*F*".h5"

        # calculate f
        f_surf[i, 2] = calculate_f(filetouse, "thermal", [Ts, meanTt, meanTe])
        f_surf[i, 3] = calculate_f(filetouse, "both", [Ts, meanTt, meanTe])
        f_surf[i, 4] = calculate_f(filetouse, "nonthermal", [Ts, meanTt, meanTe])
        i += 1
    end

    # get high and low f in tropopause
    f_tropo = Array{Float64}(undef, 3, 4) #array to store f
    f_tropo[:, 1] = [lowTt, meanTt, hiTt]
    i = 1
    for Tt in f_tropo[:, 1]
        # construct file name
        F = "temp_$(Int(meanTs))_$(Int(Tt))_$(Int(meanTe))"
        filetouse = base*F*"/converged_"*F*".h5"

        # calculate f
        f_tropo[i, 2] = calculate_f(filetouse, "thermal", [meanTs, Tt, meanTe])
        f_tropo[i, 3] = calculate_f(filetouse, "both", [meanTs, Tt, meanTe])
        f_tropo[i, 4] = calculate_f(filetouse, "nonthermal", [meanTs, Tt, meanTe])
        i += 1
    end

    # get high and low f in exobase
    f_exo = Array{Float64}(undef, 3, 4) #array to store f
    f_exo[:, 1] = [lowTe, meanTe, hiTe]
    i = 1
    for Te in f_exo[:, 1]
        # construct file name
        F = "temp_$(Int(meanTs))_$(Int(meanTt))_$(Int(Te))"
        filetouse = base*F*"/converged_"*F*".h5"

        # calculate f
        f_exo[i, 2] = calculate_f(filetouse, "thermal", [meanTs, meanTt, Te])
        f_exo[i, 3] = calculate_f(filetouse, "both", [meanTs, meanTt, Te])
        f_exo[i, 4] = calculate_f(filetouse, "nonthermal", [meanTs, meanTt, Te])
        i += 1
    end

    # get high and low f in water
    f_water = Array{Float64}(undef, 5, 4) #array to store f
    waterfolders = search_subfolders(base, r"water_\d.+")
    i = 1
    for w in waterfolders

        # get the experiment name
        waterexp = match(r"water_\d.+", w).match
        watermr = match(r"\d.+", waterexp).match
        # construct file name
        filetouse = w * "/converged_"*waterexp*".h5"

        # calculate f
        f_water[i, 1] = parse(Float64, watermr)
        f_water[i, 2] = calculate_f(filetouse, "thermal", [meanTs, meanTt, meanTe])
        f_water[i, 3] = calculate_f(filetouse, "both", [meanTs, meanTt, meanTe])
        f_water[i, 4] = calculate_f(filetouse, "nonthermal", [meanTs, meanTt, meanTe])
        i += 1
    end

    skey = L"T$_{surface}$"
    tkey = L"T$_{tropopause}$"
    ekey = L"T$_{exobase}$"

    # iterate through all the detailed files and make an array of fractionation factors by experiment
    # the mean case file is needed
    mn = base * "temp_$(meanTsint)_$(meanTtint)_$(meanTeint)/converged_temp_$(meanTsint)_$(meanTtint)_$(meanTeint).h5"

    # the f calculations for thermal 
    f_mean_thermal = calculate_f(mn, "thermal", [meanTs, meanTt, meanTe])
    f_mean_both = calculate_f(mn, "both", [meanTs, meanTt, meanTe])
    f_mean_nonthermal = calculate_f(mn, "nonthermal", [meanTs, meanTt, meanTe])

    f_thermal = DataFrame(Exp=["Surface", "Tropopause", "Exobase", "Water", ""],
                          Min=[minimum(f_surf[:, 2]), minimum(f_tropo[:, 2]), 
                               minimum(f_exo[:, 2]), minimum(f_water[:, 2]), 1],
                          Max=[maximum(f_surf[:, 2]), maximum(f_tropo[:, 2]), 
                               maximum(f_exo[:, 2]), maximum(f_water[:, 2]), 1])

    f_both = DataFrame(Exp=["Surface", "Tropopause", "Exobase", "Water", ""],
                       Min=[minimum(f_surf[:, 3]), minimum(f_tropo[:, 3]), 
                             minimum(f_exo[:, 3]), minimum(f_water[:, 3]), 1],
                       Max=[maximum(f_surf[:, 3]), maximum(f_tropo[:, 3]), 
                            maximum(f_exo[:, 3]), maximum(f_water[:, 3]), 1])

    f_nonthermal = DataFrame(Exp=["Surface", "Tropopause", "Exobase", "Water", ""],
                             Min=[minimum(f_surf[:, 4]), minimum(f_tropo[:, 4]), 
                                  minimum(f_exo[:, 4]), minimum(f_water[:, 4]), 1],
                             Max=[maximum(f_surf[:, 4]), maximum(f_tropo[:, 4]), 
                                  maximum(f_exo[:, 4]), maximum(f_water[:, 4]), 1])
 
    # Print the results for use in a table in the paper    
    println("thermal")
    println("f mean: $(f_mean_thermal)")
    println(f_thermal)
    println()
    println("thermal + nonthermal")
    println("f mean: $(f_mean_both)")
    println(f_both)
    println()
    println("nonthermal")
    println("f mean: $(f_mean_nonthermal)")
    println(f_nonthermal)

    toplot_others = Dict("Clarke+ 2019"=>[0.016, 0.047],
                   "Krasnopolsky 2002"=>[0.055, 0.082, 0.167],
                   "Krasnopolsky 2000"=>[0.016, 0.135],
                   "Yung+ 1988"=>[0.32]
                  )


    # PLOT =====================================================================
    fig, ax = subplots(figsize=(12,8))

    # STYLE SETTINGS -----------------------------------------------------------
    for side in ["left", "right"]
        ax.spines[side].set_visible(false)
    end
    ax.tick_params(which="both", labelleft=false, left=false, labeltop=true, top=true)

    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 18
    rcParams["axes.labelsize"]= 20
    rcParams["xtick.labelsize"] = 18
    rcParams["ytick.labelsize"] = 18

    # Which things to plot 
    plot_order_mine = ["Surface", "Tropopause", "Exobase", "Water", ""]
    plot_order_others = ["Clarke+ 2019", "Krasnopolsky 2002", 
                         "Krasnopolsky 2000", "Yung+ 1988"]
    dummyind_mine = [1, 2, 3, 4, 5]
    dummyind_others = [6, 7, 8, 9]

    myorange = "xkcd:faded orange"
    colors = ["#A35B24", myorange, "#FFAB6B", "#1493A3",
            #myorange, myorange, myorange, myorange, 
             "white", "#3c6843", "#682f51", "#99801d", "#485468"]
    scatcolors = ["#6B3C18", "#B8713B", "#E79B61", "#0e707c", "#fff4e0", 
                  "#3c6843", "#682f51", "#99801d", "#485468"]

    m = "D"
    lw = 5
    msize = 60 
    nameha = "right"
    nalign = 1

    # THIS STUDY ---------------------------------------------------------------
    ax.fill_between([1e-5, 1], [-0.5,-0.5], y2=[5.25, 5.25], color="#fff4e0", zorder=-10)
             
    # Plot actual values 
    for (d, e) in zip(dummyind_mine, plot_order_mine)
        x_thermal = [f_thermal[f_thermal.Exp.==e, :].Min[1], f_thermal[f_thermal.Exp.==e, :].Max[1]]
        x_both = [f_both[f_both.Exp.==e, :].Min[1], f_both[f_both.Exp.==e, :].Max[1]]
        y = fill!(similar(x_thermal), d)  # spaces them out in y axis

        # now actually plot them 
        ax.plot(x_thermal, y, linewidth=lw, color="white", zorder=4)#, solid_capstyle="round")
        ax.plot(x_thermal, y, linewidth=lw, color=colors[d], zorder=5, alpha=0.8)#, solid_capstyle="round")
        ax.scatter(x_thermal, y, color=scatcolors[d], zorder=6, marker=m, s=msize)
        ax.plot(x_both, y, linewidth=lw, color="white", zorder=4)
        ax.plot(x_both, y, linewidth=lw, color=colors[d], zorder=5, alpha=0.8)
        ax.scatter(x_both, y, color=scatcolors[d], zorder=6, marker=m, s=msize)
    end
    # mean cases
    ax.plot([f_mean_thermal, f_mean_thermal], [0, 4.5], color=myorange, zorder=-7)
    ax.scatter(f_mean_thermal, 0, color=myorange, zorder=5, marker=m, s=msize)
    ax.plot([f_mean_both, f_mean_both], [0, 4.5], color=myorange, zorder=-7) 
    ax.scatter(f_mean_both, 0, color=myorange, zorder=5, marker=m, s=msize)

    # text
    ax.text(1e-5, 5, "This study", color="black", fontsize=24, va="top")
    ax.text(2e-4, 4.75, "Thermal escape only")
    ax.text(6e-3, 4.75, "Thermal + non-thermal escape")
    ax.text(nalign, 0, "Standard atm.", color=myorange, ha=nameha)
    ax.text(f_mean_thermal-0.0003, 0, L"\overline{f}="*"$(round(f_mean_thermal, digits=3))", 
            color=myorange, ha="right", va="center")
    ax.text(f_mean_both-0.01, 0, L"\overline{f}="*"$(round(f_mean_both, digits=3))",
            color=myorange, ha="right", va="center")
    ax.text(nalign, 1, L"$\overline{T}_{surface} \pm 25$%", 
            color=colors[1], va="center", ha=nameha)
    ax.text(nalign, 2, L"$\overline{T}_{tropo} \pm 25$%", 
            color=colors[2], va="center", ha=nameha)
    ax.text(nalign, 3, L"$\overline{T}_{exobase} \pm 25$%", 
            color=colors[3], va="center", ha=nameha)
    ax.text(nalign, 4, L"Water 1-100 $\mathrm{\mu}$m", 
            color=colors[4], va="center", ha=nameha)

    # PAST STUDIES -------------------------------------------------------------
    for (d, e) in zip(dummyind_others, plot_order_others)
        x = toplot_others[e]
        y = fill!(similar(toplot_others[e]), d)
        c = colors[d]
        ax.scatter(x, y, linewidth=2, color=c, marker=m, s=msize, zorder=6)
    end

    ax.plot([1.15e-2, 6.8e-2], [6, 6], color=colors[6], linewidth=lw, linestyle=":", 
             alpha=0.4, zorder=4) 
    ax.plot(toplot_others["Krasnopolsky 2002"], fill!(similar(toplot_others["Krasnopolsky 2002"]), 7), 
            color=colors[7], linewidth=lw, alpha=0.4, zorder=4)
    ax.plot([1.02e-2, 2e-1], [8, 8], color=colors[8], linewidth=lw, linestyle=":", 
             alpha=0.4, zorder=4) 

    # text
    ax.text(1e-5, 9, "Past studies", color="black", fontsize=24, va="bottom")

    # Clarke+ 2019
    ax.text(nalign, 6, "Clarke+ 2019", color=colors[6], ha=nameha, va="center")
    ax.text(toplot_others["Clarke+ 2019"][1], 6.3, "Ls 280", color=colors[6], 
            ha="center", va="center", fontsize=14)
    ax.text(toplot_others["Clarke+ 2019"][2], 6.3, "Ls 318", color=colors[6], 
            ha="center", va="center", fontsize=14)
    ax.text(7e-2, 6, "?", color=colors[6], va="center", ha="center")
    ax.text(0.01, 6, "?", color=colors[6], va="center", ha="center")

    # Kras 2002
    ax.text(nalign, 7, "Kras. 2002", color=colors[7], ha=nameha, va="center")
    ax.text(toplot_others["Krasnopolsky 2002"][1], 7.3, "Solar min", 
            color=colors[7], ha="right", fontsize=14)
    ax.text(toplot_others["Krasnopolsky 2002"][2], 7.3, "mean", color=colors[7], 
            ha="center", fontsize=14)
    ax.text(toplot_others["Krasnopolsky 2002"][3], 7.3, "max", color=colors[7], 
            ha="center", fontsize=14)

    # Kras 2000
    ax.text(nalign, 8, "Kras. 2000", color=colors[8], ha=nameha, va="center")
    ax.text(toplot_others["Krasnopolsky 2000"][1], 8.3, "Diffusion profile 1", 
            color=colors[8], ha="center", fontsize=14)
    ax.text(toplot_others["Krasnopolsky 2000"][2], 8.3, "profile 2", 
            color=colors[8], ha="center", fontsize=14)
    ax.text(0.18, 8, "?", color=colors[8], va="center", ha="center")
    ax.text(1e-2, 8, "?", color=colors[8], va="center", ha="center")

    # Yung+ 1988
    ax.text(nalign, 9, "Yung\n1988", color=colors[9], ha=nameha, va="center")
    ax.text(0.22, 9, "Standard model (thermal escape only)", color=colors[9], ha="right", va="center", fontsize=14)
    
    # PLOT CONFIGURATION -------------------------------------------------------
    plt.margins(x=0.01, y=0.01)
    ax.set_xlabel(L"Fractionation factor $f$", fontsize=20)
    ax.set_xscale("log")
    ax.set_xticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0])
    ax.set_yticks(collect(1:1:10))
    ax.tick_params(which="both", labelsize=18)
    ax.set_xlim([1e-5, 1])
    # save to the two useful folders
    savefig(base*"f_results.png", bbox_inches="tight", dpi=300)
    savefig(results_dir*"AllResultPlots/f_results.png", bbox_inches="tight", dpi=300)
end

# Functions for analyzing the detailed cases -----------------------------------

function make_std_atmo_dict(abs_or_mr)
    #=
    Creates a dictionary object that contains various species concentrations, fluxes, and the
    fractionation factor for the standard atmosphere.

    abs_or_mr: whether to calculate the absolute abundances of a species or its mixing ratio

    returns: dictionary with the following keys and values:
        "H": absolute abundance of atomic H at the exobase
        "D": same, for atomic D
        "H2": same, for molecular hydrogen
        "HD": same, for hydrogen deuteride
        "H2MR": mixing ratio of H2 in the lower atmosphere (up to 80km) in ppm, used for comparing with krasnopolsky data
        "O2": absolute abundance of O2 at the surface
        "CO": same, but for CO
        "O3": ozone in #/cm^2 in the atmosphere (column abundance)
        "f": fractionation factor of the atmosphere
        "Hflux": escape flux of H due to loss of H, HD, and H2
        "Dflux": escape flux of D due to loss of D, HD
    =#

    # get the current array
    tfile = detailed_results_dir*"temp_$(meanTsint)_$(meanTtint)_$(meanTeint)/converged_temp_$(meanTsint)_$(meanTtint)_$(meanTeint).h5"
    ncur =  get_ncurrent(tfile)

    N0 = abs_or_mr == "abs" ? 1 : n_tot(ncur, 0, n_alt_index)
    Ntop = abs_or_mr == "abs" ? 1 : n_tot(ncur, zmax, n_alt_index)
    LA = collect(0e5:2e5:78e5)

    atmo_metric_dict = Dict("O2"=>[], "HD"=>[], "H2"=>[], "H2MR"=>[], "H"=>[], "D"=>[], "CO"=>[],
                     "Hflux"=>[], "Dflux"=>[], "f"=>[], "CO/O2"=>[], "O3"=>[],
                     "DH"=>[])
    DHprofs = Array{Any}(undef, 1, length(alt)-2)

    # Calculate the things we care about
    # H and D fluxes
    Hf, contrib_H = get_flux(:H, tfile, meantemps, therm_only=true)
    Df, contrib_D = get_flux(:D, tfile, meantemps, therm_only=true)
    append!(atmo_metric_dict["Hflux"], Hf)
    append!(atmo_metric_dict["Dflux"], Df)

    # Concentrations:
    # H at exobase
    append!(atmo_metric_dict["H"], ncur[:H][end]/Ntop)
    # D at exobase
    append!(atmo_metric_dict["D"], ncur[:D][end]/Ntop)
    # H2 at exobase
    append!(atmo_metric_dict["H2"], ncur[:H2][end]/Ntop)
    # HD at exobase
    append!(atmo_metric_dict["HD"], ncur[:HD][end]/Ntop)

    # D/H stuff:
    # D/H at 150 km
    append!(atmo_metric_dict["DH"], (ncur[:D][75]/ncur[:H][75])/(1.6e-4))
    # D/H profile
    DHprofs[1, :] = ncur[:D] ./ ncur[:H]  # altitude profile

    # Mixing ratios:
    # H2 lower atmospheric mixing ratio
    append!(atmo_metric_dict["H2MR"], sum(ncur[:H2][1:length(LA)])/sum([n_tot(ncur, h, n_alt_index) for h in LA]))
    # O2 Mixing ratio at surface
    append!(atmo_metric_dict["O2"], ncur[:O2][1]/N0)
    # CO mixing ratio
    append!(atmo_metric_dict["CO"], ncur[:CO][1]/N0)
    # CO/O2 ratio
    append!(atmo_metric_dict["CO/O2"], ncur[:CO][1]/ncur[:O2][1])
    # O3 mixing ratio
    append!(atmo_metric_dict["O3"], sum(ncur[:O3])*2e5) # gets O3 in #/cm^2

    # fractionation factor
    append!(atmo_metric_dict["f"], 2*(Df/Hf) / (ncur[:HDO][1]/ncur[:H2O][1]))

    return atmo_metric_dict
end

function analyze_water(abs_or_mr, allDbearers, make_plots=false, path=detailed_results_dir)
    #=
    Analyzes all the simulation results for water vapor experiments. Creates a dictionary object that
    contains various species concentrations, fluxes, and the fractionation factor.

    abs_or_mr: whether to store absolute concentrations or mixing ratios of species.
    allDbearers: boolean; whether to include the species OD, HDO2, and DO2
    make_plots: boolean; whether to make individual plots of the raw value of a parameter for a single simulation
    path: the directory storing the results for the detailed cases.

    returns: dictionary with the following keys and values:
        "H": absolute abundance of atomic H at the exobase
        "D": same, for atomic D
        "H2": same, for molecular hydrogen
        "HD": same, for hydrogen deuteride
        "H2MR": mixing ratio of H2 in the lower atmosphere (up to 80km) in ppm, used for comparing with krasnopolsky data
        "O2": absolute abundance of O2 at the surface
        "CO": same, but for CO
        "O3": ozone in #/cm^2 in the atmosphere (column abundance)
        "f": fractionation factor of the atmosphere
        "Hflux": escape flux of H due to loss of H, HD, and H2
        "Dflux": escape flux of D due to loss of D, HD
    =#
    # Establish parameters, filenames, etc
    wpaths = [path*"water_"*w for w in watervals_str]
    wfilelist = [path*"water_"*w*"/converged_water_"*w*".h5" for w in watervals_str]
    temps = [meanTs, meanTt, meanTe]
    oflux = 1.2e8
    q = abs_or_mr == "abs" ? " abundance" : " mixing ratio" # for labels
    mean_idx = findfirst(isequal(10), watervals) - 1
    subfolder = abs_or_mr == "abs" ? "abs/" : "mr/"

    # Establish variables to store data on simulations
    # Easier to deal with D/H profiles separately due to different array size
    DHprofs = Array{Any}(undef, length(watervals), length(alt)-2)
    wdict = Dict{String, Array}("O2"=>[], "HD"=>[], "H2"=>[], "H2MR"=>[], "H"=>[], "D"=>[], "CO"=>[],
                 "Hflux"=>[], "Dflux"=>[], "f"=>[], "CO/O2"=>[], "O3"=>[],
                 "DH"=>[], "H>Hloss"=>[], "H2>Hloss"=>[], "HD>Hloss"=>[], "D>Dloss"=>[], "HD>Dloss"=>[])
    if allDbearers
        wdict["OD"]=>[]
        wdict["HDO2"]=>[]
        wdict["DO2"]=>[]
    end

    # loop water files and collect data
    i = 1
    for (wpath, wfile) in zip(wpaths, wfilelist)
        # get the current array
        ncur = get_ncurrent(wfile)

        # if filling the absolute quantity arrays, subtract the flux loss.
        # subtract_flux is a binary value that is multiplied by the contribution value
        # to avoid writing excess if/then statements.
        subtract_flux = abs_or_mr == "abs" ? 1 : 0
        N0 = abs_or_mr == "abs" ? 1 : n_tot(ncur, 0, n_alt_index)
        Ntop = abs_or_mr == "abs" ? 1 : n_tot(ncur, 250e5, n_alt_index)
        LA = collect(0e5:2e5:78e5)

        # Calculate the things we care about
        # H and D fluxes
        Hf, contrib_H = get_flux(:H, wfile, temps, therm_only=true)
        Df, contrib_D = get_flux(:D, wfile, temps, therm_only=true)
        append!(wdict["Hflux"], Hf)
        append!(wdict["Dflux"], Df)

        append!(wdict["H>Hloss"], contrib_H[:H])
        append!(wdict["H2>Hloss"], contrib_H[:H2])
        append!(wdict["HD>Hloss"], contrib_H[:HD])
        append!(wdict["D>Dloss"], contrib_D[:D])
        append!(wdict["HD>Dloss"], contrib_D[:HD])

        # Concentrations:
        # H at exobase
        append!(wdict["H"], ncur[:H][end]/Ntop)
        # D at exobase
        append!(wdict["D"], ncur[:D][end]/Ntop)
        # H2 at exobase
        append!(wdict["H2"], ncur[:H2][end]/Ntop)
        # HD at exobase
        append!(wdict["HD"], ncur[:HD][end]/Ntop)

        # D/H stuff
        # D/H at 150 km
        append!(wdict["DH"], (ncur[:D][75]/ncur[:H][75])/(1.6e-4))
        # D/H profile
        DHprofs[i, :] = ncur[:D] ./ ncur[:H]  # altitude profile

        # Other stuff
        # H2 mixing ratio in lower atmo
        append!(wdict["H2MR"], sum(ncur[:H2][1:length(LA)])/sum([n_tot(ncur, h, n_alt_index) for h in LA])) # ppm in lower atmo
        # CO mixing ratio at surface
        append!(wdict["CO"], ncur[:CO][1]/N0)
        # O2 Mixing ratio at surface
        append!(wdict["O2"], ncur[:O2][1]/N0)
        # CO/O2 ratio at surface
        append!(wdict["CO/O2"], ncur[:CO][1]/ncur[:O2][1])
        # O3 in #/cm^2, used to convert to μm-atm later
        append!(wdict["O3"], sum(ncur[:O3])*2e5) # gets O3 in #/cm^2

        # fractionation factor
        append!(wdict["f"], 2*(Df/Hf) / (ncur[:HDO][1]/ncur[:H2O][1]))

        # Other D bearing species at exobase
        if allDbearers
            append!(wdict["OD"], ncur[:OD][end]/Ntop)
            append!(wdict["HDO2"], ncur[:HDO2][end]/Ntop)
            append!(wdict["DO2"], ncur[:DO2][end]/Ntop)
        end
        i += 1
    end

    if make_plots == true
        make_small_water_plots(watervals, wdict, DHprofs, q, mean_idx, subfolder)
    end

    return wdict
end

function analyze_Oflux(abs_or_mr, allDbearers, make_plots=false, path=detailed_results_dir)
    #=
    Analyzes all the simulation results for O flux experiments. Creates a dictionary object that
    contains various species concentrations, fluxes, and the fractionation factor.

    abs_or_mr: whether to store absolute concentrations or mixing ratios of species.
    allDbearers: boolean; whether to include the species OD, HDO2, and DO2
    make_plots: boolean; whether to make individual plots of the raw value of a parameter for a single simulation
    path: the directory storing the results for the detailed cases.

    returns: dictionary with the following keys and values:
        "H": absolute abundance of atomic H at the exobase
        "D": same, for atomic D
        "H2": same, for molecular hydrogen
        "HD": same, for hydrogen deuteride
        "H2MR": mixing ratio of H2 in the lower atmosphere (up to 80km) in ppm, used for comparing with krasnopolsky data
        "O2": absolute abundance of O2 at the surface
        "CO": same, but for CO
        "O3": ozone in #/cm^2 in the atmosphere (column abundance)
        "f": fractionation factor of the atmosphere
        "Hflux": escape flux of H due to loss of H, HD, and H2
        "Dflux": escape flux of D due to loss of D, HD
    =#
    # Establish important parameters, files, etc
    Ofluxvals = [8e7, 9e7, 1e8, 1.1e8, 1.2e8, 1.3e8, 1.4e8, 1.5e8, 1.6e8]
    Ofluxvals_str = ["8e7", "9e7", "1.0e8", "1.1e8", "1.2e8", "1.3e8", "1.4e8", "1.5e8", "1.6e8"]
    Ofilelist = [path*"Oflux_"*o*"/converged_Oflux_"*o*".h5"
                 for o in Ofluxvals_str]
    temps = [meanTs, meanTt, meanTe]
    mean_idx = findfirst(isequal(1.2e8), Ofluxvals) - 1
    q = abs_or_mr == "abs" ? " abundance " : " mixing ratio " # set label
    subfolder = abs_or_mr == "abs" ? "abs/" : "mr/"

    # Establish variables to store data on simulations
    odict = Dict("O2"=>[], "HD"=>[], "H2"=>[], "H2MR"=>[], "H"=>[], "D"=>[], "CO"=>[],
                 "Hflux"=>[], "Dflux"=>[], "f"=>[], "CO/O2"=>[], "O3"=>[],
                 "DH"=>[])
    if allDbearers
        odict["OD"]=>[]
        odict["HDO2"]=>[]
        odict["DO2"]=>[]
    end

    # Easier to deal with D/H profiles separately due to different array size
    DHprofs = Array{Any}(undef, length(Ofluxvals_str), length(alt)-2)

    # loop through O flux files
    i = 1
    for (oflux, ofile) in zip(Ofluxvals, Ofilelist)
        # get the current array
        ncur = get_ncurrent(ofile)

        N0 = abs_or_mr == "abs" ? 1 : n_tot(ncur, 0, n_alt_index)
        Ntop = abs_or_mr == "abs" ? 1 : n_tot(ncur, zmax, n_alt_index)
        LA = collect(0e5:2e5:78e5)

        # Calculate the things we care about
        # H and D fluxes
        Hf, contrib_H = get_flux(:H, ofile, temps, therm_only=true)
        Df, contrib_D = get_flux(:D, ofile, temps, therm_only=true)
        append!(odict["Hflux"], Hf)
        append!(odict["Dflux"], Df)

        # Concentrations:
        # H at exobase
        append!(odict["H"], ncur[:H][end]/Ntop)
        # D at exobase
        append!(odict["D"], ncur[:D][end]/Ntop)
        # H2 at exobase
        append!(odict["H2"], ncur[:H2][end]/Ntop)
        # HD at exobase
        append!(odict["HD"], ncur[:HD][end]/Ntop)

        # Calculate the things we care about

        # D/H stuff
        # D/H profile
        DHprofs[i, :] = ncur[:D] ./ ncur[:H]  # altitude profile
        # D/H at 150 km
        append!(odict["DH"], (ncur[:D][75]/ncur[:H][75])/(1.6e-4))

        # Other stuff
        # H2 mixing ratio in lower atmo
        append!(odict["H2MR"], sum(ncur[:H2][1:length(LA)])/sum([n_tot(ncur, h, n_alt_index) for h in LA])) # ppm in lower atmo
        # O2 Mixing ratio at surface
        append!(odict["O2"], ncur[:O2][1]/N0)
        # CO mixing ratio at surface
        append!(odict["CO"], ncur[:CO][1]/N0)
        # CO/O2 ratio at surface
        append!(odict["CO/O2"], ncur[:CO][1]/ncur[:O2][1])
        # O3 in #/cm^2, used to convert to μm-atm later
        append!(odict["O3"], sum(ncur[:O3])*2e5) # gets O3 in #/cm^2

        # fractionation factor
        append!(odict["f"], 2*(Df/Hf) / (ncur[:HDO][1]/ncur[:H2O][1]))

        # Other D bearing species abundances at exobase
        if allDbearers
            append!(odict["OD"], ncur[:OD][end]/Ntop)
            append!(odict["HDO2"], ncur[:HDO2][end]/Ntop)
            append!(odict["DO2"], ncur[:DO2][end]/Ntop)
        end
        i += 1
    end

    if make_plots == true
        make_small_Oflux_plots(Ofluxvals, Ofluxvals_str, odict, DHprofs, q, mean_idx,
                         subfolder)
    end

    return odict
end

function analyze_T(abs_or_mr, allDbearers, make_plots=false, path=detailed_results_dir)
    #=
    Analyzes all the simulation results for temperature profile experiments. Creates a dictionary object that
    contains various species concentrations, fluxes, and the fractionation factor.

    abs_or_mr: boolean; whether to store absolute concentrations or mixing ratios of species.
    allDbearers: boolean; whether to include the species OD, HDO2, and DO2
    make_plots: boolean; whether to make individual plots of the raw value of a parameter for a single simulation
    path: string; the directory storing the results for the detailed cases.

    returns: dictionary with the following keys and values:
        "H": absolute abundance of atomic H at the exobase
        "D": same, for atomic D
        "H2": same, for molecular hydrogen
        "HD": same, for hydrogen deuteride
        "H2MR": mixing ratio of H2 in the lower atmosphere (up to 80km) in ppm, used for comparing with krasnopolsky data
        "O2": absolute abundance of O2 at the surface
        "CO": same, but for CO
        "O3": ozone in #/cm^2 in the atmosphere (column abundance)
        "f": fractionation factor of the atmosphere
        "Hflux": escape flux of H due to loss of H, HD, and H2
        "Dflux": escape flux of D due to loss of D, HD
    =#

    # Set up parameters, filenames, etc
    # Dictionaries to store each experiment's data so it can be returned
    all_tdicts = Dict()

    tvals_str = Dict()
    # temppaths = Dict()
    tempfilelist = Dict()
    for k in keys(tvals)
        tvals_str[k] = [string(trunc(Int, x)) for x in tvals[k]]
        if k == "surface"
            # temppaths[k] = [path*"temp_"*t*"_$(meanTtint)_$(meanTeint)" for t in tvals_str[k]]
            tempfilelist[k] = [path*"temp_"*t*"_$(meanTtint)_$(meanTeint)"*"/converged_temp_"*t*"_$(meanTtint)_$(meanTeint).h5" for t in tvals_str[k]]  # TODO: revert if necessary
        elseif k == "tropopause"
            # temppaths[k] = [path*"temp_$(meanTsint)_"*t*"_$(meanTeint)" for t in tvals_str[k]]
            tempfilelist[k] = [path*"temp_$(meanTsint)_"*t*"_$(meanTeint)"*"/converged_temp_$(meanTsint)_"*t*"_$(meanTeint).h5" for t in tvals_str[k]]
        elseif k == "exobase"
            # temppaths[k] = [path*"temp_$(meanTsint)_$(meanTtint)_"*t for t in tvals_str[k]]
            tempfilelist[k] = [path*"temp_$(meanTsint)_$(meanTtint)_"*t*"/converged_temp_$(meanTsint)_$(meanTtint)_"*t*".h5" for t in tvals_str[k]]
        end
    end
    meanT = Dict("surface"=>meanTs, "tropopause"=>meanTt, "exobase"=>meanTe)  # nominal
    oflux = 1.2e8
    q = abs_or_mr == "abs" ? " abundance " : " mixing ratio " # set label
    subfolder = abs_or_mr == "abs" ? "abs/" : "mr/"

    # loop through which temp is varied and construct a list of datapoints
    for experiment in keys(tvals) # loop across the dictionary
        tdict = Dict("O2"=>[], "HD"=>[], "H2"=>[], "H2MR"=>[], "H"=>[], "D"=>[], "CO"=>[],
                     "Hflux"=>[], "Dflux"=>[], "f"=>[], "CO/O2"=>[], "O3"=>[],
                     "DH"=>[], "H>Hloss"=>[], "H2>Hloss"=>[], "HD>Hloss"=>[], "D>Dloss"=>[], "HD>Dloss"=>[])
        if allDbearers
            tdict["OD"]=>[]
            tdict["HDO2"]=>[]
            tdict["DO2"]=>[]
        end

        DHprofs = Array{Any}(undef, length(tvals_str[experiment]), length(alt)-2)


        # now loop through the values for each varied temp
        i = 1
        for (tv, tfile) in zip(tvals[experiment], tempfilelist[experiment])
            # set the temperature profile
            if experiment == "surface"
                temps = [tv, meanTt, meanTe]
            elseif experiment == "tropopause"
                temps = [meanTs, tv, meanTe]
            elseif experiment == "exobase"
                temps = [meanTs, meanTt, tv]
            end

            # get the current array
            ncur = get_ncurrent(tfile)
            N0 = abs_or_mr == "abs" ? 1 : n_tot(ncur, 0, n_alt_index)
            Ntop = abs_or_mr == "abs" ? 1 : n_tot(ncur, zmax, n_alt_index)
            LA = collect(0e5:2e5:78e5)

            # Calculate the things we care about
            # H and D fluxes
            Hf, contrib_H = get_flux(:H, tfile, temps, therm_only=true)
            Df, contrib_D = get_flux(:D, tfile, temps, therm_only=true)
            append!(tdict["Hflux"], Hf)
            append!(tdict["Dflux"], Df)

            # Concentrations:
            # H at exobase
            append!(tdict["H"], ncur[:H][end]/Ntop)
            # D at exobase
            append!(tdict["D"], ncur[:D][end]/Ntop)
            # H2 at exobase
            append!(tdict["H2"], ncur[:H2][end]/Ntop)
            # HD at exobase
            append!(tdict["HD"], ncur[:HD][end]/Ntop)

            # D/H stuff
            # D/H at 150 km
            append!(tdict["DH"], (ncur[:D][75]/ncur[:H][75])/(1.6e-4))
            # D/H profile
            DHprofs[i, :] = ncur[:D] ./ ncur[:H]  # altitude profile

            # Other stuff
            # H2 mixing ratio in lower atmo
            append!(tdict["H2MR"], sum(ncur[:H2][1:length(LA)])/sum([n_tot(ncur, h, n_alt_index) for h in LA])) # ppm in lower atmo
            # O2 Mixing ratio at surface
            append!(tdict["O2"], ncur[:O2][1]/N0)
            # CO mixing ratio at surface
            append!(tdict["CO"], ncur[:CO][1]/N0)
            # CO/O2 ratio at surface
            append!(tdict["CO/O2"], ncur[:CO][1]/ncur[:O2][1])
            # O3 in #/cm^2, used to convert to μm-atm later
            append!(tdict["O3"], sum(ncur[:O3])*2e5) # gets O3 in #/cm^2

            # fractionation factor
            append!(tdict["f"], 2*(Df/Hf) / (ncur[:HDO][1]/ncur[:H2O][1]))

            # Loss sources
            append!(tdict["H>Hloss"], contrib_H[:H])
            append!(tdict["H2>Hloss"], contrib_H[:H2])
            append!(tdict["HD>Hloss"], contrib_H[:HD])
            append!(tdict["D>Dloss"], contrib_D[:D])
            append!(tdict["HD>Dloss"], contrib_D[:HD])

            # Other D bearing species at exobase
            if allDbearers
                append!(tdict["OD"], ncur[:OD][end]/Ntop)
                append!(tdict["HDO2"], ncur[:HDO2][end]/Ntop)
                append!(tdict["DO2"], ncur[:DO2][end]/Ntop)
            end
            i += 1
        end

        if make_plots == true
            make_small_T_plots(tvals, tvals_str, tdict, DHprofs, experiment, q, meanT,
                         subfolder)
        end

        all_tdicts[experiment] = tdict
    end
    return all_tdicts
end

# Functions to make small supplementary plots for each detailed simulation --------

# Following two functions are subroutines common to make_small_water_plots, 
# make_small_Oflux_plots, make_small_T_plots
function DH_alt_prof_plot(DHproflist, exps, v, s, optext="", optlegend="")
    #=
    DHproflist: An array with D/H profiles in each row
    exps: a list of experiment identifiers as strings
    v: water, Oflux, or temp (for putting in correct folder)
    s: a subfolder "abs/" or "mr/" for either absolute abundances or mixing ratio
    optext: an optional extension that will be appended to v, for
            the temperature case where it goes in the temp_plots
            folder but the files also specify exobase, tropopause, etc.
    optlegend: an optional string to put into the legend
    =#
    # do the DH altitudinal profile plot
    # set up plot
    fig, ax = subplots(figsize=(6,4))
    plot_bg(ax)
    subplots_adjust(wspace=0, bottom=0.15)
    ax.set_xlabel("D/H ratio (in atomic D, H)")
    ax.set_ylabel("Altitude (km)")
    ax.set_yticks(ticks=collect(0:50:200))

    # generate colors
    c = get_grad_colors(length(exps), "viridis")

    # do actual plotting
    if typeof(exps[1]) != String
        exps = [string(x) for x in exps]
    end
    for i in range(1, length=length(exps))
        ax.plot(DHproflist[i, :], alt[2:end-1]./1e5, zorder=10, color=c[i, :],
                linewidth=3, label=optlegend*"="*exps[i])
    end

    # set savepath
    plotpath = detailed_results_dir*v*"_plots/"*s
    savepath = plotpath*v*optext*"_DH_prof.png"
    legend(fontsize=12, bbox_to_anchor=[1.01,1], loc=2, borderaxespad=0)
    savefig(savepath, bbox_inches="tight")

    # save it again but with log x-axis
    xscale("log")
    # xlim(3.5e-4,5e-4)
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
    paststudies = Dict("water"=>Dict("yung"=>15, "nair"=>[3,8.8]), # these are in pr μm
                            "Oflux"=>Dict("yung"=>5, "nair"=>9), # indices
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
    medgray = "#444444"

    # CO and O2 mixin ratio axis
    ax.plot(xvals, ydict["CO"], marker="o", zorder=10, color="navy", label="CO", alpha=0.5)
    ax.plot(xvals, ydict["O2"], marker="D", zorder=10, color="navy", label=L"O$_2$", alpha=0.5)

    ax.axvline(meanX, color=medgray, zorder=10)
    nomtexty = Dict("temp_exobase"=>8.5e-4, "temp_surface"=>2e-3,
                    "temp_tropopause"=>1.2e-3, "water"=>4e-3, "Oflux"=>6e-4)
    ax.text(meanX, nomtexty[lookupkey], "Nomimal\ncase", color=medgray)
    ax.set_facecolor("#ededed")
    ax.tick_params(axis="y", labelcolor="navy", which="both")
    ax.set_ylabel(ystr, color="navy")
    ax.set_yscale("log")
    # ax.legend()
    for side in ["top", "bottom", "left", "right"]
        ax.spines[side].set_visible(false)
    end

    # CO/O2 axis
    ax2 = ax.twinx() # set up the other axis
    ax2.axhline(0.6, color="red", linestyle="--")
    obstxtx = Dict("temp_exobase"=>170, "temp_tropopause"=>80, "temp_surface"=>250,
                   "water"=>1e-1, "Oflux"=>0)
    obstxty = Dict("temp_exobase"=>0.55, "temp_tropopause"=>0.7, "temp_surface"=>0.55,
                   "water"=>0.55, "Oflux"=>0.7)
    ax2.text(obstxtx[lookupkey], obstxty[lookupkey], "Obs.", color="red", va="top")
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
    pastred = "red"
    if pathkey=="water"
        ax2.errorbar(paststudies[lookupkey]["nair"][1], 0.14, capsize=2, zorder=10,
                     yerr=reshape([0.05; 0.43], 2, 1), color=pastred,
                     ecolor=pastred, marker="v")  # Nair 1994 low water
        ax2.errorbar(paststudies[lookupkey]["nair"][2], 0.1, capsize=2, zorder=10,
                     yerr=reshape([0.07; 0.4], 2, 1), color=pastred,
                     ecolor=pastred, marker="^")  # Nair 1994 high water
    else
        ax2.errorbar(paststudies[lookupkey]["nair"], 0.1, capsize=2, zorder=10,
                     yerr=reshape([0.07; 0.4], 2, 1), color=pastred,
                     ecolor=pastred, marker="^")  # Nair 1994 high water
    end
    ax2.scatter(paststudies[lookupkey]["yung"], 0.53, marker="*", s=100,
                color=pastred, zorder=10)  # Yung 1988

    # more ticks for exobase surface temp
    if lookupkey == "temp_exobase"
        ax.set_xticks(ticks=vcat(xvals, [325, 350]))
        [l.set_visible(false) for (i,l) in enumerate(ax.xaxis.get_ticklabels()) if i % 2 != 0]
    end

    plotpath = detailed_results_dir*pathkey*"_plots/"*s
    savepath = plotpath*lookupkey*"_CO_and_O2.png"
    savefig(savepath, bbox_inches="tight")
    close(fig)
end

function make_small_water_plots(water_x, d, DHdata, q, nom_i, s)
    #=
    Makes the plots for water vapor experiments

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
    
    rcParams["font.family"] = "sans-serif"
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 22
    rcParams["axes.labelsize"]= 24
    rcParams["xtick.labelsize"] = 22
    rcParams["ytick.labelsize"] = 22

    xlab = "Water"*q*L" (pr $\mu$m)"

    # Do the individual plots
    plots = filter!(e->e∉["CO", "O2", "CO/O2", "HD", "OD", "DO2", "HDO2",
                          "H>Hloss", "H2>Hloss", "HD>Hloss", "D>Dloss",
                          "HD>Dloss"], [k for k in keys(d)])
    for i in plots
        # set up plot
        fig, ax = subplots(figsize=(6,4))
        plot_bg(ax)
        xscale("log")
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
                        "H2MR"=>"H"*L"_2"*q,
                        "H2"=>"H"*L"_2"*q,
                        "O3"=>L"O$_3$ (#/cm$^{-2}$)")
        ylabel(ylabdict[i])

        # Certain measureables have a range not suited to log scale
        nologplease = ["Dflux", "CO/O2", "DH"]  # don't logscale everything
        if ~(i in nologplease)
            yscale("log")
            if i in ["Hflux", "HD", "H", "H2", "O3", "O2"]
                ylim(minimum(d[i])/2, maximum(d[i])*2)
            end
        end

        # set savepath
        plotpath = detailed_results_dir*"water_plots/" * s
        savepath = "water_"*i*".png"

        savefig(plotpath*savepath, bbox_inches="tight")
        close(fig)
    end

    # make special plots
    CO_O2_plot(water_x, d, xlab, "water", nom_i, s)
    DH_alt_prof_plot(DHdata, water_x, "water", s)
end

function make_small_Oflux_plots(phiO, phiO_str, d, DHdata, q, nom_i, s)
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
    rcParams["font.family"] = "sans-serif"
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 22
    rcParams["axes.labelsize"]= 24
    rcParams["xtick.labelsize"] = 22
    rcParams["ytick.labelsize"] = 22

    xlab = L"$\phi_O$ (cm$^{-2}$s$^{-1}$)"

    # Make individual plots showing change
    plots = filter!(e->e∉["CO", "O2", "CO/O2", "HD", "OD", "DO2", "HDO2"], [k for k in keys(d)])
    for i in plots
        # basic plot stuff
        fig, ax = subplots(figsize=(6,4))
        plot_bg(ax)
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
                        "H2MR"=>"H"*L"_2"*q,
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
        plotpath = detailed_results_dir*"Oflux_plots/"*s
        savepath = plotpath*"O_flux_"*i*".png"

        savefig(savepath, bbox_inches="tight")
        close(fig)
    end

    # make special plots
    CO_O2_plot(phiO_str, d, xlab, "Oflux", nom_i, s)
    DH_alt_prof_plot(DHdata, phiO_str, "Oflux", s)
end

function make_small_T_plots(T, T_str, d, DHdata, exp, q, nomT, s)
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
    rcParams["font.family"] = "sans-serif"
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 22
    rcParams["axes.labelsize"]= 24
    rcParams["xtick.labelsize"] = 22
    rcParams["ytick.labelsize"] = 22

    xlab = exp*" temperature (K)"

    # loop through the parameters of interest and plot them
    plots = filter!(e->e∉["CO", "O2", "CO/O2", "HD", "OD", "DO2", "HDO2",
                          "H>Hloss", "H2>Hloss", "HD>Hloss", "D>Dloss",
                          "HD>Dloss"], [k for k in keys(d)])
    for i in plots
        # basic plot stuff
        fig, ax = subplots(figsize=(6,4))
        plot_bg(ax)
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
                        "H2MR"=>"H"*L"_2"*q,
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
        plotpath = detailed_results_dir*"temp_plots/"*s
        savepath = plotpath*join(["temp", exp], "_")*"_"*i*".png"
        savefig(savepath, bbox_inches="tight")
        close(fig)
    end

    # make special plots
    CO_O2_plot(T[exp], d, exp*" temperature", "temp", nomT[exp], s, "_"*exp)
    DH_alt_prof_plot(DHdata, T_str[exp], "temp", s, "_"*exp,
                     latexstring("T_{$(exp[1:3])}"))
end

# Functions to make the key results plots for detailed cases -----------------------
# TODO: this needs to become 3 functions like the water and temp sections below
function make_Oflux_main_plots(output_MR, output_abs)
    #=
    output_MR:
    output_abs:
    ex: experiment type, the usual key of "water", "O flux", "surface" etc.
    SVP: whether to plot for experiments where SVP was held constant
         (SVP="const") to one temp profile, or for experiments where SVP was
         allowed to follow the temp profile (SVP="vary")
    =#
    # make plots pretty
    rcParams["font.family"] = "sans-serif"
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 22
    rcParams["axes.labelsize"]= 24
    rcParams["xtick.labelsize"] = 22
    rcParams["ytick.labelsize"] = 22

    # first plot - compare with data -------------------------------------------
    fig, ax1 = subplots(figsize=(7,6))
    plot_bg(ax1)

    c1 = ["#10007A", "#2F7965", "#e16262", "#D59A07"]  # just colors

    # calculate the relativeness of each point wrt data
    COdiff = (output_MR["CO"] .- data[1])/s[1]
    O2diff = (output_MR["O2"] .- data[2])/s[2]
    # H2diff = (ABSdictvar["H2"] .- data[3])/s[3]  # absolute abundance
    H2diff = (output_MR["H2MR"]./1e-6 .- data[3])/s[3] # ppm
    O3diff = (areadensity_to_micron_atm(output_abs["O3"]) .- data[4])/s[4]


    # plot the actual stuff
    ax1.plot(Oflux_vals, COdiff, marker="o", ms=sz, color=c1[1], zorder=10)
    ax1.plot(Oflux_vals, O2diff, marker="x", ms=sz, color=c1[2], zorder=10)
    ax1.plot(Oflux_vals, H2diff, marker="*", ms=sz, color=c1[3], zorder=10)
    ax1.plot(Oflux_vals, O3diff, marker="^", ms=sz, color=c1[4], zorder=10)
    # plot the mean values
    ax1.axvline(nom_i_py["Oflux"], color=medgray, zorder=5)
    ax1.axhline(0, color="black")

    # Text on plots
    ax1.text(8.7, 1, "Global\nmean", color=medgray, ha="right")

    ax1.set_xlabel(L"\Phi_O"*L" (cm$^{-2}s^{-1}$)")
    ax1.set_ylabel(L"($X_{model}$-$X_{obs}$)/$\sigma$")
    savefig(detailed_results_dir*"output_vs_data_oflux.png", bbox_inches="tight")
    savefig(results_dir*"AllResultPlots/output_vs_data_oflux.png", bbox_inches="tight")


    # second plot: H2, HD, H, D, Hflux, Dflux ----------------------------------
    fig, ax2 = subplots(figsize=(7,6))

    c2 = ["#6270d9","#95b034","#d14f58","#5ca85c","#ce6d30", "#0d186d"]

    # to do the division/normalization we need to re-find the index of the value
    # against which to normalize. To do this, look in the xvals by experiment,
    # then index that according to whether we cut it or not. There doesn't seem
    # to be a cleaner way to do the indexing of cut.
    normidx = findfirst(isequal("1.2e8"), Oflux_vals[1:1:length(Oflux_vals)])

    HDdiff = normalize(output_MR["HD"], normidx)
    H2diff = normalize(output_MR["H2"], normidx)
    Hdiff = normalize(output_MR["H"], normidx)
    Ddiff = normalize(output_MR["D"], normidx)
    Hfdiff = normalize(output_MR["Hflux"], normidx)
    Dfdiff = normalize(output_MR["Dflux"], normidx)

    # plot the actual stuff
    ax2.plot(Oflux_vals, Hdiff, marker="x", ms=sz,  color=c2[2], zorder=10, label="H")
    ax2.plot(Oflux_vals, Ddiff, marker="*", ms=sz,  color=c2[5], zorder=10, label="D")
    ax2.plot(Oflux_vals, H2diff, marker="o", ms=sz, color=c2[1], zorder=10, label="H2")
    ax2.plot(Oflux_vals, HDdiff, marker="o", ms=sz, color=c2[6], zorder=10, label="HD")
    ax2.plot(Oflux_vals, Hfdiff, marker="^", ms=sz, color=c2[4], zorder=10, label=L"\phi_H")
    ax2.plot(Oflux_vals, Dfdiff, marker="D", ms=sz, color=c2[3], zorder=10, label=L"\phi_D")
    # nominal value
    ax2.axvline(nom_i_py["Oflux"], color=medgray, zorder=5)
    ax2.axhline(1, color=medgray, zorder=5)

    # Text on plots
    ax2.text(8.7, 1, "Global\nmean", color=medgray, ha="right")

    # Other species that can be plotted, or species at other locations
    # HDtopdiff = normalize(MRdictvar["HDtop"], normidx[ex])
    # ODtopdiff = normalize(MRdictvar["HDtop"], normidx[ex])
    # HDO2topdiff = normalize(MRdictvar["HDtop"], normidx[ex])
    # DO2topdiff = normalize(MRdictvar["HDtop"], normidx[ex])
    # ax2.plot(xvals[ex], HDtopdiff, marker="o", color="#7D3403", zorder=10)
    # ax2.plot(xvals[ex], ODtopdiff, marker="o", color="red", zorder=10)
    # ax2.plot(xvals[ex], HDO2topdiff, marker="o", color="blue", zorder=10)
    # ax2.plot(xvals[ex], DO2topdiff, marker="o", color="purple", zorder=10)

    ax2.set_xlabel(L"\Phi_O"*L" (cm$^{-2}s^{-1}$)")
    ax2.set_ylabel(L"X/X$_{nominal}$")
    #ax2.legend(bbox_to_anchor=(1.05, 1))
    savefig(detailed_results_dir*"compare_nominal_Oflux.png", bbox_inches="tight")
    savefig(results_dir*"AllResultPlots/compare_nominal_Oflux.png", bbox_inches="tight")

    # # third plot: f ----------------------------------------------------------
    # only for water and O flux - the others are done in a panel
    fig, ax = subplots(figsize=(7,6))
    plot_bg(ax)

    ax.tick_params(rotation=45, axis="x")
    ax.text(8.7, 0.00112, "Global\nmean", color=medgray, ha="right")


    ax.set_ylabel(L"$f$", color="black")
    ax.set_xlabel(L"\Phi_O"*L" (cm$^{-2}s^{-1}$)")
    linecol = "cornflowerblue"

    ax.plot(Oflux_vals, output_MR["f"], marker="o", color=linecol, zorder=10)
    ax.axvline(nom_i_py["Oflux"], color=medgray, zorder=5)


    savefig(detailed_results_dir*"f_tradeoff_Oflux.png", bbox_inches="tight")
    savefig(results_dir*"AllResultPlots/f_tradeoff_Oflux.png", bbox_inches="tight")
end

# Key results plots for detailed cases, water vapor:
function make_water_Hspecies_plot(output_dict, abs_or_mr)
    #=
    Plot showing H, D, HD, H2, ΦH and ΦD as relative to the global mean profile
    =#

    # to do the division/normalization we need to re-find the index of the value
    # against which to normalize. To do this, look in the xvals by experiment,
    # then index that according to whether we cut it or not. There doesn't seem
    # to be a cleaner way to do the indexing of cut.
    normidx = findfirst(isequal(10), watervals[1:1:length(watervals)])

    # make plots pretty
    rcParams["font.family"] = "sans-serif"
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 18
    rcParams["axes.labelsize"]= 20
    rcParams["xtick.labelsize"] = 18
    rcParams["ytick.labelsize"] = 18

    # Make the plot
    fig, ax = subplots(figsize=(7,6))
    plot_bg(ax)

    HDdiff = normalize(output_dict["HD"], normidx)
    H2diff = normalize(output_dict["H2"], normidx)
    Hdiff = normalize(output_dict["H"], normidx)
    Ddiff = normalize(output_dict["D"], normidx)
    Hfdiff = normalize(output_dict["Hflux"], normidx)
    Dfdiff = normalize(output_dict["Dflux"], normidx)

    # plot the actual stuff

    # colors for each line
    c = [Hcolor, Dcolor, H2color, HDcolor, Hflux_color, Dflux_color]

    ax.plot(watervals, Hdiff, marker="x", ms=sz,  color=c[1], zorder=10, label="H")
    ax.plot(watervals, Ddiff, marker="*", ms=sz+10,  color=c[2], zorder=10,
            label="D")#, alpha=0.7)

    ax.plot(watervals, H2diff, marker="o", ms=sz, color=c[3], zorder=9, label="H2")
    ax.plot(watervals, HDdiff, marker="D", ms=sz, color=c[4], zorder=10,
            label="HD")#, alpha=0.7)

    ax.plot(watervals, Hfdiff, marker="^", ms=sz, color=c[5], zorder=9, label=L"\phi_H")
    ax.plot(watervals, Dfdiff, marker="v", ms=sz, color=c[6], zorder=10,
            label=L"\phi_D")#, alpha=0.7)
    # nominal value
    ax.axvline(nom_i_py["water"], color=medgray, zorder=5)
    ax.axhline(1, color=medgray, zorder=5)

    # Text on plots
    ax.set_xscale("log")
    ax.text(11, 0.9, "Standard\ntemperature", color=medgray, ha="left", va="top")
    ax.text(1, 0.99, "H", color=c[1], ha="left", va="top")
    ax.text(1.1, 0.85, "D", color=c[2], ha="left", va="top")
    ax.text(1, 1.02, L"H_2", color=c[3], ha="left", va="top")
    ax.text(1.3, 0.87, "HD", color=c[4], ha="left", va="top")
    ax.text(1.5, 0.995, L"$\phi_H$", color=c[5], ha="left", va="top")
    ax.text(1.7, 0.9, L"$\phi_D$", color=c[6], ha="left", va="top")

    # Other species that can be plotted, or species at other locations
    # HDtopdiff = normalize(MRdictvar["HDtop"], normidx[ex])
    # ODtopdiff = normalize(MRdictvar["HDtop"], normidx[ex])
    # HDO2topdiff = normalize(MRdictvar["HDtop"], normidx[ex])
    # DO2topdiff = normalize(MRdictvar["HDtop"], normidx[ex])
    # ax2.plot(xvals[ex], HDtopdiff, marker="o", color="#7D3403", zorder=10)
    # ax2.plot(xvals[ex], ODtopdiff, marker="o", color="red", zorder=10)
    # ax2.plot(xvals[ex], HDO2topdiff, marker="o", color="blue", zorder=10)
    # ax2.plot(xvals[ex], DO2topdiff, marker="o", color="purple", zorder=10)

    ax.set_xlabel(L"total atmospheric water (pr $\mu$m)")
    ax.set_ylabel(L"X/X(\overline{T}_{standard})")
    savefig(detailed_results_dir*"compare_nominal_water_"*abs_or_mr*".png", bbox_inches="tight", dpi=300)
    savefig(results_dir*"AllResultPlots/compare_nominal_water_"*abs_or_mr*".png", bbox_inches="tight", dpi=300)
end

function make_water_output_vs_data(output_MR, output_abs)
    #=TODO: Fill me in=#
    rcParams["font.family"] = "sans-serif"
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 18
    rcParams["axes.labelsize"]= 20
    rcParams["xtick.labelsize"] = 18
    rcParams["ytick.labelsize"] = 18


    # first plot - compare with data -------------------------------------------
    fig, ax = subplots(figsize=(7,6))
    plot_bg(ax)

    c1 = [COcolor, O2color, H2color, O3color] 

    # calculate the relativeness of each point wrt data
    COdiff = (output_MR["CO"] .- data[1])/s[1]
    O2diff = (output_MR["O2"] .- data[2])/s[2]
    # H2diff = (ABSdictvar["H2"] .- data[3])/s[3]  # absolute abundance
    H2diff = (output_MR["H2MR"]./1e-6 .- data[3])/s[3] # ppm
    O3diff = (areadensity_to_micron_atm(output_abs["O3"]) .- data[4])/s[4]

    # plot the actual stuff
    ax.plot(watervals, COdiff, marker="s", ms=sz, color=c1[1], zorder=10)
    ax.plot(watervals, O2diff, marker="+", ms=sz, color=c1[2], zorder=10)
    ax.plot(watervals, H2diff, marker="o", ms=sz, color=c1[3], zorder=10)
    ax.plot(watervals, O3diff, marker="p", ms=sz, color=c1[4], zorder=10)
    # plot the mean values
    ax.axvline(nom_i_py["water"], color=medgray, zorder=5)
    ax.axhline(0, color="black")

    # Text on plots
    ax.set_xscale("log")
    ax.text(10.2, 2.5, "Standard\ntemperature", color=medgray, ha="left", va="top")
    ax.text(1, -2.4, "CO", color=c1[1], ha="left", va="top")
    ax.text(20, -2.5, L"O$_2$", color=c1[2], ha="left", va="top")
    ax.text(1, -0.5, L"H$_2$", color=c1[3], ha="left", va="top")
    ax.text(1, 2.5, L"O$_3$", color=c1[4], ha="left", va="top")
    ax.set_xlabel(L"total atmospheric water (pr $\mu$m)")
    ax.set_ylabel(L"($X_{model}$-$X_{obs}$)/$\sigma$")
    savefig(detailed_results_dir*"output_vs_data_water.png", bbox_inches="tight", dpi=300)
    savefig(results_dir*"AllResultPlots/output_vs_data_water.png", bbox_inches="tight", dpi=300)
end

function make_water_f_plot(output_MR)
    #=
    f as a function of water vapor
    =#

    # make plots pretty
    rcParams["font.family"] = "sans-serif"
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 18
    rcParams["axes.labelsize"]= 20
    rcParams["xtick.labelsize"] = 18
    rcParams["ytick.labelsize"] = 18

    fig, ax = subplots(figsize=(7,6))
    plot_bg(ax)
    
    ax.plot(watervals, output_MR["f"], marker="o", ms=sz, color=f_main_color, zorder=10)
    # ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(20))
    # ax.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(10))
    ax.set_ylim(0.0015, 0.002)

    # ax.set_yticks([0.0015, 0.0016, 0.0017, 0.0018, 0.0019, 0.002])
    ax.axvline(nom_i_py["water"], color=medgray, zorder=5)
    ax.text(9, 0.00198, "Standard\ntemperature", color=medgray, ha="right", va="top")

    ax.set_xscale("log")
    ax.tick_params(axis="y", labelcolor=f_main_color, which="both")
    ax.set_ylabel(L"$f$", color="black")
    ax.set_xlabel(L"total atmospheric water (pr $\mu$m)")

    # plot the secondary axis showing increase
    nomdict = make_std_atmo_dict("mr")
    fdiff = normalize_val(output_MR["f"], nomdict["f"])
    ax_2 = ax.twinx()
    ax_2.set_ylabel(L"$f$/$f(\overline{T}_{standard})$", color=f_sec_color, fontsize=24)
    ax_2.plot(watervals, fdiff, marker="o", ms=sz, color=f_sec_color)
    ax_2.tick_params(axis="y", labelcolor=f_sec_color, which="both")
    for side in ["top", "bottom", "left", "right"]
        ax_2.spines[side].set_visible(false)
    end
    ax_2.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.05))
    ax_2.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(0.01))
    ax_2.set_ylim(0.8, 1.1)

    savefig(detailed_results_dir*"f_vs_water.png", bbox_inches="tight", dpi=300)
    savefig(results_dir*"AllResultPlots/f_vs_water.png", bbox_inches="tight", dpi=300)
end

function make_water_loss_contribs_plot(output_abs)
    #=
    Make a plot showing how much H, D, HD, and H2 contribute to loss of H, D in the water
    experiments.
    =#

    # make plots pretty
    rcParams["font.family"] = "sans-serif"
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 18
    rcParams["axes.labelsize"]= 20
    rcParams["xtick.labelsize"] = 18
    rcParams["ytick.labelsize"] = 18

    mpltic = pyimport("matplotlib.ticker")

    # FIGURE: Just exobase =======================================================================
    fig, ax = subplots(1, 2, sharex=true, sharey=true, figsize=(14, 5))

    for i in 1:2
        plot_bg(ax[i])
        ax[i].tick_params(which="both", right=true)
        ax[i].yaxis.set_major_locator(mpltic.MultipleLocator(0.2))
        ax[i].yaxis.set_minor_locator(mpltic.MultipleLocator(0.05))
        ax[i].set_yticks(0:0.2:1, minor=false)
        ax[i].set_ylim(0, 1)
    end

    # H loss -----------------------------------------------------
    total_H_flux = output_abs["Hflux"]
    HDfrac = output_abs["HD>Hloss"]./total_H_flux
    H2frac = output_abs["H2>Hloss"]./total_H_flux
    Hfrac = output_abs["H>Hloss"]./total_H_flux

    ax[1].fill_between(watervals, 0, (1 .- Hfrac-H2frac-HDfrac), color=HDcolor, zorder=10, alpha=0.6)  # HD
    ax[1].fill_between(watervals, (1 .- Hfrac-H2frac), (1 .- Hfrac), color=H2color, zorder=10, alpha=0.6)  # H2
    ax[1].fill_between(watervals, (1 .- Hfrac), 1, color=Hcolor, zorder=10, alpha=0.6)  # H
    ax[1].text(0.1, 0.8, "atomic H", color=Hcolor_dark, transform=ax[1].transAxes, fontsize=24, zorder=20)

    # D Loss -----------------------------------------------------
    total_D_flux = output_abs["Dflux"]
    Dfrac = output_abs["D>Dloss"]./total_D_flux
    HDfrac = output_abs["HD>Dloss"]./total_D_flux
    println(output_abs["HD"])

    ax[2].fill_between(watervals, (1 .- Dfrac), 1, color=Dcolor, zorder=10, alpha=0.6)  # D
    ax[2].fill_between(watervals, 0, (1 .- Dfrac), color=HDcolor, zorder=10, alpha=0.6) # HD
    ax[2].text(0.1, 0.8, "atomic D", color=Dcolor, transform=ax[2].transAxes, fontsize=24, zorder=20)

    # axis format ------------------------------------------------
    mpltic = pyimport("matplotlib.ticker")
    ax[1].set_xlabel(L"total atmospheric water (pr $\mu$m)", fontsize=24)
    ax[1].set_ylabel("Fraction of total loss", fontsize=24)


    # plot panel labels ------------------------------------------
    ax[1].text(0.8, 0.1, L"H_2", color=H2color_dark, transform=ax[1].transAxes, va="top", fontsize=24, zorder=20)
    ax[2].text(0.8, 0.1, "HD", color=HDcolor, transform=ax[2].transAxes, va="top", fontsize=24, zorder=20)
    # ax[2].plot([335, 345], [0.06, 0.02], color=HDcolor)

    ax[1].text(0.01, 1, "a", color="black", fontsize=26, weight="bold", va="top", zorder=20, transform=ax[1].transAxes)
    ax[2].text(0.01, 1, "b", color="black", fontsize=26, weight="bold", va="top", zorder=20, transform=ax[2].transAxes)
    savefig(detailed_results_dir*"loss_contributions_water.png", bbox_inches="tight", dpi=300)
end

# Key results plots for detailed cases, temperatures:
function make_T_Hspecies_plot(output_dict, abs_or_mr)
    #=
    Plot showing H, D, HD, H2, ΦH and ΦD as relative to the global mean profile
    =#

    # where to place text, basically.
    exps = ["surface", "tropopause", "exobase"]
    gmean_txt_loc = Dict("surface"=>[219, 1.4], "tropopause"=>[131, 1.9],
                         "exobase"=>[210, 100])
    # ylims = Dict("surface"=>[-10, 190], "tropopause"=>[-4, 4.2], "exobase"=>[-4, 4.2])
    xt_args = Dict("surface"=>150:20:280, "tropopause"=>100:10:160,
                   "exobase"=>150:50:350)
    if abs_or_mr == "mr"
        linelbls = DataFrame(Exp=["surface", "tropopause", "exobase"],
                              H=[[150, 1.8],     [100, 2.4],  [150, 10]],
                              D=[[260, 0.5],     [100, 1.5],  [340, 0.1]],
                              H2=[[190, 2],      [105, 3.2], [150, 2.2]],
                              HD=[[150, 1.5],    [100, 2.9], [340, 0.8]],
                              fluxH=[[150, 1],  [100, 0.9], [340, 2.2]],
                              fluxD=[[150, 0.7], [100, 0.4],  [155, 0.4]])
    else
        linelbls = DataFrame(Exp=["surface", "tropopause", "exobase"],
                              H=[[150, 1.1],   [150, 1.12],  [150, 30]],
                              D=[[175, 0.75],   [108, 0.75],  [340, 0.11]],
                              H2=[[195, 1.4],   [110, 1.35], [340, 0.19]],
                              HD=[[165, 1.2],  [100, 1.15], [340, 0.3]],
                              fluxH=[[150, 0.97],  [155, 1.15], [150, 0.8]],
                              fluxD=[[160, 0.7], [100, 0.7],  [145, 0.1]])
    end

     # colors for each
     # order: H, D, H2, HD, ΦH, ΦD.
    c = [Hcolor, Dcolor, H2color, HDcolor, Hflux_color, Dflux_color] # H species #3, colorblind safe
    # to do the division/normalization we need to re-find the index of the value
    # against which to normalize. To do this, look in the xvals by experiment,
    # then index that according to whether we cut it or not. There doesn't seem
    # to be a cleaner way to do the indexing of cut.
    normidx = Dict("exobase"=>findfirst(isequal(meanTe), tvals["exobase"][1:1:length(tvals["exobase"])]),
                   "tropopause"=>findfirst(isequal(meanTt), tvals["tropopause"][1:1:length(tvals["tropopause"])]),
                   "surface"=>findfirst(isequal(meanTs), tvals["surface"][1:1:length(tvals["surface"])]),
                   )
    nomdict = make_std_atmo_dict(abs_or_mr) # get info for nom case which has values not % by 10

    # make plots pretty
    rcParams["font.family"] = "sans-serif"
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 18
    rcParams["axes.labelsize"]= 20
    rcParams["xtick.labelsize"] = 18
    rcParams["ytick.labelsize"] = 18

    # Make the plot
    fig, ax = subplots(1, 3, sharex=false, sharey=false, figsize=(21, 5))
    subplots_adjust(wspace=0.3)

    for i in 1:3
        plot_bg(ax[i])
        ex = exps[i]

        Hdiff = normalize_val(output_dict[ex]["H"], nomdict["H"])
        Ddiff = normalize_val(output_dict[ex]["D"], nomdict["D"])
        H2diff = normalize_val(output_dict[ex]["H2"], nomdict["H2"])
        HDdiff = normalize_val(output_dict[ex]["HD"], nomdict["HD"])
        Hfdiff = normalize_val(output_dict[ex]["Hflux"], nomdict["Hflux"])
        Dfdiff = normalize_val(output_dict[ex]["Dflux"], nomdict["Dflux"])

        # plot the actual stuff
        ax[i].plot(tvals[ex], Hdiff, marker="x", ms=sz, color=c[1], zorder=10, label="H")
        ax[i].plot(tvals[ex], Ddiff, marker="*", ms=sz, color=c[2], zorder=10, label="D")
        ax[i].plot(tvals[ex], H2diff, marker="o", ms=sz, color=c[3], zorder=10, label="H2")
        ax[i].plot(tvals[ex], HDdiff, marker="D", ms=sz, color=c[4], zorder=10, label="HD")
        ax[i].plot(tvals[ex], Hfdiff, marker="^", ms=sz, color=c[5], zorder=9, label=L"\phi_H")
        ax[i].plot(tvals[ex], Dfdiff, marker="v", ms=sz, color=c[6], zorder=9, label=L"\phi_D")
        # nominal value
        ax[i].axvline(nom_i_py[ex], color=medgray, zorder=5)
        ax[i].text(gmean_txt_loc[ex][1], gmean_txt_loc[ex][2], "Standard\ntemperature",
                   color=medgray, ha="left", va="top")
        ax[i].axhline(1, color=medgray, zorder=5)

        # text on plots
        dfentry = linelbls[linelbls.Exp.==ex, :]
        ax[i].text(dfentry.H[1][1], dfentry.H[1][2], "H", color=c[1], ha="left", va="top")
        ax[i].text(dfentry.D[1][1], dfentry.D[1][2], "D", color=c[2], ha="left", va="top")
        ax[i].text(dfentry.H2[1][1], dfentry.H2[1][2], "H2", color=c[3], ha="left", va="top")
        ax[i].text(dfentry.HD[1][1], dfentry.HD[1][2], "HD", color=c[4], ha="left", va="top")
        ax[i].text(dfentry.fluxH[1][1], dfentry.fluxH[1][2], L"\phi_H", color=c[5], ha="left", va="top")
        ax[i].text(dfentry.fluxD[1][1], dfentry.fluxD[1][2], L"\phi_D", color=c[6], ha="left", va="top")


        # Various configurations and such
        if ex=="surface"
            ax[i].xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(20))
            ax[i].xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(10))
        elseif ex=="tropopause"
            ax[i].xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(10))
            # ax[i].yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(1))
        elseif ex=="exobase"
            ax[i].set_yscale("log")
        end

        ax[i].set_xlabel(ex*" temperature (K)")
    end
    ax[1].set_ylabel(L"X/X(\overline{T}_{standard})")
    ax[1].text(150, 1.4, "a", color="black", weight="bold", va="top",
               ha="center", fontsize=26)#
    ax[2].text(100, 1.9, "b", color="black", weight="bold", va="top",
               ha="center", fontsize=26)
    ax[3].text(150, 100, "c", color="black", weight="bold", va="top",
               ha="center", fontsize=26)

    # Other species that can be plotted, or species at other locations
    # HDtopdiff = normalize(MRdictvar["HDtop"], normidx[ex])
    # ODtopdiff = normalize(MRdictvar["HDtop"], normidx[ex])
    # HDO2topdiff = normalize(MRdictvar["HDtop"], normidx[ex])
    # DO2topdiff = normalize(MRdictvar["HDtop"], normidx[ex])
    # ax.plot(xvals[ex], HDtopdiff, marker="o", color="#7D3403", zorder=10)
    # ax.plot(xvals[ex], ODtopdiff, marker="o", color="red", zorder=10)
    # ax.plot(xvals[ex], HDO2topdiff, marker="o", color="blue", zorder=10)
    # ax.plot(xvals[ex], DO2topdiff, marker="o", color="purple", zorder=10)

    savefig(detailed_results_dir*"compare_nominal_temps_"*abs_or_mr*".png", bbox_inches="tight", dpi=300)
    savefig(results_dir*"AllResultPlots/compare_nominal_temps_"*abs_or_mr*".png", bbox_inches="tight", dpi=300)
end

function make_T_output_vs_data(output_MR, output_abs)
    #=
    Makes the plots comparing simulation output with observational data.

    =#

    # set up
    exps = ["surface", "tropopause", "exobase"]
    c1 = [COcolor, O2color, H2color, O3color]  # CO, O2, H2, O3

    nom_i_py = Dict("exobase"=>meanTe, "tropopause"=>meanTt, "surface"=>meanTs)


    # where to place text, basically.
    mean_text = Dict("surface"=>[218, 100], "tropopause"=>[131, 4], "exobase"=>[208, 4])
    ylims = Dict("surface"=>[-10, 190], "tropopause"=>[-6.3, 4.2], "exobase"=>[-6.3, 4.2])
    xt_args = Dict("surface"=>150:20:280, "tropopause"=>100:10:160, "exobase"=>150:50:350)
    linelbls = DataFrame(Exp=["surface", "tropopause", "exobase"],
                          CO=[[155, -3],  [100, -5],   [175, -5]],
                          O2=[[220, -1.2],  [115, -1.8],   [150, -2.3]],
                          H2=[[180, 0.8],   [100, 2.3],  [150, 4.1]],
                          O3=[[155, 110], [100, -0.8], [155, -0.8]])

    # make plots pretty
    rcParams["font.family"] = "sans-serif"
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 18
    rcParams["axes.labelsize"]= 20
    rcParams["xtick.labelsize"] = 18
    rcParams["ytick.labelsize"] = 18

    fig, ax = subplots(1, 3, sharex=false, sharey=false, figsize=(21, 5))
    subplots_adjust(wspace=0.15)


    for i in 1:3
        plot_bg(ax[i])
        ex = exps[i]

        # calculate the relativeness of each point wrt data
        COdiff = (output_MR[ex]["CO"] .- data[1])/s[1]
        O2diff = (output_MR[ex]["O2"] .- data[2])/s[2]
        # H2diff = (ABSdictvar["H2"] .- data[3])/s[3]  # absolute abundance
        H2diff = (output_MR[ex]["H2MR"]./1e-6 .- data[3])/s[3] # ppm
        O3diff = (areadensity_to_micron_atm(output_abs[ex]["O3"]) .- data[4])/s[4]

        # plot
        ax[i].plot(tvals[ex], COdiff, marker="s", ms=sz, color=c1[1], zorder=10)
        ax[i].plot(tvals[ex], O2diff, marker="+", ms=sz, color=c1[2], zorder=10)
        ax[i].plot(tvals[ex], H2diff, marker="o", ms=sz, color=c1[3], zorder=10)
        ax[i].plot(tvals[ex], O3diff, marker="p", ms=sz, color=c1[4], zorder=10)
        # plot the mean values
        ax[i].axvline(nom_i_py[ex], color=medgray, zorder=5)
        ax[i].axhline(0, color="black")

        # text
        ax[i].text(mean_text[ex][1], mean_text[ex][2], "Standard\ntemperature", color=medgray, ha="left", va="top")
        dfentry = linelbls[linelbls.Exp.==ex, :]
        ax[i].text(dfentry.CO[1][1], dfentry.CO[1][2], "CO", color=c1[1], ha="left", va="top")
        ax[i].text(dfentry.O2[1][1], dfentry.O2[1][2], L"O$_2$", color=c1[2], ha="left", va="top")
        ax[i].text(dfentry.H2[1][1], dfentry.H2[1][2], L"H$_2$", color=c1[3], ha="left", va="top")
        ax[i].text(dfentry.O3[1][1], dfentry.O3[1][2], L"O$_3$", color=c1[4], ha="left", va="top")

        # labels and such
        ax[i].set_xlabel(ex*" temperature (K)")
        ax[i].set_xticks(xt_args[ex])
        ax[i].set_ylim(ylims[ex][1], ylims[ex][2])
    end


    # SURFACE PLOT CONFIG
    ax[1].xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(10))
    ax[1].set_yscale("symlog", linthresh=0.2)
    ax[1].set_yticks([-10, -1, 0, 1, 10, 100])
    ax[1].set_yticklabels([-10, -1, 0, 1, 10, 100])
    ax[1].set_ylabel(L"($X_{model}-X_{obs}$)/$\sigma$")

    # TROPOPAUSE PLOT CONFIG
    ax[2].xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(10))
    ax[2].yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(1))

    # Extra clarifying text
    ax[3].text(365, 4, "Model value \n"*L"$>$ observations", va="top")
    ax[3].text(365, 0, "Model value \n"*L"$\approx$observations", va="center")
    ax[3].text(365, -4, "Model value \n"*L"$<$observations")

    ax[1].text(270, 80, "a", color="black", fontsize=26, weight="bold", va="bottom")
    ax[2].text(160, 4, "b", color="black", fontsize=26, weight="bold", va="top")
    ax[3].text(350, 4, "c", color="black", fontsize=26, weight="bold", va="top")

    savefig(detailed_results_dir*"output_vs_data_temps.png", bbox_inches="tight", dpi=300)
    savefig(results_dir*"AllResultPlots/output_vs_data_temps.png", bbox_inches="tight", dpi=300)
end

function make_T_f_plot(output_MR)
    #=
    makes a 3-panel plot of the tradeoff of f.
    =#

    # make plots pretty
    rcParams["font.family"] = "sans-serif"
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 18
    rcParams["axes.labelsize"]= 20
    rcParams["xtick.labelsize"] = 18
    rcParams["ytick.labelsize"] = 18

    # Dictionaries for where to put things for each experiment
    exps = ["surface", "tropopause", "exobase"]
    gmean_txt_loc = Dict("surface"=>[218, 0.8], "tropopause"=>[131, 0.8], "exobase"=>[208, 0.8])
    xt_args = Dict("surface"=>[150, 20, 275], "tropopause"=>[70, 30, 165], "exobase"=>[150, 50, 350])
    nom_i_py = Dict("exobase"=>meanTe, "tropopause"=>meanTt, "surface"=>meanTs)
    ylims = Dict("tropopause"=>[1e-5, 1e0], "exobase"=>[1e-5, 1e0], "surface"=>[1e-5, 1e0])
    minor_mult = Dict("exobase"=>25, "tropopause"=>10, "surface"=>10)
    major_mult = Dict("exobase"=>50, "tropopause"=>30, "surface"=>20)

    # Make the figure
    fig, ax = subplots(1, 3, sharex=false, sharey=false, figsize=(21, 5))
    subplots_adjust(wspace=0.3)

    for i in 1:3
        plot_bg(ax[i])
        ex = exps[i]

        # plot
        ax[i].plot(tvals[ex], output_MR[ex]["f"], marker="o", ms=sz, color=f_main_color, zorder=10)
        ax[i].axvline(nom_i_py[ex], color=medgray, zorder=5)
        ax[i].text(gmean_txt_loc[ex][1], gmean_txt_loc[ex][2], "Standard\ntemperature",
                   color=medgray, ha="left", va="top")

        # axis format
        mpltic = pyimport("matplotlib.ticker")
        ax[i].set_xlabel(ex*" temperature (K)", fontsize=24)
        ax[i].xaxis.set_major_locator(mpltic.MultipleLocator(major_mult[ex]))
        ax[i].xaxis.set_minor_locator(mpltic.MultipleLocator(minor_mult[ex]))
        ax[i].set_xticks(xt_args[ex][1]:xt_args[ex][2]:xt_args[ex][3], minor=false)
        ax[i].tick_params(axis="y", labelcolor=f_main_color, which="both")
        # ax[i].tick_params(axis="x", which="both")

        ax[i].set_yscale("log")
        ax[i].set_ylim(ylims[ex][1], ylims[ex][2])
        ax[i].set_yticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1])

        # plot the secondary axis showing increase
        nomdict = make_std_atmo_dict("mr")
        fdiff = normalize_val(output_MR[ex]["f"], nomdict["f"])
        ax_2 = ax[i].twinx()
        if i==3
            ax_2.set_ylabel(L"$f$/$f(\overline{T}_{standard})$", color=f_sec_color, fontsize=24)
        end
        ax_2.plot(tvals[ex], fdiff, marker="o", ms=sz, color=f_sec_color)
        ax_2.tick_params(axis="y", labelcolor=f_sec_color, which="both")
        for side in ["top", "bottom", "left", "right"]
            ax_2.spines[side].set_visible(false)
        end
    end



    ax[1].set_ylabel(L"Fractionation factor $f$", color=f_main_color, fontsize=24)
    ax[1].text(150, 1, "a", color="black", fontsize=26, weight="bold", va="top")
    ax[2].text(100, 1, "b", color="black", fontsize=26, weight="bold", va="top")
    ax[3].text(150, 1, "c", color="black", fontsize=26, weight="bold", va="top")


    savefig(detailed_results_dir*"f_vs_temps.png", bbox_inches="tight", dpi=300)
    savefig(results_dir*"/AllResultPlots/"*"f_vs_temps.png", bbox_inches="tight", dpi=300)
end

function make_T_loss_contribs_plot(output_abs)
    #=
    Make a plot showing how much H, D, HD, and H2 contribute to loss of H, D.
    First plot: Six panels. Row: (1) H loss and (2) D loss. Column: surface, tropopause, exobase.
    Second plot: just exobase temperatures.
    =#

    # make plots pretty
    rcParams["font.family"] = "sans-serif"
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 18
    rcParams["axes.labelsize"]= 20
    rcParams["xtick.labelsize"] = 18
    rcParams["ytick.labelsize"] = 18

    mpltic = pyimport("matplotlib.ticker")

    # Dictionaries for where to put things for each experiment
    exps = ["surface", "tropopause", "exobase"]
    xt_args = Dict("surface"=>[150, 20, 275], "tropopause"=>[100, 10, 160], "exobase"=>[150, 50, 350])
    # ylims = Dict("tropopause"=>[1e-5, 1e0], "exobase"=>[1e-5, 1e0], "surface"=>[1e-5, 1e0])
    minor_mult = Dict("exobase"=>25, "tropopause"=>10, "surface"=>10)
    major_mult = Dict("exobase"=>50, "tropopause"=>30, "surface"=>20)

    # FIGURE: all experiments ==================================================================
    fig, ax = subplots(2, 3, sharex=false, sharey=true, figsize=(21, 10))
    subplots_adjust(wspace=0.05)

    for i in 1:3
        plot_bg(ax[1, i])
        plot_bg(ax[2, i])
        ex = exps[i]

        # plot the information
        # H loss
        total_H_flux = output_abs[ex]["Hflux"]
        HDfrac = output_abs[ex]["HD>Hloss"]./total_H_flux
        H2frac = output_abs[ex]["H2>Hloss"]./total_H_flux
        Hfrac = output_abs[ex]["H>Hloss"]./total_H_flux

        ax[1, i].fill_between(tvals[ex], 0, (1 .- Hfrac-H2frac-HDfrac), color="#dedeff", zorder=10, alpha=0.6)  # HD is smallest contributor
        ax[1, i].fill_between(tvals[ex], (1 .- Hfrac-H2frac), (1 .- Hfrac), color="#0068f9", zorder=10, alpha=0.6)  # H2 middle contributor
        ax[1, i].fill_between(tvals[ex], (1 .- Hfrac), 1, color="#949ffd", zorder=10, alpha=0.6)  # H biggest contributor
        ax[1, i].text(0.1, 0.8, "atomic H", color="#0068f9", transform=ax[1,i].transAxes, fontsize=24)

        # D Loss
        total_D_flux = output_abs[ex]["Dflux"]
        Dfrac = output_abs[ex]["D>Dloss"]./total_D_flux
        HDfrac = output_abs[ex]["HD>Dloss"]./total_D_flux
        println(output_abs[ex]["HD"])

        ax[2, i].fill_between(tvals[ex], (1 .- Dfrac), 1, color="#7ab19f", zorder=10, alpha=0.6)  # D contribution
        ax[2, i].fill_between(tvals[ex], 0, (1 .- Dfrac), color="#007f65", zorder=10, alpha=0.6) # HD contrib
        ax[2, i].text(0.1, 0.8, "atomic D", color="#007f65", transform=ax[2,i].transAxes, fontsize=24)

        # axis format
        ax[2, i].set_xlabel(ex*" temperature (K)", fontsize=24)
        ax[1, i].xaxis.set_major_locator(mpltic.MultipleLocator(major_mult[ex]))
        ax[1, i].xaxis.set_minor_locator(mpltic.MultipleLocator(minor_mult[ex]))
        ax[2, i].xaxis.set_major_locator(mpltic.MultipleLocator(major_mult[ex]))
        ax[2, i].xaxis.set_minor_locator(mpltic.MultipleLocator(minor_mult[ex]))
        ax[1, i].set_xticks(xt_args[ex][1]:xt_args[ex][2]:xt_args[ex][3], minor=false)
        ax[2, i].set_xticks(xt_args[ex][1]:xt_args[ex][2]:xt_args[ex][3], minor=false)
        ax[1, i].set_ylim(0.001, 1)
        ax[2, i].set_ylim(0.001, 1)
        ax[1, i].set_yscale("log")
        ax[2, i].set_yscale("log")
        ax[1, 1].set_ylabel("Fraction of total loss", fontsize=24)
        ax[2, 1].set_ylabel("Fraction of total loss", fontsize=24)
    end

    # plot panel labels
    ax[1, 3].text(0.8, 0.1, L"H_2", color="#0057e8", transform=ax[1,3].transAxes, va="top", fontsize=24)
    ax[2, 3].text(0.8, 0.1, "HD", color="#006e54", transform=ax[2,3].transAxes, va="top", fontsize=24)
    ax[1, 1].text(150, 1, "a", color="black", fontsize=26, weight="bold", va="top", zorder=15)
    ax[1, 2].text(100, 1, "b", color="black", fontsize=26, weight="bold", va="top", zorder=15)
    ax[1, 3].text(150, 1, "c", color="black", fontsize=26, weight="bold", va="top", zorder=15)
    ax[2, 1].text(150, 1, "d", color="black", fontsize=26, weight="bold", va="top", zorder=15)
    ax[2, 2].text(100, 1, "e", color="black", fontsize=26, weight="bold", va="top", zorder=15)
    ax[2, 3].text(150, 1, "f", color="black", fontsize=26, weight="bold", va="top", zorder=15)
    savefig(detailed_results_dir*"loss_contributions_all.png", bbox_inches="tight", dpi=300)

    # FIGURE: Just exobase =======================================================================
    fig, ax = subplots(1, 2, sharex=true, sharey=true, figsize=(14, 5))
    subplots_adjust(wspace=0.05)
    ex = "exobase"

    for i in 1:2
        plot_bg(ax[i])
        ax[i].tick_params(which="both", right=true)
        ax[i].xaxis.set_major_locator(mpltic.MultipleLocator(major_mult[ex]))
        ax[i].xaxis.set_minor_locator(mpltic.MultipleLocator(minor_mult[ex]))
        ax[i].set_xticks(xt_args[ex][1]:xt_args[ex][2]:xt_args[ex][3], minor=false)
        ax[i].yaxis.set_major_locator(mpltic.MultipleLocator(0.2))
        ax[i].yaxis.set_minor_locator(mpltic.MultipleLocator(0.05))
        ax[i].set_yticks(0:0.2:1, minor=false)
        ax[i].set_ylim(0, 1)
    end

    # H loss -----------------------------------------------------
    total_H_flux = output_abs[ex]["Hflux"]
    HDfrac = output_abs[ex]["HD>Hloss"]./total_H_flux
    H2frac = output_abs[ex]["H2>Hloss"]./total_H_flux
    Hfrac = output_abs[ex]["H>Hloss"]./total_H_flux

    ax[1].fill_between(tvals[ex], 0, (1 .- Hfrac-H2frac-HDfrac), color=HDcolor, zorder=10, alpha=0.6)  # HD
    ax[1].fill_between(tvals[ex], (1 .- Hfrac-H2frac), (1 .- Hfrac), color=H2color, zorder=10, alpha=0.6)  # H2
    ax[1].fill_between(tvals[ex], (1 .- Hfrac), 1, color=Hcolor, zorder=10, alpha=0.6)  # H
    ax[1].text(0.1, 0.8, "atomic H", color=Hcolor_dark, transform=ax[1].transAxes, fontsize=24, zorder=20)

    # D Loss -----------------------------------------------------
    total_D_flux = output_abs[ex]["Dflux"]
    Dfrac = output_abs[ex]["D>Dloss"]./total_D_flux
    HDfrac = output_abs[ex]["HD>Dloss"]./total_D_flux

    ax[2].fill_between(tvals[ex], (1 .- Dfrac), 1, color=Dcolor, zorder=10, alpha=0.6)  # D
    ax[2].fill_between(tvals[ex], 0, (1 .- Dfrac), color=HDcolor, zorder=10, alpha=0.6) # HD
    ax[2].text(0.1, 0.8, "atomic D", color=Dcolor, transform=ax[2].transAxes, fontsize=24, zorder=20)

    # axis format ------------------------------------------------
    mpltic = pyimport("matplotlib.ticker")
    ax[1].set_xlabel(ex*" temperature (K)", fontsize=24)
    ax[1].set_ylabel("Fraction of total loss", fontsize=24)


    # plot panel labels ------------------------------------------
    ax[1].text(0.8, 0.1, L"H_2", color=H2color_dark, transform=ax[1].transAxes, va="top", fontsize=24, zorder=20)
    ax[2].text(0.8, 0.1, "HD", color=HDcolor, transform=ax[2].transAxes, va="top", fontsize=24, zorder=20)
    ax[2].plot([335, 345], [0.06, 0.02], color=HDcolor)

    ax[1].text(0.01, 1, "a", color="black", fontsize=26, weight="bold", va="top", zorder=20, transform=ax[1].transAxes)
    ax[2].text(0.01, 1, "b", color="black", fontsize=26, weight="bold", va="top", zorder=20, transform=ax[2].transAxes)
    savefig(detailed_results_dir*"loss_contributions_by_species_exobase.png", bbox_inches="tight", dpi=300)
    savefig(detailed_results_dir*"loss_contributions_by_species_exobase.png", bbox_inches="tight", dpi=300)
end

# MAIN ==========================================================================

# Make the primary results plot (f with comparison to other studies) ------------
plot_results_caltech_together(results_dir*"MainCases/")


# Make the results plots for the detailed cases ----------------------------------
println("Enter a folder to use, no slashes (e.g. Research/Results/DetailedCases/<folder>): ")
append_me = readline(stdin)
detailed_results_dir = results_dir*"DetailedCases/"*append_me*"/"

makeplots = false   # whether to make the individiual small plots, one per species/measurable,
                    # that go with each experiment
other_deuterated = false  # set to true to look at minor D-carrying species. 
                          # WARNING: not implemented for main results plots, only for 
                          # the small supplemental plots!
write_new_files = false  # set to true if running for first time after new simulations

# 0. small supplemental plots (not used in paper)
# 1. f as a function of parameter (e.g. temperature or water vapor)
# 2. atomic/molecular H/D and fluxes compared across the simulations 
#    MR = some data stored as mixing ratio
#    abs = data stored as absolute quantities; both are used in results
# 3. Comparison of model output of CO, O2, O3, H2 to data in lit
# 4. Plots showing the contribution of molecular and atomic H, D to loss
#    Only the exobase case of this plot is used in the paper in Supplementary Info.

# WATER
println("Analyzing water model output, building dicts, making small plots")
water_data_mr = analyze_water("mr", other_deuterated, makeplots)
water_data_abs = analyze_water("abs", other_deuterated, makeplots)
println("Making results plots for the detailed water cases")
make_water_f_plot(water_data_mr)
make_water_Hspecies_plot(water_data_mr, "mr")
make_water_Hspecies_plot(water_data_abs, "abs")
make_water_output_vs_data(water_data_mr, water_data_abs)
make_water_loss_contribs_plot(water_data_abs)
println()


# O2  TODO: update this
# println()
# println("Analyzing O flux model output, building dicts, making supporting plots")
# o_data_abs = analyze_Oflux("abs", other_deuterated, makeplots)
# o_data_mr = analyze_Oflux("mr", other_deuterated, makeplots)


# TEMPERATURE
println("Analyzing temp model output, building dicts, making small plots")
T_data_mr = analyze_T("mr", other_deuterated, makeplots)
T_data_abs = analyze_T("abs", other_deuterated, makeplots)
println("Making results plots for the detailed temperature cases")
make_T_f_plot(T_data_mr)
make_T_Hspecies_plot(T_data_mr, "mr")
make_T_Hspecies_plot(T_data_abs, "abs")
make_T_output_vs_data(T_data_mr, T_data_abs)
make_T_loss_contribs_plot(T_data_abs)


# write files if being written for the first time
if write_new_files
    println("Writing jld storage files")
    wd_mr = jldopen(detailed_results_dir*"water_MR_data.jld", "w")
    @write wd_mr water_data_mr
    close(wd_mr)
    wd_abs = jldopen(detailed_results_dir*"water_abs_data.jld", "w")
    @write wd_abs water_data_abs
    close(wd_abs)

    # O_mr = jldopen(detailed_results_dir*"O_MR_data.jld", "w")
    # @write O_mr o_data_mr
    # close(O_mr)
    # O_abs = jldopen(detailed_results_dir*"O_abs_data.jld", "w")
    # @write O_abs o_data_abs
    # close(O_abs)

    T_mr = jldopen(detailed_results_dir*"T_MR_data.jld", "w")
    @write T_mr T_data_mr
    close(T_mr)
    T_abs = jldopen(detailed_results_dir*"T_abs_data.jld", "w")
    @write T_abs T_data_abs
    close(T_abs)
end
