################################################################################
# plot_profiles.jl
# TYPE: (1) Model files - required
# DESCRIPTION: plot the temperature profiles and water profile used in the model
# 
# Eryn Cangi
# Created May 2019
# Last edited: 21 July 2020
# Currently tested for Julia: 1.4.1
################################################################################

using PyCall
using PyPlot
using HDF5
using PlotUtils
using LaTeXStrings
using Analysis

include("PARAMETERS.jl")

function plot_T_6panel(base, profile_array)
    #=
    Plots a 6-panel figure of temperature profiles, the ones where there is 
    some perturbance from mean.

    base: Location of main results
    profile_array: list of the 6 temperature profiles other than the mean.
    =#
    

    titles = [L"$T_{surf}=0.75\overline{T}_{surf}$", L"$T_{tropo}=0.75\overline{T}_{tropo}$", 
              L"$T_{exo}=0.75\overline{T}_{exo}$", L"$T_{surf}=1.25\overline{T}_{surf}$", 
              L"$T_{tropo}=1.25\overline{T}_{tropo}$", L"$T_{exo}=1.25\overline{T}_{exo}$"]

    rcparams = PyCall.PyDict(matplotlib."rcParams")
    rcparams["font.sans-serif"] = ["Louis George Caf?"]
    rcparams["font.monospace"] = ["FreeMono"]
    rcparams["font.size"] = 18
    rcparams["axes.labelsize"]= 20
    rcparams["xtick.labelsize"] = 18
    rcparams["ytick.labelsize"] = 18

    alt = (0:2e5:zmax)

    fig, ax = subplots(2, 3, sharex=true, sharey=true, figsize=(10,8))
    subplots_adjust(wspace=0.05, hspace=0.4)

    i = 1
    j = 1
    modcol = "navy"
    meancol = "xkcd:cerulean blue"
    for (el, tit) in zip(profile_array, titles)
        ax[j, i].plot([Tpiecewise(a, meanTs, meanTt, meanTe) for a in alt], 
                      alt/1e5, color=meancol, linestyle="--", label="Standard profile")
        ax[j, i].plot([Tpiecewise(a, el[1], el[2], el[3]) for a in alt], 
                       alt/1e5, color=modcol, label="Alternate profiles") 
        ax[j, i].set_title(tit, fontsize=18)
        ax[j, i].set_xticks([100, 150, 200, 250])
        ax[j, i].tick_params(axis="x", labelbottom=true, labelsize=18)
        plot_bg(ax[j, i])

        if i==1
            ax[j, i].set_ylabel("Altitude (km)")
            ax[j, i].set_yticks([0, 50, 100, 150, 200, 250])
            ax[j, i].set_yticklabels([0, 50, 100, 150, 200, 250])
            if j==2
                ax[j, i].set_xlabel("Temperature (K)")
            end
        end
        
        i+=1
        if i == 4
            i = 1
            j = 2
        else
            continue
        end
    end

    ax[1, 1].text(meanTe+8, 220, L"T_{exo}")
    ax[1, 1].text(meanTt+7, 80, L"T_{tropo}")
    ax[1, 1].text(lowTs-65, -2, L"T_{surf-}")
    ax[1, 1].plot(meanTe, 250, marker="o", color=modcol)
    ax[1, 1].plot(meanTt, 80, marker="o", color=modcol)
    ax[1, 1].plot(lowTs, 0, marker="o", color=modcol)
    legend(bbox_to_anchor=(1.05, -0.07), fontsize=16)
    
    savefig(base*"ALL STUDY PLOTS/temp_profiles.png", bbox_inches="tight")
    savefig(base*"VarWaterTemp/temp_profiles.png", bbox_inches="tight")
end

function plot_one_profile(Tp, titletext, base)
    #=
    Tp: a list of the temps [T_surface, T_tropo, T_exo]
    titletext: exactly what it sounds like
    base: where to save the figure

    This is similar to plot_temp_prof but more useful for plotting a profile 
    during analysis.
    =#
    # make plots pretty
    rcParams = PyCall.PyDict(matplotlib."rcParams")
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 18
    rcParams["axes.labelsize"]= 20
    rcParams["xtick.labelsize"] = 18
    rcParams["ytick.labelsize"] = 18

    alt = (0:dz:zmax)

    # plot
    fig, ax = subplots(figsize=(4,6))
    plot_bg(ax)
    bloo = "xkcd:cerulean blue"
    plot([Tpiecewise(a, Tp[1], Tp[2], Tp[3]) for a in alt], alt/1e5, color=bloo)
    ax.plot(meanTe, 250, marker="o", color=bloo)
    ax.plot(meanTt, 85, marker="o", color=bloo)
    ax.plot(meanTs, 0, marker="o", color=bloo)

    # text
    ax.text(Tp[1]-12, 250, L"\overline{T}_{exo}="*"$(Int(Tp[3]))K", ha="right", va="center")
    ax.text(Tp[2]+5, 85, L"\overline{T}_{tropo}="*"$(Int(Tp[2]))K", va="center")
    ax.text(Tp[3]-15, 0, L"\overline{T}_{surf}="*"$(Int(Tp[1]))K", ha="right", va="center")

    # ticks and labels 
    ax.set_xlabel("Temperature (K)")
    ax.set_xticks([125,150,175,200,225])
    ax.set_ylabel("Altitude (km)")
    ax.set_yticks([0, 50, 100, 150, 200, 250])
    ax.set_yticklabels([0,50,100,150,200,250])
    title(titletext)

    # save it 
    savename = "temp_profile_"*join(split(titletext), "_")*".png"
    savefig(base*"ALL STUDY PLOTS/"*savename, bbox_inches="tight")
    savefig(base*"VarWaterTemp/"*savename, bbox_inches="tight")
end

function plot_all_water_profs(Tp, hygropause_alt, base, resultsfolder; plot_HDO=true)
    #=
    Plots the various water profiles all together. Unfortunately some text has 
    to be copied directly from converge_new_file.jl to set up the water profiles.

    Tp: a vector of the three temperatures
    hygropause_alt: altitude of the hygropause in cm
    base: Results folder (general)
    resultsfolder: actual folder containing results, located within "base"
    plot_HDO: whether to plot HDO on the plot. default true.
    =#

    Temp(z::Float64) = Tpiecewise(z, Tp[1], Tp[2], Tp[3])
    n_current = get_ncurrent("converged_250km_atmosphere.h5")

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
    n_current[:JHDO2toHDOpO1D] =  zeros(length(alt))
    n_current[:JODtoO1DpD]  =  zeros(length(alt))

    # Create the saturation profile 
    H2Osat = map(x->Psat(x), map(Temp, alt)) # array in #/cm^3 by altitude
    H2Osatfrac = H2Osat./map(z->n_tot(n_current, z), alt)  # get SVP as fraction of total atmo

    # Make the plot 
    fig, ax = subplots(figsize=(7,5))
    plot_bg(ax)
    rcParams = PyDict(matplotlib."rcParams")
    rcParams["font.size"] = 16
    rcParams["axes.labelsize"]= 20
    rcParams["xtick.labelsize"] = 18
    rcParams["ytick.labelsize"] = 18
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    fs = Dict("smtxt"=>14, "lbl"=>20, "tick"=>18, "mid"=>16)

    waterfolders = search_subfolders(base*resultsfolder, r"water_\d.+")
    waters = sort([parse(Float64, match(r"\d.+", w).match) for w in waterfolders])
    cols = ["#74a9cf", "#3690c0", "#0570b0","#034e7b", "#013554"] # colors for each water profile
    prum = [100, 50, 25, 10, 1] # It's backwards because the colors are 
    j = 1
    for MR in waters
        # set H2O SVP fraction to minimum for all alts above first time min is reached
        H2Oinitfrac = H2Osatfrac[1:something(findfirst(isequal(minimum(H2Osatfrac)), H2Osatfrac), 0)]
        H2Oinitfrac = [H2Oinitfrac;   # ensures no supersaturation
                       fill(minimum(H2Osatfrac), length(alt)-2-length(H2Oinitfrac))]

        # make profile constant in the lower atmosphere (well-mixed)
        H2Oinitfrac[findall(x->x<hygropause_alt, alt)] .= MR
 
        # restrict the initial water fraction to be below the saturation vapor pressure curve
        for i in [1:length(H2Oinitfrac);]
            H2Oinitfrac[i] = H2Oinitfrac[i] < H2Osatfrac[i+1] ? H2Oinitfrac[i] : H2Osatfrac[i+1]
        end
        
        # Plot it 
        ax.semilogx(H2Oinitfrac, alt[2:end-1]./1e5, color=cols[j], linewidth=2)
        if plot_HDO
            HDOinitfrac = H2Oinitfrac * DH  # initial profile for HDO
            ax.semilogx(HDOinitfrac, alt[2:end-1]./1e5, color=cols[j], linewidth=2, linestyle=":")
        end
        j += 1
    end

    if plot_HDO
        withHDO = "_with_HDO"
        ax.set_xlim(left=1e-11, right=1e-2)
        ax.set_xticks([1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2])
        text_x = 1e-11
        ax.text(10.0^(-7.9), -20, "A", color=cols[1], ha="center", fontsize=fs["smtxt"])
        ax.text(10.0^(-6.8), -20, "B", color=cols[2], ha="right", fontsize=fs["smtxt"])
        ax.text(10.0^(-6.3), -20, "C", color=cols[3], ha="right", fontsize=fs["smtxt"])
        ax.text(10.0^(-6.1), -20, "D", color=cols[4], ha="center", fontsize=fs["smtxt"])
        ax.text(10.0^(-5.8), -20, "E", color=cols[5], fontsize=fs["smtxt"])
        ax.text(10^(-10.5), 150, "HDO", color=medgray)
    else
        withHDO = ""
        text_x = 1e-7
        ax.set_xlim(left=1e-8, right=1e-2)
        ax.set_xticks([1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2])
    end
    
    # Place text indicating profiles and other annotations
    ax.text(10.0^(-4.7), -20, "A", color=cols[1], ha="right", fontsize=fs["smtxt"])
    ax.text(10.0^(-3.7), -20, "B", color=cols[2], ha="right", fontsize=fs["smtxt"])
    ax.text(10.0^(-3.2), -20, "C", color=cols[3], ha="right", fontsize=fs["smtxt"])
    ax.text(10.0^(-3), -20, "D", color=cols[4], ha="center", fontsize=fs["smtxt"])
    ax.text(10.0^(-2.8), -20, "E", color=cols[5], fontsize=fs["smtxt"])
    ax.text(text_x, -20, "Water profile:", fontsize=fs["mid"], color=medgray)
    ax.text(1e-6, 23, "Saturation curve", color=medgray, rotation=-15)
    ax.text(10^(-7.5), 150, "H"*L"_2"*"O", color=medgray)

    # formatting
    ax.set_xlabel("Water mixing ratio", fontsize=fs["lbl"])
    ax.set_ylabel("Altitude (km)", fontsize=fs["lbl"])
    ax.set_ylim(-25, 255)
    ax.set_yticks([0,50,100,150,200,250])

    ax.tick_params(which="both", labelsize=fs["tick"])
    for lbl in ax.xaxis.get_ticklabels()[2:2:end]
        lbl.set_visible(false)
    end
    
    # save
    savefig(base*"ALL STUDY PLOTS/water_profiles"*withHDO*".png", bbox_inches="tight")
    savefig(base*resultsfolder*"water_profiles"*withHDO*".png", bbox_inches="tight")
end

function plot_temp_3panel(base)
    #=
    base: main results folder (general)

    Makes a 3-subfigure plot of all the temperature profiles used for the 
    detailed temperature variation experiments.
    =#

    rcParams = PyCall.PyDict(matplotlib."rcParams")
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 18
    rcParams["axes.labelsize"]= 20
    rcParams["xtick.labelsize"] = 18
    rcParams["ytick.labelsize"] = 18

    Ts = collect(150:10:270)
    Tt = collect(lowTt:10:hiTt)
    Te = collect(150:25:350)
    
    numsurf = length(Ts)
    numtropo = length(Tt)
    numexo = length(Te)
    totallen = length(Ts) + length(Tt) + length(Te)

    fig, ax = subplots(1, 3, figsize=(15, 3))
    for axob in ax
        axob.set_xlabel("Temperature (K)")
        axob.set_yticks(collect(0:50:zmax/1e5))
        axob.tick_params(axis="x", labelrotation=45)
        plot_bg(axob)
    end
  
    # Surface axis
    ax[1].set_xticks(collect(meanTt:20:hiTs))
    ax[1].set_xticklabels(collect(meanTt:20:hiTs))
    ax[1].set_ylabel("Altitude (km)")
    ax[1].set_title(L"Varying T$_{surf}$")

    # Tropopause axis
    ax[2].set_xticks(collect(lowTt:20:220))
    ax[2].set_xticklabels(collect(lowTt:20:220))
    ax[2].set_title(L"Varying T$_{tropo}$")

    # Exobase axis 
    ax[3].set_xticks(collect(125:25:350))
    ax[3].set_xticklabels(collect(125:25:350))
    ax[3].set_title(L"Varying T$_{exo}$")

    for k in range(1, length=totallen)
        cols = get_colors(numsurf, "plasma")
        if k <= numsurf  # surface
            Tprof = map(h->Tpiecewise(h, Ts[k], meanTt, meanTe, "surf"), alt)
            ax[1].plot(Tprof, alt./1e5, color=cols[k, :])
        end
        
        cols = get_colors(numtropo, "plasma")
        if numsurf+1 <= k <= numsurf+numtropo  # mesosphere
            Tprof = map(h->Tpiecewise(h, meanTs, Tt[k-numsurf], meanTe, "tropo"), alt)
            ax[2].plot(Tprof, alt./1e5, color=cols[k-numsurf, :])
        end
        
        cols = get_colors(numexo, "plasma")
        if numsurf+numtropo+1 <= k  # exobase
            Tprof = map(h->Tpiecewise(h, meanTs, meanTt, Te[k-(numsurf+numtropo)], "exo"), alt)
            ax[3].plot(Tprof, alt./1e5, color=cols[k-(numsurf+numtropo), :])
        end
    end
    savefig(base*"ALL STUDY PLOTS/tradeoff_temp_profiles.png", bbox_inches="tight")
    savefig(base*"TradeoffPlots/Main Results/tradeoff_temp_profiles.png", bbox_inches="tight")
end


Tprofs = [[lowTs, meanTt, meanTe], [meanTs, lowTt, meanTe], [meanTs, meanTt, lowTe], 
                 [hiTs, meanTt, meanTe], [meanTs, hiTt, meanTe], [meanTs, meanTt, hiTe]]
plot_T_6panel(results_dir, Tprofs)

plot_temp_3panel(results_dir)

plot_one_profile([meanTs, meanTt, meanTe], "Standard temperature", results_dir)

plot_all_water_profs(meantemps, hygropause_alt, results_dir, "MainCases/", plot_HDO=true)