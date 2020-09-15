################################################################################
# plot_flux_contributions.jl
# TYPE: (2) Analysis files - optional
# DESCRIPTION: plots the total thermal and non-thermal flux, as well as 
# contributions from each H- or D-bearing species, to understand how tempreature
# at tropopause or exobase affects f.
# 
# Eryn Cangi
# Created 26 August 2020
# Last edited: 
# Currently tested for Julia: 1.4.1
################################################################################

using Analysis
using PyPlot
using PyCall
using PlotUtils
using HDF5, JLD
using DataFrames

include("PARAMETERS.jl")

function make_figure(temps, dict_t, dict_nt, exptype)
	#=
	Makes a figure showing flux quantities for loss of H and D.
	Plots total loss in thermal and thermal+nonthermal escape, contributions
	from H, H2, HD, and D, and the fractionation factor for comparison.

	temps: An array of the values for one temperature parameter (i.e. T_tropo) for which to plot.
	       i.e. if studying the tropopause temperature experiments, this array should run from
	       100:160.
    dict_t: Dictionary with three keys: H, D, and f. each entry for H and D is also a dictionary of 
            flux values for thermal escape across all the experiments in temps. f is just a list of
            fractionation factor values across ecah temperature experiment.
    dict_nt: the same as dict_t but for non-thermal escape only.
    exptype: either "tropopause" or "exobase", just used for filling in the saved filename.

    output: a three-panel plot .png.
	=#

	# make plots pretty
	rcParams["font.family"] = "sans-serif"
	rcParams["font.sans-serif"] = ["Louis George Caf?"]
	rcParams["font.monospace"] = ["FreeMono"]
	rcParams["font.size"] = 16
	rcParams["axes.labelsize"]= 18
	rcParams["xtick.labelsize"] = 14
	rcParams["ytick.labelsize"] = 14

	totcol = "red"
	atomiccol = "blue"
	moleccol = "xkcd:kelly green"

	thermstyle = "--"
	ntstyle = ":"
	tntstyle = "-"

	linelbltxt = Dict("Ht"=>"H thermal",
	                  "Hnt"=>"H non-thermal",
	                  "Dt"=>"D thermal",
	                  "Dnt"=>"D non-thermal",
	                  "H2t"=>L"H$_2$ thermal",
	                  "H2nt"=>L"H$_2$ non-thermal",
	                  "HDt"=>"HD thermal",
	                  "HDnt"=>"HD non-thermal",
	                  "totalt"=>"total thermal",
	                  "totaltnt"=>"total thermal+non-thermal",
	                  "ft"=>"f thermal",
	                  "ftnt"=>"f thermal+non-thermal")

	linelbls = DataFrame(Exp=["tropopause", "exobase"],
                          Ht=[[102, 2e8],  [178, 1.8e8]],
                          Hnt=[[102, 3.7e7],  [210, 3e7]],
                          H2t=[[100, 1e6], [160, 1e4]],
                          H2nt=[[100, 8e6], [155, 9e6]],
                          Dt=[[122, 1.2e2],  [180, 2e1]],
                          Dnt=[[142, 3e3],  [215, 5e3]],
                          HDt=[[141, 3e-1], [215, 2e-1]],
                          HDnt=[[135, 1e3], [148, 1.5e3]],
                          totaltH=[[135, 2e8], [205, 1.8e8]],
                          totaltntH=[[135, 3.7e8],  [200, 5e8]],
                          totaltD=[[100, 7e1], [152, 4e0]],
                          totaltntD=[[105, 1.5e4],  [150, 3.1e4]],
                          ft=[[140, 4e-3], [150, 2e-4]],
                          ftnt=[[100, 3e-2], [205, 1.3e-1]])

	# H fig ==========================================
	fig, ax = subplots(2, 1, sharex=true, sharey=false, figsize=(7, 12))
	subplots_adjust(hspace=0.1, wspace=1)
	s = "H"

	# H species panel -----------------------
	plot_bg(ax[1])
	ax[1].plot(temps, dict_nt[s]["total"]+dict_t[s]["total"], marker="o", color=totcol, zorder=10, linestyle=tntstyle, label="Total, t+nt")
	ax[1].plot(temps, dict_t[s]["total"], marker="o", color=totcol, zorder=10, linestyle=thermstyle, label="Total, thermal")

	# H
	ax[1].plot(temps, dict_t[s][:H], marker="x", color=atomiccol, zorder=10, linestyle=thermstyle, label="H contrib.")
	ax[1].plot(temps, dict_nt[s][:H], marker="x", color=atomiccol, zorder=10, linestyle=ntstyle) # note this is non-thermal only!

	# H2
	ax[1].plot(temps, dict_t[s][:H2], marker="D", color=moleccol, zorder=10, linestyle=thermstyle, label="H2 contrib.")
	ax[1].plot(temps, dict_nt[s][:H2], marker="D", color=moleccol, zorder=10, linestyle=ntstyle) # note this is non-thermal only!

	# line labels 
	dfentry = linelbls[linelbls.Exp.==exptype, :]
	ax[1].text(dfentry.Ht[1][1], dfentry.Ht[1][2], linelbltxt["Ht"], color=atomiccol, ha="left", va="top")
	ax[1].text(dfentry.Hnt[1][1], dfentry.Hnt[1][2], linelbltxt["Hnt"], color=atomiccol, ha="left", va="top")
    ax[1].text(dfentry.H2t[1][1], dfentry.H2t[1][2], linelbltxt["H2t"], color=moleccol, ha="left", va="top")
	ax[1].text(dfentry.H2nt[1][1], dfentry.H2nt[1][2], linelbltxt["H2nt"], color=moleccol, ha="left", va="top")
    ax[1].text(dfentry.totaltH[1][1], dfentry.totaltH[1][2], linelbltxt["totalt"], color=totcol, ha="left", va="top")
    ax[1].text(dfentry.totaltntH[1][1], dfentry.totaltntH[1][2], linelbltxt["totaltnt"], color=totcol, ha="left", va="top")

	# axis format
	ax[1].set_yscale("log")
	ax[1].set_title("$(s) loss")#("$(s) fluxes, total")
	ax[1].set_xticks(temps)
	ax[1].set_ylabel(L"Flux (# cm$^{-2}$s$^{-1}$)")

	# D species panel -------------------------------------
	s = "D"
	plot_bg(ax[2])
	ax[2].set_xticks(temps)

	# totals
	ax[2].plot(temps, dict_t[s]["total"], marker="o", color=totcol, linestyle=thermstyle, label="Total, thermal ")
	ax[2].plot(temps, dict_nt[s]["total"]+dict_t[s]["total"], marker="o", color=totcol, linestyle=tntstyle, label="Total, t+nt ")   # thermal + nonthermal total escape!

	# D
	ax[2].plot(temps, dict_t[s][:D], marker="x", color=atomiccol, linestyle=thermstyle, label="D")
	ax[2].plot(temps, dict_nt[s][:D], marker="x", color=atomiccol, linestyle=ntstyle) # nonthermal only!

	# HD
	ax[2].plot(temps, dict_t[s][:HD], marker="D", color=moleccol, linestyle=thermstyle, label="HD")
	ax[2].plot(temps, dict_nt[s][:HD], marker="D", color=moleccol, linestyle=ntstyle) # nonthermal only!

	# separate axis for fractionation factor
	ax4 = ax[2].twinx()
	ax4.set_yscale("log")
	ax4.set_ylabel("fractionation factor", color=medgray)
	ax4.tick_params(axis="y", labelcolor=medgray, which="both")
	for side in ["top", "bottom", "left", "right"]
	    ax4.spines[side].set_visible(false)
	end
	ax4.plot(temps, dict_t["f"], marker="o", color=medgray, zorder=5, linestyle=thermstyle, label="f thermal")
	ax4.plot(temps, dict_nt["f"], marker="o", color=medgray, zorder=5, linestyle=tntstyle, label="f both")
	
	# axis labels and such 
	ax[2].set_title("$(s) loss")
	ax[2].set_yscale("log")
	ax[2].set_xlabel("$(exptype) temperature (K)")
	ax[2].set_ylabel(L"Flux (# cm$^{-2}$s$^{-1}$)")

	# line labels 
	dfentry = linelbls[linelbls.Exp.==exptype, :]

	ax[2].text(dfentry.Dt[1][1], dfentry.Dt[1][2], linelbltxt["Dt"], color=atomiccol, ha="left", va="top")
	ax[2].text(dfentry.Dnt[1][1], dfentry.Dnt[1][2], linelbltxt["Dnt"], color=atomiccol, ha="left", va="top")
    ax[2].text(dfentry.HDt[1][1], dfentry.HDt[1][2], linelbltxt["HDt"], color=moleccol, ha="left", va="top")
	ax[2].text(dfentry.HDnt[1][1], dfentry.HDnt[1][2], linelbltxt["HDnt"], color=moleccol, ha="left", va="top")
    ax[2].text(dfentry.totaltD[1][1], dfentry.totaltD[1][2], linelbltxt["totalt"], color=totcol, ha="left", va="top")
    ax[2].text(dfentry.totaltntD[1][1], dfentry.totaltntD[1][2], linelbltxt["totaltnt"], color=totcol, ha="left", va="top")
    ax4.text(dfentry.ft[1][1], dfentry.ft[1][2], linelbltxt["ft"], color=medgray, ha="left", va="top")
    ax4.text(dfentry.ftnt[1][1], dfentry.ftnt[1][2], linelbltxt["ftnt"], color=medgray, ha="left", va="top")
	
	savefig(basepath*"$(exptype)_contributions.png")
	savefig(results_dir*"/ALL STUDY PLOTS/"*"$(exptype)_contributions.png")
end

# ===============================================================================

basepath = results_dir * det_cases_dir

# Temperature ranges which can be plotted (i.e. temperatures for which we can estimate the nonthermal flux)
Tt = [100, 110, 120, 130, 140, 150, 160]
Te = [150, 205, 250]

# Dictionaries to store the data across simulations
tropo_t = Dict("H"=>Dict("total"=>[], :H=>[], :H2=>[], :HD=>[]), 
	            "D"=>Dict("total"=>[], :D=>[], :HD=>[]),
	            "f"=>[])
tropo_t_nt = Dict("H"=>Dict("total"=>[], :H=>[], :H2=>[], :HD=>[]),
			   "D"=>Dict("total"=>[], :D=>[], :HD=>[]),
			   "f"=>[])

exo_t = Dict("H"=>Dict("total"=>[], :H=>[], :H2=>[], :HD=>[]), 
	            "D"=>Dict("total"=>[], :D=>[], :HD=>[]),
	            "f"=>[])
exo_t_nt = Dict("H"=>Dict("total"=>[], :H=>[], :H2=>[], :HD=>[]),
			   "D"=>Dict("total"=>[], :D=>[], :HD=>[]),
			   "f"=>[])

# These sections open the file and collect the fluxes and contributions from each species for thermal and non-thermal escape.
# The first section is for the tropopause temperatures.
for t in Tt
	f_t, f_nt, c_t, c_nt = get_flux(:H, basepath*"temp_216_$(t)_205/converged_temp_216_$(t)_205.h5", 1.2e8, [216., Float64(t), 205.], repro=false, therm_only=false)
	append!(tropo_t["H"]["total"], f_t)   # Total H thermal ecsape flux from H, H2, and HD species
	append!(tropo_t["H"][:H], c_t[:H])    # Contribution of atomic H loss to thermal escape of H
	append!(tropo_t["H"][:H2], c_t[:H2])  # Contribution of molecular H
	append!(tropo_t["H"][:HD], c_t[:HD])  # Contribution of molecular hydrogen deuteride
	append!(tropo_t_nt["H"]["total"], f_nt)   # Total H non-thermal escape flux from H, H2, HD. Non-thermal ONLY.
	append!(tropo_t_nt["H"][:H], c_nt[:H])    # Contribution of atomic H loss to non-thermal escape of H
	append!(tropo_t_nt["H"][:H2], c_nt[:H2])  # Contrib. of molec. H
	append!(tropo_t_nt["H"][:HD], c_nt[:HD])  # Contrib. of molec. H deuteride

	# This is the same information but for loss of D via D and HD.
	f_t, f_nt, c_t, c_nt = get_flux(:D, basepath*"temp_216_$(t)_205/converged_temp_216_$(t)_205.h5", 1.2e8, [216., Float64(t), 205.], repro=false, therm_only=false)
	append!(tropo_t["D"]["total"], f_t)
	append!(tropo_t["D"][:D], c_t[:D])
	append!(tropo_t["D"][:HD], c_t[:HD])
	append!(tropo_t_nt["D"]["total"], f_nt) # non-thermal escape flux, NOT thermal + nonthermal.
	append!(tropo_t_nt["D"][:D], c_nt[:D])
	append!(tropo_t_nt["D"][:HD], c_nt[:HD])

	# fractionation factor. Here we retrieve f for thermal only and thermal+non-thermal as in the main results.
	append!(tropo_t["f"], calculate_f(basepath*"temp_216_$(t)_205/converged_temp_216_$(t)_205.h5", "thermal", [216., Float64(t), 205.], 1.2e8))
	append!(tropo_t_nt["f"], calculate_f(basepath*"temp_216_$(t)_205/converged_temp_216_$(t)_205.h5", "both", [216., Float64(t), 205.], 1.2e8))
end

# same thing but for exobase temperatures.
for t in Te
	f_t, f_nt, c_t, c_nt = get_flux(:H, basepath*"temp_216_130_$(t)/converged_temp_216_130_$(t).h5", 1.2e8, [216., 130., Float64(t)], repro=false, therm_only=false)
	append!(exo_t["H"]["total"], f_t)
	append!(exo_t["H"][:H], c_t[:H])
	append!(exo_t["H"][:H2], c_t[:H2])
	append!(exo_t["H"][:HD], c_t[:HD])
	append!(exo_t_nt["H"]["total"], f_nt) # note this is the total non-thermal escape flux, NOT thermal + nonthermal.
	append!(exo_t_nt["H"][:H], c_nt[:H])
	append!(exo_t_nt["H"][:H2], c_nt[:H2])
	append!(exo_t_nt["H"][:HD], c_nt[:HD])

	f_t, f_nt, c_t, c_nt = get_flux(:D, basepath*"temp_216_130_$(t)/converged_temp_216_130_$(t).h5", 1.2e8, [216., 130., Float64(t)], repro=false, therm_only=false)
	append!(exo_t["D"]["total"], f_t)
	append!(exo_t["D"][:D], c_t[:D])
	append!(exo_t["D"][:HD], c_t[:HD])
	append!(exo_t_nt["D"]["total"], f_nt) # note this is the total non-thermal escape flux, NOT thermal + nonthermal.
	append!(exo_t_nt["D"][:D], c_nt[:D])
	append!(exo_t_nt["D"][:HD], c_nt[:HD])

	# fractiontion factor
	append!(exo_t["f"], calculate_f(basepath*"temp_216_130_$(t)/converged_temp_216_130_$(t).h5", "thermal", [216., 130., Float64(t)], 1.2e8))
	append!(exo_t_nt["f"], calculate_f(basepath*"temp_216_130_$(t)/converged_temp_216_130_$(t).h5", "both", [216., 130., Float64(t)], 1.2e8))
end


make_figure(Tt, tropo_t, tropo_t_nt, "tropopause")
make_figure(Te, exo_t, exo_t_nt, "exobase")