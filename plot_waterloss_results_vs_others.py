################################################################################
# plot_waterloss_results_vs_others.py
# TYPE: (2) Analysis - required
# DESCRIPTION: Simple script to make a plot of water loss results vs others
# It's in python. sorry. I was lazy and the plotting methods I use in Julia just
# call matplotlib anyway. So when I can, I just use python.
#
# Eryn Cangi
# Created 23 October 2019
# Last edited: 21 July 2020
# Currently tested for Python: 3+
################################################################################

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
plt.style.use('default')
plt.rc('text', usetex=False)
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Louis George Caf?'
plt.rcParams['font.monospace'] = 'FreeMono'
plt.rcParams['font.size'] = 18
plt.rcParams['axes.labelsize'] = 22
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18

research_dir = "/home/emc/GDrive-CU/Research-Modeling/FractionationFactor/Code/"
results_dir = research_dir+"Results/"

studies = ["Yung+ 1988", "Krasnopolsky 2000", "Lammer+ 2003", "Kurokawa+ 2014", "Villanueva+ 2015", "Alsaeed+ 2019", "This work"]
dummyind = [0, 1, 2, 3, 4, 5, 6]
my_water = [66, 122]  # TODO: fill in my water results here.
water = [[3.4, 3.4], [65, 120], [14, 34], [51, 152], [137, 165], [20, 220], my_water]
col = ["#10007A", "#10007A", "#2F7965", "#2F7965", "#2F7965", "#e16262", "#10007A"]  
lb = ["1D photochemical modeling", "1D photochemical modeling",
      "Observations", "Observations", "Observations",
      "1D box model", "1D photochemical modeling"]

fig, ax = plt.subplots(figsize=(7, 6))
ax.set_facecolor("#ededed")
ax.grid(zorder=0, color="white", which="major")
for side in ["top", "bottom", "left", "right"]:
    ax.spines[side].set_visible(False)


for d, w in zip(dummyind, water):
    ax.plot(w, [d, d], linewidth=15, color=col[d])

p1 = mpatches.Patch(color="#10007A", label="1D photochemical modeling")
p2 = mpatches.Patch(color="#2F7965", label="observations")
p3 = mpatches.Patch(color="#e16262", label="1D box model")

plt.legend(handles=[p1, p2, p3], loc=[0.52, 0.01], fontsize=12)#bbox_to_anchor=(1.01, 1), 

ax.scatter(3.4, 0, marker="|", s=150, zorder=20, color="#10007A")
plt.yticks(dummyind, labels=studies)
plt.title("Study comparison of lost water")
plt.xlabel("Water lost (m GEL)")

plt.savefig(results_dir+"AllResultPlots/water_loss_comparison.png", bbox_inches="tight", dpi=300)
plt.savefig(results_dir+"VarWaterTemp/water_loss_comparison.png", bbox_inches="tight", dpi=300)