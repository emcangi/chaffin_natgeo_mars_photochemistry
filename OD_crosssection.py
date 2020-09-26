################################################################################
# HDO_crosssection.py
# TYPE: (0) Setup files
# DESCRIPTION: Takes a bunch of data from Cheng+1999 and Cheng+2004 to get HDO
# cross sections by wavelength. The data is not complete for the wavelengths 
# that we need, so extrapolation and fudge factors are employed to make it work.
#
# Eryn Cangi
# Finalized 11 October 2019
# Currently tested for Python: 3.7
################################################################################

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams.update({'font.size': 16})
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

# load the OH cross sections
oh = np.loadtxt(open(research_dir+"uvxsect/binnedOH.csv", "rb"), delimiter=",", skiprows=4)

# load what I got for OD
od = np.loadtxt(research_dir+"uvxsect/OD with whole nums.dat")

c1 = "cornflowerblue"
c2 = "xkcd:bright orange"

fig = plt.figure(figsize=(10, 8))
plt.semilogy(oh[:, 0], oh[:, 1], color=c1, label="OH (Barfield 1972)")
plt.semilogy(od[:, 0], od[:, 1], color=c2, label="OD (Nee & Lee 1984)")
ax = plt.gca()

ax.set_facecolor("#ededed")
ax.grid(zorder=0, color="white", which="both")
for side in ["top", "bottom", "left", "right"]:
    ax.spines[side].set_visible(False)

#ax.axvline(145)
#ax.set_xticks(np.arange(120, 220, 10))
plt.title("OH and OD cross sections", y=1.02)
plt.ylabel("Cross section (cm^2)")
plt.xlabel("Wavelength (nm)")
plt.legend()
plt.savefig(results_dir+"AllResultPlots/ODxsects.png", bbox_inches="tight", dpi=300)