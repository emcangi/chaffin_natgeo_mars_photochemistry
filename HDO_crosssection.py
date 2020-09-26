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
from matplotlib.lines import Line2D

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

# Fit the 250K data on the red edge and extrapolate to 220 nm ==================

# Do the line fit
data250 = np.loadtxt("uvxsect/HDO_250K.dat")
wvs_to_fit = data250[-30:][:, 0]
xsect_to_fit = data250[-30:][:, 1]
xsect_to_fit_log = np.log10(xsect_to_fit)
z = np.polyfit(wvs_to_fit, xsect_to_fit_log, 1)
f = np.poly1d(z)
tryfit = 10**f(wvs_to_fit)

# Plot the original HDO data with the attempted line fit
# fig = plt.figure(figsize=(10, 8))
# plt.semilogy(data250[:, 0], data250[:, 2], color="cornflowerblue", label="HDO orig")
# plt.semilogy(wvs_to_fit, tryfit, color="blue", label="HDO try fit")
# plt.legend()
# plt.title("250K HDO tryfit", y=1.08)
# plt.ylabel("Cross section (cm^2)")
# plt.xlabel("Wavelength (nm)")
# plt.show()

#extrapolate to longer wavelengths.
newwv_250 = np.arange(182, 220, 0.2)
new_xsects_250 = f(newwv_250)
new_xsects_250 = 10**new_xsects_250
combined_wvs_250K = np.hstack((data250[:, 0], newwv_250))
combined_xsects_250K = np.hstack((data250[:, 1], new_xsects_250))

# Plot the extrapolated HDO cross sections
# fig = plt.figure(figsize=(10, 8))
# plt.semilogy(data250[:, 0], data250[:, 2], color="cornflowerblue", label="HDO orig")
# plt.semilogy(newwv_250, new_xsects_250, color="blue", label="HDO extrap")
# plt.semilogy(combined_wvs_250K, combined_xsects_250K, color="red", alpha=0.3, linewidth=8, label="put into one array")
# plt.legend()
# plt.title("250K HDO extrapolated by me", y=1.08)
# plt.ylabel("Cross section (cm^2)")
# plt.xlabel("Wavelength (nm)")
# plt.show()

# Now work on 298 K ============================================================
data298 = np.loadtxt("uvxsect/HDO_298K.dat")
data298 = data298[:-5]  # get rid of the empty rows

# Show that it is reasonable to approximate the cross section with a 1D polynomial in 
# logspace after ~178 nm.
# fig = plt.figure(figsize=(10, 8))
# plt.rcParams.update({'font.size': 16})
# plt.semilogy(a[:, 0], a[:, 2], color="cornflowerblue")
# plt.title("298K HDO cross section, data from Cheng+ 1999", y=1.02)
# plt.ylabel("Cross section (cm^2)")
# plt.xlabel("Wavelength (nm)")
# plt.show()

# fit function from 178 nm to 195 ----------------------------------------------
wvs_298 = data298[:, 0]
xsects_298 = data298[:, 2]

# slice the data to just the linear-in-logspace part
i = np.where(wvs_298 == 178)[0][0]
wvs_forfit = wvs_298[i:]
xsects_forfit = xsects_298[i:]
logxsects = np.log10(xsects_forfit)

# do the fit
z = np.polyfit(wvs_forfit, logxsects, 1)
f = np.poly1d(z)

# create new arrays for extrapolated data
newwv_298 = np.arange(195.2, 219.6, 0.2)
new_xsects_298 = 10**f(newwv_298)

# plot
# fig = plt.figure(figsize=(10, 8))
# plt.semilogy(data298[:, 0], data298[:, 2], "r-", label="Cheng+1999 data")
# plt.semilogy(newwv_298, new_xsects_298, "b-", label="extrapolated")
# plt.title("298K HDO cross section", y=1.02)
# plt.ylabel("Cross section (cm^2)")
# plt.xlabel("Wavelength (nm)")
# plt.legend()
# plt.show()
extrap_red_edge = np.column_stack((newwv_298, new_xsects_298))

#np.savetxt("HDO_extrapolated_298K.dat", extrap_red_edge, fmt=['%.1f', '%.7e'] )

# Use Cheng+ 2004 data to get cross sections between 125 to 140 nm -------------
cheng2004 = np.loadtxt("uvxsect/HDO_cheng2004.dat")

# We need the values at the half-wavelengths, so we interpolate.
newx = np.arange(125.5, 140, 1)
yinterp = np.interp(newx, cheng2004[:, 0], cheng2004[:, 1])
interpdata = np.column_stack((newx, yinterp))

# Plot the interpolated data
# fig = plt.figure(figsize=(10,8))
# plt.semilogy(cheng2004[:,0], cheng2004[:,1], color="cornflowerblue", label="Cheng+2004 data")
# plt.semilogy(newx, yinterp, 'rx', label="interpolated values")
# plt.title("298K HDO cross section, with interpolation of Cheng+2004", y=1.02)
# plt.ylabel("Cross section (cm^2)")
# plt.xlabel("Wavelength (nm)")
# plt.show()

#np.savetxt("HDOinterp.dat", interpdata, fmt=['%.1f', '%.7e'] )

# Now extrapolate on the blue edge to the wavelength we need, 121.5 nm. --------

# fit function for λ = 125, 126 because those are the only reasonable points....
wvs_fit = cheng2004[:, 0][0:2]
xsects_fit = cheng2004[:, 1][0:2]
logxsects = np.log10(xsects_fit)

# do the fit
z = np.polyfit(wvs_fit, logxsects, 1)
f = np.poly1d(z)

# create new arrays for extrapolated data
newwv_blu_298 = np.arange(121.5, 125.5, 1)
new_xsects_blu_298 = 10**f(newwv_blu_298)

# plot
# fig = plt.figure(figsize=(10, 8))
# plt.semilogy(cheng2004[:, 0], cheng2004[:, 1], "r-", label="Cheng+2004 data")
# plt.semilogy(newwv_blu_298, new_xsects_blu_298, "b-", label="extrapolated")
# plt.title("298K HDO cross section blue edge", y=1.02)
# plt.ylabel("Cross section (cm^2)")
# plt.xlabel("Wavelength (nm)")
# plt.legend()
# plt.show()
# extrapdata = np.column_stack((newwv_blu_298, new_xsects_blu_298))

#np.savetxt("HDO_blueedge_extrap_298K.dat", extrapdata, fmt=['%.1f', '%.7e'] )

# Now interpolate on the extrapolated red edge data ----------------------------

# We need the values at the half-wavelengths, so we interpolate.
newx = np.arange(195.5, 220.5, 1)
yinterp = np.interp(newx, extrap_red_edge[:, 0], extrap_red_edge[:, 1])
interpdata = np.column_stack((newx, yinterp))

# fig = plt.figure(figsize=(10, 8))
# plt.semilogy(extrap_red_edge[:, 0], extrap_red_edge[:, 1], 
#              color="cornflowerblue", label="extrapolated data")
# plt.semilogy(newx, yinterp, 'rx', label="interpolated values")
# plt.title("298K HDO red edge cross section, with interpolation", y=1.02)
# plt.ylabel("Cross section (cm^2)")
# plt.xlabel("Wavelength (nm)")
# plt.show()

# np.savetxt("HDO_rededge_interp.dat", interpdata, fmt=['%.1f', '%.7e'])

# Make the final figure ========================================================

# First put all the wavelengths and cross sections together from data, 
# extrapolations, and so forth for the 250K data.
all_250_wvs = np.hstack((newwv_blu_298, cheng2004[:, 0], data250[:, 0], newwv_250))
all_250_xsects = np.hstack((new_xsects_blu_298*0.5, cheng2004[:, 1]*0.5, data250[:, 1], new_xsects_250))

# interpolate all the 250K data
x_interp_250 = np.arange(121.5, 220, 1)  # wavelengths binned in 1 nm centered on whole integer wavelengths
y_interp_250 = np.interp(x_interp_250, all_250_wvs, all_250_xsects)
interpdata_250 = np.column_stack((x_interp_250, y_interp_250))

# Final figure of HDO cross sections at 250 and 298K
fig = plt.figure(figsize=(7, 5))
ax = plt.gca()
ax.set_facecolor("#ededed")
ax.grid(zorder=0, color="white", which="major")
for side in ["top", "bottom", "left", "right"]:
    ax.spines[side].set_visible(False)

c1 = "cornflowerblue"
plt.semilogy(newwv_blu_298, new_xsects_blu_298, color=c1, linestyle=":")#, label="extrapolated")
plt.semilogy(cheng2004[:,0], cheng2004[:,1], color=c1, linestyle="--", label="Cheng+2004 data (300K)")
plt.semilogy(data298[:,0], data298[:,2], color=c1, linestyle="-", label="Cheng+1999 data (298K)")
plt.semilogy(extrap_red_edge[:,0], extrap_red_edge[:,1], color=c1, linestyle=":")#, label="extrapolated")
# c2 = "xkcd:bright orange"
# plt.semilogy(np.hstack((newwv_blu_298, cheng2004[:,0])), np.hstack((new_xsects_blu_298*0.5,cheng2004[:,1]*0.5)), color=c2, linestyle="-.")
# plt.semilogy(newwv_250, new_xsects_250, color=c2, linestyle=":")#, label="extrapolated")
# plt.semilogy(data250[:,0], data250[:,2], color=c2, label="Cheng+1999 (250K)")
plt.title("HDO photodissociation cross section", y=1.02)
plt.ylabel("Cross section (cm^2)")
plt.xlabel("Wavelength (nm)")
plt.xlim(right=200)
plt.ylim(bottom=1e-23)

custom_lines = [Line2D([0], [0], color=c1, linestyle="-"),
                # Line2D([0], [0], color=c2, linestyle="-"), 
                Line2D([0], [0], color=c1, linestyle="--"), 
                Line2D([0], [0], color=c1, linestyle=":")]
                # Line2D([0], [0], color=c2, linestyle="-.")]

plt.legend(custom_lines, ["Cheng+1999 (298K)",
                          # "Cheng+1999 (250K)",
                          "Cheng+2004 (300K)", "extrapolated"], loc="lower left")
                          # "0.5*(Cheng+2004 data and extrapolation)"], loc="lower left")
plt.savefig(results_dir+"AllResultPlots/HDO_xsects.png", bbox_inches="tight", dpi=300)

np.savetxt("uvxsect/HDO_xsects.dat", interpdata, fmt=['%.1f', '%.7e'] )
