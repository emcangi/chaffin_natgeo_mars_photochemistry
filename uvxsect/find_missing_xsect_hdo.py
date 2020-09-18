"""
PYTHON code to interpolate the cross-section for HDO at wavelength = 152.6 nm, temp = 298K.
Data is from BM Cheng et al, 1999, DOI 10.1029/1999GL008367, but for some reason there is no data
recorded for the cross-section at 298K and 152.6 nm. 
Data from paper is found at http://140.110.206.16/wp-content/uploads/2014/02/H2O-HOD-and-D2O.pdf.
If this source is ever lost, try:
https://web.archive.org/web/20010303004339/http://ams-bmc.srrc.gov.tw:80/database.htm
What follows is a basic interpolation.

Author: Eryn Cangi
Date: November 2017
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp

# load the data, get rid of the 250K column 
a = np.loadtxt("HDO.dat")
b = np.delete(a, 1, 1)

# get rid of the bad line for fitting purposes
badline = np.where(b[:,0]==152.6)[0][0]
b = np.delete(b, badline, 0)

# check it
plt.plot(b[:,0], b[:,1])
plt.show()

# do the interpolation
f = interp.interp1d(b[:,0], b[:,1])

# plot
fig = plt.figure(figsize=(8,6))
ax = plt.gca()
ax.axvline(152.6) # missing datapoint was wavelength = 152.6 nm
xnew = np.loadtxt("HDO_wavelengths_corrected.txt")  # list of wavelengths with 152.6 included.
ynew = f(xnew)
plt.plot(b[:,0], b[:,1], 'o', xnew, ynew, '-')
plt.show()

# produce the missing value, which by the way is 1.855e-18.
i = list(xnew).index(152.6)
print(i)
print(list(ynew)[i])
