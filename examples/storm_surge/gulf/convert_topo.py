#!/usr/bin/env python
# encoding: utf-8

import numpy as np
import scipy.interpolate as interp
import matplotlib.pyplot as plt

plot = False

raw_bathy_file = 'bathy/SWG_surveys.xyz'

data = np.loadtxt(raw_bathy_file,skiprows=1)

x = data[:,0]
y = data[:,1]
z = data[:,2]

# Create interpolation function
bathy = interp.interp2d(x[0:100],y[0:100],z[0:100])

# Extract info about data
lower_right = (np.min(x),np.min(y))
upper_left = (np.max(x),np.max(y))
print "Extents of bathymetry data:"
print "   Lower Right Extent = (%s,%s)" % (lower_right[0],lower_right[1])
print "   Upper Left Extent =  (%s,%s)" % (upper_left[0],upper_left[1])

# Plot
N = 100
X,Y = np.meshgrid(np.linspace(lower_right[0],upper_left[0],N),
                  np.linspace(lower_right[1],upper_left[1],N))
fig = plt.figure(1)
axes = fig.add_subplot(111)
axes.pcolor(X,Y,bathy(X,Y))
fig.add_colorbar()

plt.show()
