#!/usr/bin/env python
# encoding: utf-8

r"""Calculate refinement resolutions given ratios provided"""

import sys

import numpy as np

# Generic, spheriod based conversion
R_earth = 6378.1 * 1000.0
deg2meters = lambda theta,lat:R_earth * theta * np.pi / 180.0 * np.cos(lat * np.pi / 180.0)
meters2deg = lambda d,lat:d / (R_earth * np.pi / 180.0 * np.cos(lat * np.pi / 180.0))

# Based at lat = 24º
long2meters = lambda degree_resolution:degree_resolution * 100950.05720513177
lat2meters = lambda degree_resolution:degree_resolution * 110772.87259559495

def calculate_resolution(ratios,base_resolutions=[0.25,0.25]):
    r"""Given *ratios* and starting resolutions, calculate level resolutions"""
    num_levels = len(ratios) + 1

    degree_resolutions = np.empty((num_levels,2))
    meter_resolutions = np.empty((num_levels,2))
    degree_resolutions[0,:] = base_resolutions
    meter_resolutions[0,0] = long2meters(base_resolutions[0])
    meter_resolutions[0,1] = lat2meters(base_resolutions[1])
    for level in xrange(1,num_levels):
        degree_resolutions[level,:] = degree_resolutions[level-1,:] / ratios[level-1]
        meter_resolutions[level,0] = long2meters(degree_resolutions[level,0])
        meter_resolutions[level,1] = lat2meters(degree_resolutions[level,1])

    print "Resolutions:"
    for level in xrange(num_levels):
        print " Level %s - (%sº,%sº) - (%s m, %s m)" % (str(level+1),
                                                        degree_resolutions[level,0],
                                                        degree_resolutions[level,1],
                                                        meter_resolutions[level,0],
                                                        meter_resolutions[level,1])
    

if __name__ == "__main__":
    if len(sys.argv) > 1:
        ratios = [int(ratio) for ratio in sys.argv[1:]]
    else:
        ratios = [1]
    calculate_resolution(ratios)