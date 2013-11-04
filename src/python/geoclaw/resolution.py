#!/usr/bin/env python
# encoding: utf-8

r"""Calculate refinement resolutions given ratios provided"""

import sys
import argparse

import numpy as np

# Generic, spheriod based conversion
R_earth = 6378.1 * 1000.0
deg2meters = lambda theta,lat:R_earth * theta * np.pi / 180.0 * np.cos(lat * np.pi / 180.0)
meters2deg = lambda d,lat:d / (R_earth * np.pi / 180.0 * np.cos(lat * np.pi / 180.0))

# Based at lat = 24º
long2meters = lambda degree_resolution:degree_resolution * 100950.05720513177
lat2meters = lambda degree_resolution:degree_resolution * 110772.87259559495

def calculate_resolution(ratios, base_resolutions=[0.25,0.25], lat_long=True):
    r"""Given *ratios* and starting resolutions, calculate level resolutions"""
    num_levels = len(ratios) + 1

    # Calculate resolution on each level
    natural_resolutions = np.empty((num_levels, 2))
    natural_resolutions[0,:] = base_resolutions
    for level in xrange(1,num_levels):
        natural_resolutions[level, :] = natural_resolutions[level-1, :] / ratios[level-1]

    # Print out and convert to meters if applicable
    if lat_long:
        meter_resolutions = np.empty((num_levels,2))
        for level in xrange(num_levels):
            meter_resolutions[level,0] = long2meters(natural_resolutions[level,0])
            meter_resolutions[level,1] = lat2meters(natural_resolutions[level,1])


    print "Resolutions:"
    for level in xrange(num_levels):
        level_string = " Level %s" % str(level + 1)
        if lat_long:
            level_string = " - ".join((level_string, "(%sº,%sº)" % (natural_resolutions[level,0], natural_resolutions[level,1])))
            level_string = " - ".join((level_string, "(%s m,%s m)" % (meter_resolutions[level,0], meter_resolutions[level,1])))
        else:
            level_string = " - ".join((level_string, "(%s m,%s m)" % (natural_resolutions[level,0], natural_resolutions[level,1])))
        print level_string
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="resolutions", 
                   description="Compute the effective resolution of refinement")
    parser.add_argument('ratios', metavar="ratios", type=int, nargs="+",
                        help="Ratios used in refinement")
    parser.add_argument('--latlong','-l', dest='lat_long', action='store_true',
                        help="Computer assuming lat-long base grid")
    parser.add_argument('--base', metavar="resolution", dest='base_resolutions', 
                        action='store', nargs=2, default=[0.25,0.25], 
                        help="Base resolutions")
    args = parser.parse_args()
    
    calculate_resolution(args.ratios, base_resolutions=args.base_resolutions,
                         lat_long=args.lat_long)