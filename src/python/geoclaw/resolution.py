#!/usr/bin/env python
# encoding: utf-8

r"""Calculate refinement resolutions given ratios provided"""

from __future__ import absolute_import
from __future__ import print_function
import argparse

import numpy as np

from . import topotools
from six.moves import range

def calculate_resolution(ratios, base_resolutions=[0.25,0.25], 
                                 lat_long=True,
                                 latitude=24.0):
    r"""Given *ratios* and starting resolutions, calculate level resolutions"""
    num_levels = len(ratios) + 1

    # Calculate resolution on each level
    natural_resolutions = np.empty((num_levels, 2))
    natural_resolutions[0,:] = base_resolutions
    for level in range(1,num_levels):
        natural_resolutions[level, :] = natural_resolutions[level-1, :] / ratios[level-1]

    # Print out and convert to meters if applicable
    if lat_long:
        meter_resolutions = np.empty((num_levels,2))
        for level in range(num_levels):
            meter_resolutions[level,:] = topotools.dist_latlong2meters(
                                                   natural_resolutions[level,0],
                                                   natural_resolutions[level,1],
                                                   latitude)


    print("Resolutions:")
    for level in range(num_levels):
        level_string = " Level %s" % str(level + 1)
        if lat_long:
            level_string = " - ".join((level_string, "(%sº,%sº)" % (natural_resolutions[level,0], natural_resolutions[level,1])))
            level_string = " - ".join((level_string, "(%s m,%s m)" % (meter_resolutions[level,0], meter_resolutions[level,1])))
        else:
            level_string = " - ".join((level_string, "(%s m,%s m)" % (natural_resolutions[level,0], natural_resolutions[level,1])))
        print(level_string)
    

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
    parser.add_argument('--lat', dest='latitude', action='store', default=24.0,
                        help="Latitude to use in degrees to meters conversion")
    args = parser.parse_args()
    
    calculate_resolution(args.ratios, base_resolutions=args.base_resolutions,
                         lat_long=args.lat_long, latitude=args.latitude)