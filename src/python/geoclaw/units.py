#!/usr/bin/env python
# encoding: utf-8
r"""
Defines common units used throughout GeoClaw along with functions for
converting between them.
"""

from __future__ import print_function
from __future__ import absolute_import


import sys

from clawpack.geoclaw.data import LAT2METER

# Define baseline units so that we can convert to these and then back to
# requested.
standard_units = {'time': 's',
                  'spherical_length': 'lat-long',
                  'length': 'm',
                  'speed': 'm/s',
                  'radius': 'm',
                  'pressure': 'Pa'}

# Unit conversion definitions - handles conversions to standard units above
# The dictionary contains a list of values, the first is the conversion factor
# and the second a human readable version of the unit.
conversion_factor = {}

# Length
conversion_factor['m'] = [1.0, 'meters']
conversion_factor['cm'] = [1e-2, 'centimeters']
conversion_factor['km'] = [1e3, 'kilometers']
conversion_factor['miles'] = [1.60934e3, 'miles']
conversion_factor['nm'] = [1852.0, 'nautical miles']
conversion_factor['lat-long'] = [LAT2METER, 'longitude-latitude']
length_units = ['m', 'cm', 'km', 'miles', 'nm', 'lat-long']

# Pressure - Rigidity
conversion_factor['Pa'] = [1.0, "pascals"]
conversion_factor['KPa'] = [1e3, "kilopascals"]
conversion_factor['MPa'] = [1e6, "megapascals"]
conversion_factor['GPa'] = [1.e9, "gigapascals"]
conversion_factor['mbar'] = [1e2, "millibar"]
conversion_factor['dyne/cm^2'] = [0.1, "Dynes/cm^2"]
conversion_factor['dyne/m^2'] = [1.e-5, "Dynes/m^2"]
pressure_units = ['Pa', 'KPa', 'MPa', 'GPa', 'mbar', 'dyne/cm^2', 'dyne/m^2']

# Speeds
conversion_factor['m/s'] = [1.0, 'meters/second']
conversion_factor['knots'] = [0.51444444, 'knots (nm / hour)']
speed_units = ['m/s', 'knots']

# Moments
conversion_factor['N-m'] = [1.0, "Newton-Meters"]
conversion_factor['dyne-cm'] = [1.e-7, "Dynes - Centimeter"]
moment_units = ['N-m', 'dyne-cm']


# Unit conversion function
def convert_units(value, old_units, new_units, verbose=False):
    r""""""

    if old_units not in conversion_factor:
        raise ValueError("Units %s not found in list of supported ",
                         "conversions." % str(old_units))

    if new_units not in conversion_factor:
        raise ValueError("Units %s not found in list of supported ",
                         "conversions." % str(new_units))

    return value * conversion_factor[old_units][0] /    \
                   conversion_factor[new_units][0]


if __name__ == '__main__':
    # Add commandline unit conversion capability
    pass
