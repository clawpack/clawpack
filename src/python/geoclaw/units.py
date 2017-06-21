#!/usr/bin/env python
# encoding: utf-8
r"""
Defines common units used throughout GeoClaw along with functions for
converting between them.
"""

from __future__ import print_function
from __future__ import absolute_import

import six
from six.moves import range

import sys

import numpy

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
conversion_factor['miles'] = [None, 'miles']
conversion_factor['nm'] = [1852.0, 'nautical miles']
conversion_factor['lat-long'] = [LAT2METER, 'longitude-latitude']

# Pressure - Rigidity
unit_conversion_factor['Pa'] = [1.0, "pascals"]
unit_conversion_factor['KPa'] = [1e3, "kilopascals"]
unit_conversion_factor['MPa'] = [1e6, "megapascals"]
unit_conversion_factor['GPa'] = [1.e9, "gigapascals"]
unit_conversion_factor['mbar'] = [1e2, "millibar"]
unit_conversion_factor['dyne/cm^2'] = [0.1, "Dynes/cm^2"]
unit_conversion_factor['dyne/m^2'] = [1.e-5, "Dynes/m^2"]

# Speeds
unit_conversion_factor['m/s'] = [1.0, 'meters/second']
unit_conversion_factor['nm/h'] = [0.51444444, 'knots']

# Moments
unit_conversion_factor['N-m'] = [1.0, "Newton-Meters"]
unit_conversion_factor['dyne-cm'] = [1.e-7, "Dynes - Centimeter"]


# Unit conversion function
def convert_units(value, old_units, new_units, verbose=False):
    r""""""

    if old_units not in conversion_factor:
        raise ValueError("Units %s not found in list of supported ",
                         "conversions." % str(old_units))

    converted_value = value * conversion_factor[old_units] / conversion_factor[new_units]

    raise NotImplementedError("This has not yet been implemented.")
    return converted_value


def test_conversions():

    raise NotImplementedError("Test is not yet implemented.")

    # Check that these are consistent:
    check = [unit_conversion_factor[standard_units[param]] == 1. for param in
             standard_units.keys()]
    if not numpy.alltrue(check):
        raise ValueError("Conversion factors should be 1 for all standard_units")

if __name__ == '__main__':
    # Add commandline unit conversion capability
    pass
