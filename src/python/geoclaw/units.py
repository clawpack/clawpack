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


def units_available():
    r"""
    Constructs a string suitable for reading detailing the units available.
    """
    output = ""
    output = "\n".join((output, "Length"))
    for unit in length_units:
        output = "\n".join((output, "  %s - %s" % (conversion_factor[unit][1],
                                                   unit)))
    output = "\n".join((output, "Pressure - Rigidity"))
    for unit in pressure_units:
        output = "\n".join((output, "  %s - %s" % (conversion_factor[unit][1],
                                                   unit)))
    output = "\n".join((output, "Speeds"))
    for unit in speed_units:
        output = "\n".join((output, "  %s - %s" % (conversion_factor[unit][1],
                                                   unit)))
    output = "\n".join((output, "Moments"))
    for unit in moment_units:
        output = "\n".join((output, "  %s - %s" % (conversion_factor[unit][1],
                                                   unit)))
    return output


# Unit conversion function
def convert(value, old_units, new_units, verbose=False):
    r"""Convert *value* from *old_units* to *new_units*

    :Note:
    Currently this function only handles multiplicative conversions.  The
    reasoning behind not just returning this conversion factor as this function
    in the future will also support more complex unit conversions, e.g.
    converting between temperature scales.

    :Input:
     - *value* (ndarray or float) The value(s) to be converted.
     - *old_units* (string) Type of units that value comes in as.
     - *new_units* (string) Type of units that value should be converted to.
     - *verbose* (bool) Verbose output (default is False)

    :Output:
     - (ndarray or float) The converted value(s)
    """

    if old_units not in conversion_factor:
        raise ValueError("Units %s not found in list of supported ",
                         "conversions." % str(old_units))

    if new_units not in conversion_factor:
        raise ValueError("Units %s not found in list of supported ",
                         "conversions." % str(new_units))

    if verbose:
        print("Convert %s %s to %s." % (value, old_units, new_units))

    return value * conversion_factor[old_units][0] /    \
                   conversion_factor[new_units][0]


if __name__ == '__main__':
    # Add commandline unit conversion capability
    if len(sys.argv) == 4:
        convert(float(sys.argv[1], sys.argv[2], sys.argv[3]))
    else:
        # Usage and available units
        print("Usage:  Convert value in units to new units.")
        print("  units <value> <old units> <new units>")
        print("where <old units> and <new units> are one of the available")
        print("units listed below.")
        print("")
        print("Available Units:")
        print(units_available())
