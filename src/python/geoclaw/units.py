#!/usr/bin/env python
# encoding: utf-8
r"""
Defines common units used throughout GeoClaw along with functions for
converting between them.
"""

from __future__ import print_function
from __future__ import absolute_import

import sys
import collections

from clawpack.geoclaw.data import LAT2METER

# Define baseline units so that we can convert to these and then back to
# requested.
standard_units = {'time': 's',
                  'spherical_length': 'lat-long',
                  'length': 'm',
                  'speed': 'm/s',
                  'radius': 'm',
                  'pressure': 'Pa',
                  'temperature': 'C'}

# Unit conversion definitions - handles conversions to standard units above
# The dictionary contains a list of values, the first is the conversion factor
# and the second a human readable version of the unit.
conversion_func = {}
units = {}

# Length
conversion_func['m'] = [lambda L: L,
                        lambda L: L]
conversion_func['cm'] = [lambda L: L * 1e-2,
                         lambda L: L / 1e-2]
conversion_func['km'] = [lambda L: L * 1e3,
                         lambda L: L / 1e3]
conversion_func['miles'] = [lambda L: L * 1.60934e3,
                            lambda L: L / 1.60934e3]
conversion_func['nmi'] = [lambda L: L * 1852.0,
                          lambda L: L / 1852.0]
conversion_func['lat-long'] = [lambda L: L * LAT2METER,
                               lambda L: L / LAT2METER]
units['length'] = collections.OrderedDict({'m': 'meters', 'cm': 'centimeters', 
                                           'km': 'kilometers', 'miles': 'miles',
                                           'nmi': 'nautical miles', 
                                           'lat-long': 'longitude-latitude'})

# Pressure - Rigidity
conversion_func['Pa'] = [lambda P: P, lambda P: P]
conversion_func['hPa'] = [lambda P: P * 1e2, 
                          lambda P: P / 1e2]
conversion_func['KPa'] = [lambda P: P * 1e3,
                          lambda P: P / 1e3]
conversion_func['MPa'] = [lambda P: P * 1e6,
                          lambda P: P / 1e6]
conversion_func['GPa'] = [lambda P: P * 1.e9,
                          lambda P: P / 1.e9]
conversion_func['mbar'] = [lambda P: P * 1e2,
                           lambda P: P / 1e2]
conversion_func['dyne/cm^2'] = [lambda P: P * 0.1,
                                lambda P: P / 0.1]
conversion_func['dyne/m^2'] = [lambda P: P * 1.e-5,
                               lambda P: P / 1.e-5]
units['pressure'] = {'Pa': 'pascals', 'hPa': 'hectopascals', 
                     'KPa': 'kilopascals', 'MPa': 'megapascals', 
                     'GPa': 'gigapascals', 'mbar': 'millibar', 
                     'dyne/cm^2': 'Dynes/cm^2', 'dyne/m^2': 'Dynes/m^2'}

# Speeds
conversion_func['m/s'] = [lambda v: v, lambda v: v]
conversion_func['knots'] = [lambda v: v * 0.51444444, 
                            lambda v: v / 0.51444444]
units['speed'] = {'m/s': 'meters/second', 'knots': 'knots (nm / hour)'}

# Moments
conversion_func['N-m'] = [lambda M: M, lambda M: M]
conversion_func['dyne-cm'] = [lambda M: M * 1.e-7,
                              lambda M: M / 1.e-7]
units['moment'] = {'N-m': "Newton-Meters", 'dyne-cm': "Dynes-Centimeter"}

# Temperature
conversion_func['C'] = [lambda temp: temp, lambda temp: temp]
conversion_func['F'] = [lambda temp: (temp - 32.0) * 5.0 / 9.0, 
                        lambda temp: temp * 9.0 / 5.0 + 32.0]
conversion_func['K'] = [lambda temp: temp + 273.15,
                        lambda temp: temp - 273.15]
units['temperature'] = {'C': "Celsius", 'F': "Fahrenheit", 'K': "Kelvin"}


def units_available():
    r"""
    Constructs a string suitable for reading detailing the units available.
    """
    output = ""
    for (measurement_type, measurement_units) in units.items():
        output = "\n".join((output, measurement_type.capitalize()))
        for (abbrv, full_name) in measurement_units.items():
            output = "\n".join((output, "  %s (%s)" % (abbrv, full_name)))

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

    found_type = None
    for (measurement_type, measurement_units) in units.items():
        if old_units in measurement_units:
            found_type = measurement_type
            break

    if found_type is None:
        raise ValueError("Units %s not found in list of " % str(old_units),
                         "supported conversions." )

    if new_units not in units[found_type].keys():
        raise ValueError("Units %s not found in list of " % str(new_units),
                         "supported conversions of %s type." % found_type)

    if verbose:
        print("Convert %s %s to %s." % (value, old_units, new_units))

    return conversion_func[new_units][1](conversion_func[old_units][0](value))


if __name__ == '__main__':
    # Add commandline unit conversion capability
    if len(sys.argv) == 4:
        convert(float(sys.argv[1]), sys.argv[2], sys.argv[3])
    else:
        # Usage and available units
        print("Usage:  Convert value in units to new units.")
        print("  units <value> <old units> <new units>")
        print("where <old units> and <new units> are one of the available")
        print("units listed below.")
        print("")
        print("Available Units:")
        print("  First value is the abbreviation that should be used as an")
        print("  input unit while the second is the full name of the unit.")
        print(units_available())
