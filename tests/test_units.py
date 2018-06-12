#!/usr/bin/env python
# encoding: utf-8

from __future__ import absolute_import
from __future__ import print_function
import six
from six.moves import range

import numpy

from clawpack.geoclaw.units import units, convert


def test_conversions(verbose=False):
    r"""Test unit conversions."""

    for (measurement_type, measurement_units) in units.items():
        value = numpy.pi
        units_list = list(units[measurement_type].keys())
        for i in range(len(units_list)):
            if verbose:
                print("%s (%s) -> (%s)" % (value, units_list[i - 1], 
                                                  units_list[i]))
            value = convert(value, units_list[i - 1], units_list[i])
        numpy.testing.assert_allclose([value], [numpy.pi],
                      err_msg="Measurement tyep %s failed." % measurement_type)

if __name__ == '__main__':
    test_conversions(verbose=True)
