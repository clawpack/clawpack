#!/usr/bin/env python
# encoding: utf-8

from __future__ import absolute_import
from __future__ import print_function
import six
from six.moves import range

import numpy

from clawpack.geoclaw.units import *


def test_conversions():
    r"""Test unit conversions."""

    value = numpy.pi
    for i in range(1, len(length_units)):
        value = convert(value, length_units[i - 1], length_units[i])
    value = convert(value, length_units[-1], length_units[0])
    numpy.testing.assert_allclose([value], [numpy.pi],
                                  err_msg="Length conversions failed.")

    value = numpy.pi
    for i in range(1, len(pressure_units)):
        value = convert(value, pressure_units[i - 1], pressure_units[i])
    value = convert(value, pressure_units[-1], pressure_units[0])
    numpy.testing.assert_allclose([value], [numpy.pi],
                                  err_msg="Pressure conversions failed.")

    value = numpy.pi
    for i in range(1, len(speed_units)):
        value = convert(value, speed_units[i - 1], speed_units[i])
    value = convert(value, speed_units[-1], speed_units[0])
    numpy.testing.assert_allclose([value], [numpy.pi],
                                  err_msg="Speed conversions failed.")

    value = numpy.pi
    for i in range(1, len(moment_units)):
        value = convert(value, moment_units[i - 1], moment_units[i])
    value = convert(value, moment_units[-1], moment_units[0])
    numpy.testing.assert_allclose([value], [numpy.pi],
                                  err_msg="Moment conversions failed.")

if __name__ == '__main__':
    test_conversions()
