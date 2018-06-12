#!/usr/bin/env python

r"""Multilayer Shallow Water Test Case

To create new regression data use
    `python regression_tests.py True`
"""

from __future__ import absolute_import
import os
import sys
import unittest

import numpy

import clawpack.geoclaw.test as test
import clawpack.geoclaw.topotools as topotools


def transform_p2c(x, y, x0, y0, theta):
    return ( x*numpy.cos(theta) + y*numpy.sin(theta) - x0,
            -x*numpy.sin(theta) + y*numpy.cos(theta) - y0)

# Bathymetry
def bathy_step(x, y, location=0.15, angle=0.0, left=-1.0, right=-0.2):
    x_c,y_c = transform_p2c(x, y, location, 0.0, angle)
    return ((x_c <= 0.0) * left
          + (x_c >  0.0) * right)


class MultilayerTest(test.GeoClawRegressionTest):

    r"""Multilayer plane-wave regression test for GeoClaw

    initial condition angle = numpy.pi / 4.0
    bathy_angle = numpy.pi / 8.0

    """

    def setUp(self):

        super(MultilayerTest, self).setUp()

        # Make topography
        topo_func = lambda x, y: bathy_step(x, y, location=0.15,
                                                  angle=numpy.pi / 8.0,
                                                  left=-1.0, right=-0.2)
        topo = topotools.Topography(topo_func=topo_func)
        topo.x = numpy.linspace(-1.16, 2.16, 166)
        topo.y = numpy.linspace(-1.16, 2.16, 166)
        topo.write(os.path.join(self.temp_path, "jump_topo.topotype2"))

    def runTest(self, save=False):
        r"""Test multi-layer basic plane-waves."""

        # Load and write data, change init-condition's starting angle
        self.load_rundata()
        self.rundata.qinit_data.angle = numpy.pi / 4.0
        self.write_rundata_objects()

        # Run code and check surface heights
        self.run_code()
        self.check_gauges(save=save, gauge_id=0, indices=(6, 7), atol=1e-5)
        self.check_gauges(save=save, gauge_id=1, indices=(6, 7), atol=1e-5)
        self.check_gauges(save=save, gauge_id=2, indices=(6, 7), atol=1e-5)
        self.check_gauges(save=save, gauge_id=3, indices=(6, 7), atol=1e-5)
        self.check_gauges(save=save, gauge_id=4, indices=(6, 7), atol=1e-5)

        # If we have gotten here then we do not need to copy the run results
        self.success = True


if __name__ == "__main__":
    if len(sys.argv) > 1:
        if bool(sys.argv[1]):
            # Fake the setup and save out output
            test = MultilayerTest()
            try:
                test.setUp()
                test.runTest(save=True)
            finally:
                test.tearDown()
            sys.exit(0)
    unittest.main()
