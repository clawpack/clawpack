#!/usr/bin/env python

r"""Bowl-Slosh regression test for GeoClaw

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


class BowlSloshTest(test.GeoClawRegressionTest):

    r"""Bowl-Slosh regression test for GeoClaw"""

    def setUp(self):

        super(BowlSloshTest, self).setUp()

        # Make topography
        a = 1.
        h0 = 0.1
        topo_func = lambda x,y: h0 * (x**2 + y**2) / a**2 - h0

        topo = topotools.Topography(topo_func=topo_func)
        topo.topo_type = 2
        topo.x = numpy.linspace(-2.0, 2.0, 200)
        topo.y = numpy.linspace(-2.0, 2.0, 200)
        topo.write(os.path.join(self.temp_path, "bowl.topotype2"), \
                topo_type=2, Z_format="%22.15e")

        from . import make_fgmax_grid 
        make_fgmax_grid.make_fgmax_grid1(self.temp_path)


    def runTest(self, save=False, indices=(2, 3)):
        r"""Test bowl-slosh example

        Note that this stub really only runs the code and performs no tests.

        """

        # Write out data files
        self.load_rundata()
        self.write_rundata_objects()

        # Run code
        self.run_code()

        # Perform tests
        self.check_gauges(save=save, gauge_id=1, indices=(2, 3))
        self.check_fgmax(save=save)
        self.success = True


if __name__=="__main__":
    if len(sys.argv) > 1:
        if bool(sys.argv[1]):
            # Fake the setup and save out output
            test = BowlSloshTest()
            try:
                test.setUp()
                test.runTest(save=True)
            finally:
                test.tearDown()
            sys.exit(0)
    unittest.main()
