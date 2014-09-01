#!/usr/bin/env python

r"""Bowl-Slosh regression test for GeoClaw

To create new regression data use
    `python regression_tests.py True`
"""

import os
import sys
import unittest

import numpy

import clawpack.geoclaw.tests as tests
import clawpack.geoclaw.topotools as topotools


class BowlSloshTest(tests.GeoClawTest):

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
        topo.write(os.path.join(self.temp_path, "bowl.topotype2"))


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
        self.check_gauges(save=save, indices=(2, 3))


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