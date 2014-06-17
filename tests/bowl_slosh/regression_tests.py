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


    def runTest(self, save=False):

        # Run code
        super(BowlSloshTest, self).runTest()

        # Get gauge data
        data = numpy.loadtxt(os.path.join(self.temp_path, 'fort.gauge'))
        data_sum = [data[:,2].sum(), data[:,3].sum()]

        # Get (and save) regression comparison data
        regression_data_file = os.path.join(self.test_path, "regression_data.txt")
        if save:
            numpy.savetxt(regression_data_file, data)
        regression_data = numpy.loadtxt(regression_data_file)
        regression_sum = [regression_data[:,2].sum(), regression_data[:,3].sum()]
        # regression_sum = regression_data

        # Compare data
        tolerance = 1e-14
        assert numpy.allclose(data_sum, regression_sum, tolerance), \
                "\n data: %s, \n expected: %s" % (data_sum, regression_sum)
        assert numpy.allclose(data, regression_data, tolerance), \
                "Full gauge match failed."

        # If we have gotten here then we do not need to copy the run results
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