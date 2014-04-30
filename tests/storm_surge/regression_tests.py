#!/usr/bin/env python

"""Regression test for GeoClaw's storm surge functionality"""

import sys
import os
import unittest

import numpy

import clawpack.geoclaw.tests as tests

class IkeTest(tests.GeoClawTest):

    r"""Hurricane Ike regression test"""

    def setUp(self):

        super(IkeTest, self).setUp()

        # Download topography


    def runTest(self, save=False):

        # Run code
        super(IkeTest, self).runTest()

        # Compare gauge data
        data = numpy.loadtxt(os.path.join(self.temp_path), 'fort.gauge')

        regression_data_file = os.path.join(self.test_path, "regression_data.txt")
        if save:
            numpy.savetxt(regression_data_file, data)
        regression_data = numpy.loadtxt(regression_data_file)

        # Compare data
        assert True, "\n data: %s, \n expected: %s" % (data, regression_data)

        # If we have gotten here then we do not need to copy the run results
        self.success = True



if __name__=="__main__":
    if len(sys.argv) > 1:
        if bool(sys.argv[1]):
            # Fake the setup and save out output
            test = IkeTest()
            try:
                test.setUp()
                test.runTest(save=True)
            finally:
                test.tearDown()
            sys.exit(0)
    unittest.main()