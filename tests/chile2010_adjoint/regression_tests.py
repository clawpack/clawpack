#!/usr/bin/env python

r"""chile2010_adjoint regression test for GeoClaw

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

thisfile = os.path.realpath(__file__)
testdir = os.path.split(thisfile)[0]


class Chile2010AdjointTest(test.GeoClawRegressionTest):

    r"""Chile2010AdjointTest regression test for GeoClaw"""

    def setUp(self):

        super(Chile2010AdjointTest, self).setUp()

        # Make topography
        os.system('make -s topo')


    def runTest(self, save=False, indices=(2, 3)):
        r"""Test chile2010_adjoint example

        Note that this stub really only runs the code and performs no tests.

        """


        # Run adjoint problem
        adjointdir = testdir + '/adjoint'

        # Running the adjoint problem
        os.chdir(adjointdir)
        os.system('make -s topo') # also make qinit for adjoint
        os.system('make -s new')
        os.system('make .output > /dev/null')
        os.chdir(testdir)


        # Write out data files
        self.load_rundata()
        self.write_rundata_objects()

        # Run code
        self.run_code()

        # Perform tests
        self.check_gauges(save=save, gauge_id=1, indices=(2, 3))
        self.success = True


if __name__=="__main__":
    if len(sys.argv) > 1:
        if bool(sys.argv[1]):
            # Fake the setup and save out output
            test = Chile2010AdjointTest()
            try:
                test.setUp()
                test.runTest(save=True)
            finally:
                test.tearDown()
            sys.exit(0)
    unittest.main()
