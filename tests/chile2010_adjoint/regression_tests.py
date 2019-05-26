#!/usr/bin/env python

r"""chile2010_adjoint regression test for GeoClaw

To create new regression data use
    `python regression_tests.py True`
"""

from __future__ import absolute_import
import os
import sys
import unittest
import shutil

import numpy

import clawpack.geoclaw.test as test
import clawpack.geoclaw.topotools as topotools
from clawpack.clawutil.test import wip


try:
    CLAW = os.environ['CLAW']
except:
    raise Exception("*** Must first set CLAW enviornment variable")

# Scratch directory for storing topo and dtopo files:
scratch_dir = os.path.join(CLAW, 'geoclaw', 'scratch')

thisfile = os.path.realpath(__file__)
testdir = os.path.split(thisfile)[0]


class Chile2010AdjointTest(test.GeoClawRegressionTest):

    r"""Chile2010AdjointTest regression test for GeoClaw"""

    def setUp(self):

        super(Chile2010AdjointTest, self).setUp()

        start_dir = os.getcwd()
        test_adjoint_path = os.path.join(self.test_path, 'adjoint')
        temp_adjoint_path = os.path.join(self.temp_path, 'adjoint')
        #print('+++ test_adjoint_path = ',test_adjoint_path)
        #print('+++ temp_adjoint_path = ',temp_adjoint_path)

        shutil.copytree(test_adjoint_path, temp_adjoint_path)

        # run adjoint code
        os.chdir(temp_adjoint_path)
        #print('+++ Running adjoint in directory ',os.getcwd())
        #print('+++   contents: ', os.listdir('.'))
        os.system('make -s topo')
        os.system('make -s data')
        os.system('make -s new')
        os.system('make .output > output.txt')
        #print('+++   contents of _output: ', os.listdir('_output'))

        # set up forward code
        shutil.copy(os.path.join(self.test_path, "maketopo.py"),
                                 self.temp_path)
        os.chdir(self.temp_path)
        #print('+++ Running forward in directory ',os.getcwd())
        #print('+++   contents: ', os.listdir())
        os.system('python maketopo.py')
        #print('+++   scratch directory: ', scratch_dir)
        #print('+++   contents of scratch: ', os.listdir(scratch_dir))
        os.chdir(start_dir)


    def runTest(self, save=False, indices=(2, 3)):
        r"""Test chile2010_adjoint example

        Note that this stub really only runs the code and performs no tests.

        """

        # Write out data files
        self.load_rundata()
        temp_adjoint_path = os.path.join(self.temp_path, 'adjoint')
        self.rundata.adjointdata.adjoint_outdir = \
             os.path.join(temp_adjoint_path, '_output')
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
