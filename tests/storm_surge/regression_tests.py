#!/usr/bin/env python

"""Regression test for GeoClaw's storm surge functionality"""

import sys
import os
import unittest
import glob
import shutil

import numpy

import clawpack.geoclaw.tests as tests

class IkeTest(tests.GeoClawTest):

    r"""Hurricane Ike regression test"""

    def setUp(self):

        super(IkeTest, self).setUp()

        # Download topography
        self.get_remote_file("http://users.ices.utexas.edu/~kyle/bathy/" + \
                                     "gulf_caribbean.tt3.tar.bz2")
        # Download storm data
        # Eventually probably want to do this
        # self.get_remote_file("http://ftp.nhc.noaa.gov/atcf/archive/2008/" + \
        #                      "aal082008.dat.gz")
        #
        shutil.copy(os.path.join(self.test_path, 'ike.storm'), self.temp_path)


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