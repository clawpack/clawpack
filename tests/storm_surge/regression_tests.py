#!/usr/bin/env python

"""Regression test for GeoClaw's storm surge functionality"""

import sys
import os
import unittest
import shutil

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


    def runTest(self, save=False, indices=(2, 3)):
        r"""Storm Surge Regression Test

        :Input:
         - *save* (bool) - If *True* will save the output from this test to 
           the file *regresion_data.txt*.  Passed to *check_gauges*.  Default is
           *False*.
         - *indices* (tuple) - Contains indices to compare in the gague 
           comparison and passed to *check_gauges*.  Defaults to *(2, 3)*.

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
            test = IkeTest()
            try:
                test.setUp()
                test.runTest(save=True)
            finally:
                test.tearDown()
            sys.exit(0)
    unittest.main()