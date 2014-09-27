#!/usr/bin/env python

"""Regression test for GeoClaw's storm surge functionality"""

import sys
import os
import unittest
import gzip

import numpy

import clawpack.geoclaw.tests as tests
import clawpack.geoclaw.topotools

class IkeTest(tests.GeoClawTest):

    r"""Hurricane Ike regression test"""

    def setUp(self):

        super(IkeTest, self).setUp()

        # Download storm data
        remote_url = "http://ftp.nhc.noaa.gov/atcf/archive/2008/bal092008.dat.gz"
        path = self.get_remote_file(remote_url, unpack=False)
        storm_path = os.path.join(os.path.dirname(path), 'ike.storm')

        # Need to additionally deal with the fact the file is gzipped
        with gzip.GzipFile(path, 'r') as gzip_file:
            file_content = gzip_file.read()
        
        with open(storm_path, 'w') as out_file:
            out_file.write(file_content)

        # Download file
        self.get_remote_file(
           "http://www.columbia.edu/~ktm2132/bathy/gulf_caribbean.tt3.tar.bz2")

        # Create synthetic bathymetry - needs more work
        topo = clawpack.geoclaw.topotools.Topography()
        topo.x = numpy.linspace(-100, -69, 124)
        topo.y = numpy.linspace(7.0, 33.0, 104)
        topo.Z = 25.0 * ((topo.X + 84.5)**2 + (topo.Y - 20.0)**2) - 4000.0
        topo.write(os.path.join(self.temp_path, 'gulf_caribbean.tt3'))


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
