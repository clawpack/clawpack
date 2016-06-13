#!/usr/bin/env python

r"""Regression tests.  Execute via:
    python regression_tests.py
to test, or
    python regression_tests.py True
to create new regression data for archiving.
"""

import os
import sys
import unittest
import shutil

import numpy

import clawpack.geoclaw.test as test
import clawpack.geoclaw.topotools as topotools
import clawpack.geoclaw.dtopotools as dtopotools

class DTopoTests(test.GeoClawRegressionTest):
    
    def setUp(self):

        super(DTopoTests, self).setUp()

        # Make topography
        h0 = 1000.0
        topo_func = lambda x,y: -h0*(1 + 0.5 * numpy.cos(x - y))
        topo = topotools.Topography(topo_func=topo_func)
        topo.topo_type = 2
        topo.x = numpy.linspace(-10.0, 10.0, 201)
        topo.y = numpy.linspace(-10.0, 10.0, 201)
        topo.write(os.path.join(self.temp_path, "topo1.topotype2"), \
                topo_type=2, Z_format="%22.15e")

        h0 = 1000.0
        topo_func = lambda x,y: -h0*(1. + numpy.exp(x+y))
        topo = topotools.Topography(topo_func=topo_func)
        topo.topo_type = 2
        topo.x = numpy.linspace(-0.5, -0.3, 21)
        topo.y = numpy.linspace(-0.1, 0.4, 51)
        topo.write(os.path.join(self.temp_path, "topo2.topotype2"), \
                topo_type=2, Z_format="%22.15e")

        # Make dtopography
        subfault_path = os.path.join(self.test_path, "dtopo1.csv")
        input_units = {'slip': 'm', 'depth': 'km', 'length': 'km', 'width': 'km'}
        fault = dtopotools.CSVFault()
        fault.read(subfault_path, input_units=input_units, 
                coordinate_specification="top center")
        fault.rupture_type = 'dynamic'
        times = numpy.linspace(0.0, 1.0, 25)
        x = numpy.linspace(-0.4,0.6,151)
        y = numpy.linspace(-0.4,0.4,121)
        dtopo = fault.create_dtopography(x,y,times=times)
        dtopo.write(os.path.join(self.temp_path, "dtopo1.tt3"), dtopo_type=3)

        subfault_path = os.path.join(self.test_path, "dtopo2.csv")
        input_units = {'slip': 'm', 'depth': 'km', 'length': 'km', 'width': 'km'}
        fault = dtopotools.CSVFault()
        fault.read(subfault_path, input_units=input_units, 
                    coordinate_specification="top center")
        fault.rupture_type = 'dynamic'
        times = numpy.linspace(0.5, 1.2, 25)
        x = numpy.linspace(-0.9,0.1,201)
        y = numpy.linspace(-0.4,0.4,161)
        dtopo = fault.create_dtopography(x,y,times=times)
        dtopo.write(os.path.join(self.temp_path, "dtopo2.tt3"), dtopo_type=3)

        # copy existing file:
        shutil.copy(os.path.join(self.test_path, "dtopo3.tt1"),
                                 self.temp_path)


    def runTest(self, save=False, indices=(2, 3)):
        r"""DTopography basic regression test

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
            test = DTopoTests()
            try:
                test.setUp()
                test.runTest(save=True)
            finally:
                test.tearDown()
            sys.exit(0)
    unittest.main()

