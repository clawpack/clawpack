#!/usr/bin/env python

r"""Bowl-Slosh regression test for GeoClaw

To create new regression data use
    `python regression_tests.py True`
"""

import os
import sys
import unittest
import subprocess
import tempfile
import time

import numpy

import clawpack.geoclaw.test as test
import clawpack.geoclaw.topotools as topotools

class NetCDFBowlSloshTest(test.GeoClawRegressionTest):

    r"""NetCDF regression test for GeoClaw based on the bowl-slosh example"""

    def __init__(self, methodName="runTest"):

        super(NetCDFBowlSloshTest, self).__init__(methodName=methodName)

        self.netcdf_passed = False


    def setUp(self):

        self.temp_path = tempfile.mkdtemp()

        self.stdout = open(os.path.join(self.temp_path, "run_output.txt"), "w")
        self.stdout.write("Output from Test %s\n" % self.__class__.__name__)
        # TODO - Should change this to use the time module's formatting 
        # apparatus
        tm = time.localtime()
        year = str(tm[0]).zfill(4)
        month = str(tm[1]).zfill(2)
        day = str(tm[2]).zfill(2)
        hour = str(tm[3]).zfill(2)
        minute = str(tm[4]).zfill(2)
        second = str(tm[5]).zfill(2)
        date = 'Started %s/%s/%s-%s:%s.%s\n' % (year,month,day,hour,minute,second)
        self.stdout.write(date)
        self.stdout.write(("="*80 + "\n"))

        self.stderr = open(os.path.join(self.temp_path, "error_output.txt"), "w")
        self.stderr.write("Errors from Test %s\n" % self.__class__.__name__)
        self.stderr.write(date)
        self.stderr.write(("="*80 + "\n"))

        self.stdout.flush()
        self.stderr.flush()

        self.stdout.write("Paths:")
        self.stdout.write("  %s" % self.temp_path)
        self.stdout.write("  %s" % self.test_path)
        self.stdout.flush()

        # Make topography
        a = 1.
        h0 = 0.1
        topo_func = lambda x,y: h0 * (x**2 + y**2) / a**2 - h0

        topo = topotools.Topography(topo_func=topo_func)
        topo.x = numpy.linspace(-2.1, 2.1, 210)
        topo.y = numpy.linspace(-2.1, 2.1, 210)
        try:
            topo.write(os.path.join(self.temp_path, "bowl.nc"))
        
        except ImportError:
            # Assume that NetCDF is not installed and move on
            self.netcdf_passed = False
            self.success = True
            self.stdout.write("NetCDF topography test skipped due to failure" + 
                              "to build test program.")

        else:
            self.build_executable()

    def build_executable(self):
        try:
            self.stdout.write("Teting NetCDF output:\n")
            subprocess.check_call("cd %s ; make netcdf_test " % self.test_path, 
                                                                stdout=self.stdout,
                                                                stderr=self.stderr,
                                                                shell=True)
        except subprocess.CalledProcessError:

            self.stdout.write("NetCDF topography test skipped due to failure" + 
                              "to build test program.")
            self.success = True
            self.netcdf_passed = False

        else:
            # Force recompilation of topo_module to add NetCDF flags
            mod_path = os.path.join(os.environ["CLAW"], "geoclaw", "src", "2d",
                                    "shallow", "topo_module.mod")
            obj_path = os.path.join(os.environ["CLAW"], "geoclaw", "src", "2d",
                                    "shallow", "topo_module.o")
            if os.path.exists(mod_path):
                os.remove(mod_path)
            if os.path.exists(obj_path):
                os.remove(obj_path)

            self.netcdf_passed = True
            super(NetCDFBowlSloshTest, self).build_executable()


    def runTest(self, save=False, indices=(2, 3)):
        r"""Test NetCDF topography support

        """

        # Check to see if NetCDF has been built
        if self.netcdf_passed:
            # Write out data files
            self.load_rundata()
            self.write_rundata_objects()

            # Run code
            self.run_code()

            # Perform tests
            self.check_gauges(save=save, indices=(2, 3), tolerance=1e-4)
            self.success = True


    def check_gauges(self, save=False, gauge_num=1, indices=[0], 
                           regression_data_path="regression_data.txt",
                           tolerance=1e-14):
        r"""Basic test to assert gauge equality

        :Input:
         - *save* (bool) - If *True* will save the output from this test to 
           the file *regresion_data.txt*.  Default is *False*.
         - *indices* (tuple) - Contains indices to compare in the gague 
           comparison.  Defaults to *(2, 3)*.
         - *regression_data_path* (path) - Path to the regression test data.
           Defaults to 'regression_data.txt'.
        """

        # Get gauge data
        data = numpy.loadtxt(os.path.join(self.temp_path, 'fort.gauge'))
        data_sum = []
        gauge_numbers = numpy.array(data[:, 0], dtype=int)
        gauge_indices = numpy.nonzero(gauge_num == gauge_numbers)[0]
        for index in indices:
            data_sum.append(data[gauge_indices, index].sum())

        # Get (and save) regression comparison data
        regression_data_file = os.path.join(self.test_path,
                                            regression_data_path)
        if save:
            numpy.savetxt(regression_data_file, data)
        regression_data = numpy.loadtxt(regression_data_file)
        regression_sum = []
        gauge_numbers = numpy.array(data[:, 0], dtype=int)
        gauge_indices = numpy.nonzero(gauge_num == gauge_numbers)[0]
        for index in indices:
            regression_sum.append(regression_data[gauge_indices, index].sum())
        # regression_sum = regression_data

        # Compare data
        assert numpy.allclose(data_sum, regression_sum, tolerance), \
                "\n data: %s, \n expected: %s" % (data_sum, regression_sum)
        assert numpy.allclose(data, regression_data, tolerance), \
                "Full gauge match failed."


    def tearDown(self):
        r"""Tear down test infrastructure.

        This version removes source that may have been compiled with NetCDF
        turned on so that it does not conflict with subsequent tests.

        """

        # Force recompilation of topo_module to add NetCDF flags
        mod_path = os.path.join(os.environ["CLAW"], "geoclaw", "src", "2d",
                                "shallow", "topo_module.mod")
        obj_path = os.path.join(os.environ["CLAW"], "geoclaw", "src", "2d",
                                "shallow", "topo_module.o")
        if os.path.exists(mod_path):
            os.remove(mod_path)
        if os.path.exists(obj_path):
            os.remove(obj_path)

        super(NetCDFBowlSloshTest, self).tearDown()



if __name__=="__main__":
    if len(sys.argv) > 1:
        if bool(sys.argv[1]):
            # Fake the setup and save out output
            test = NetCDFBowlSloshTest()
            try:
                test.setUp()
                test.runTest(save=True)
            finally:
                test.tearDown()
            sys.exit(0)
    unittest.main()
