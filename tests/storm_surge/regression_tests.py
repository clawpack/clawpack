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

    def runTest(self, save=False):

        # Run code
        super(IkeTest, self).runTest()

        # Could compare AMR behavior...
        # num_levels = 2
        # file_list = glob.glob(os.path.join(self.temp_path,
                                                # "_output", 
                                                # "fort.q*"))

        # time = numpy.empty(len(file_list), dtype=float)
        # num_grids = numpy.zeros((time.shape[0], num_levels), dtype=int)
        # num_cells = numpy.zeros((time.shape[0], num_levels), dtype=int)

        # for (n,path) in enumerate(file_list):
        #     # Read t file
        #     t_path = path[:-5] + "t" + path[-4:]
        #     with open(t_path, 'r') as t_file:
        #         time[n] = seconds2days(float(t_file.readline().split()[0]) - landfall)
        #         t_file.readline()
        #         t_file_num_grids = int(t_file.readline().split()[0])

        #     # Read q_file
        #     with open(path, 'r') as q_file:
        #         line = "\n"
        #         while line != "":
        #             line = q_file.readline()
        #             if "grid_number" in line:
        #                 # print "grid number:", int(line.split()[0])
        #                 level = int(q_file.readline().split()[0])
        #                 num_grids[n,level - 1] += 1 
        #                 mx = int(q_file.readline().split()[0])
        #                 my = int(q_file.readline().split()[0])
        #                 num_cells[n,level - 1] += mx * my

        # Compare gauge data
        data = numpy.loadtxt(os.path.join(self.temp_path, 'fort.gauge'))
        data_sum = [data[:,2].sum(), data[:,3].sum()]

        regression_data_file = os.path.join(self.test_path, "regression_data.txt")
        if save:
            numpy.savetxt(regression_data_file, data)
        regression_data = numpy.loadtxt(regression_data_file)
        regression_sum = [regression_data[:,2].sum(), regression_data[:,3].sum()]

        # Compare data
        tolerance = 1e-14
        assert numpy.allclose(data_sum, regression_sum, tolerance), \
                "\n data: %s, \n expected: %s" % (data_sum, regression_sum)
        # assert numpy.allclose(data, regression_data, tolerance), \
        #         "Full gauge match failed."

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