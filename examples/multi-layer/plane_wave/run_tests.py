#!/usr/bin/env python
# encoding: utf-8

r"""Contains plane wave tests for the multi-layer SWE code

Can run a particular test by giving the number of the test at the comand line.

:Note:

This test suite requires the `batch` package.

"""

# ============================================================================
#  Distributed under the terms of the Berkeley Software Distribution (BSD) 
#  license
#                     http://www.opensource.org/licenses/
# ============================================================================

from __future__ import absolute_import
from __future__ import print_function
import sys
import numpy

import batch

import clawpack.geoclaw.topotools as tt

import setrun

class PlaneWaveTest(batch.Job):

    r"""Base plane wave test case.

    """

    def __init__(self, angle=0.0, bathy_angle=0.0, location=[-0.1, 0.0], 
                       bathy_location=0.15):

        super(PlaneWaveTest, self).__init__()

        self.executable = 'xgeoclaw'
        self.type = "multilayer"
        self.name = "planewave"

        # Convert angle to degrees for the label
        self.prefix = "ml_2d_ia%s_ba%s" % (int(angle * 180.0 / numpy.pi),
            int(bathy_angle * 180.0 / numpy.pi))

        # Data objects
        self.rundata = setrun.setrun()

        self.rundata.qinit_data.angle = angle
        self.rundata.qinit_data.init_location = location

        self.bathy_location = bathy_location
        self.bathy_angle = bathy_angle

        # Add gauges down perpendicular
        self.rundata.gaugedata.gauges = []
        gauge_locations = [-0.1,0.0,0.1,0.2,0.3]
        for (i,x_c) in enumerate(gauge_locations):
            # y0 = (self.run_data.clawdata.yupper - self.run_data.clawdata.ylower) / 2.0
            # x_p,y_p = transform_c2p(x_c,0.0,location[0],location[1],angle)
            x_p = x_c * numpy.cos(angle)
            y_p = x_c * numpy.sin(angle)
            # print "+=====+"
            # print x_c,0.0
            # print x_p,y_p
            if (self.rundata.clawdata.lower[0] < x_p < self.rundata.clawdata.upper[0] and
                    self.rundata.clawdata.lower[1] < y_p < self.rundata.clawdata.upper[1]):
                self.rundata.gaugedata.gauges.append([i, x_p, y_p, 0.0, 1e10])
                # print "Gauge %s: (%s,%s)" % (i,x_p,y_p)
        # print "+=====+"

        # self.setplot = lambda plotdata:setplot.setplot(
        #                                         plotdata, 
        #                                         bathy_location=bathy_location, 
        #                                         bathy_angle=bathy_angle)

    def __str__(self):
        output = super(PlaneWaveTest, self).__str__()
        output += "  Angle = %s\n" % self.rundata.qinit_data.angle
        output += "  Bathy Angle = %s\n" % self.bathy_angle
        return output

    def write_data_objects(self):
        super(PlaneWaveTest, self).write_data_objects()

        # Write out bathy file
        mx = self.rundata.clawdata.num_cells[0]
        my = self.rundata.clawdata.num_cells[1]
        xlower = self.rundata.clawdata.lower[0]
        xupper = self.rundata.clawdata.upper[0]
        ylower = self.rundata.clawdata.lower[1]
        yupper = self.rundata.clawdata.upper[1]
        dx = (xupper - xlower) / mx
        dy = (yupper - ylower) / my
        d = min(dx,dy)
        mx = int((xupper - xlower) / d) + 8
        my = int((yupper - ylower) / d) + 8
        
        xlower = xlower - d*4.0
        ylower = ylower - d*4.0
        xupper = xupper + d*4.0
        yupper = yupper + d*4.0

        step = lambda x,y: setrun.bathy_step(x, y, location=self.bathy_location,
                                            angle=self.bathy_angle)
        
        tt.topo2writer('./topo.tt2', step, xlower, xupper, ylower, yupper, 
                                     mx, my, nodata_value=-99999)

        # Write out simple bathy geometry file for communication to the plotting
        with open("./bathy_geometry.data", 'w') as bathy_geometry_file:
            bathy_geometry_file.write("%s\n%s" % (self.bathy_location,
                                                  self.bathy_angle) )


class BubbleTest(PlaneWaveTest):

    r"""Jump bathymetry with gaussian bump initial conditions."""

    def __init__(self, eigen_method=1):

        super(BubbleTest, self).__init__()

        self.type = "multilayer"
        self.name = "bubble"
        self.prefix = "ml_e%s" % eigen_method

        # Data objects
        self.rundata = setrun.setrun()

        self.rundata.clawdata.lower = [0.0, 0.0]
        self.rundata.clawdata.upper = [1.0, 1.0]
        self.rundata.clawdata.num_cells = [100, 100]
        self.rundata.multilayer_data.eigen_method = eigen_method
        self.rundata.multilayer_data.eta = [0.0, -0.5]

        self.rundata.qinit_data.qinit_type = 7
        self.rundata.qinit_data.init_location = [-0.25,0.0]
        self.rundata.qinit_data.wave_family = 0
        self.rundata.qinit_data.epsilon = 0.4
        self.rundata.qinit_data.sigma = 0.08

        self.bathy_location = 0.8


    def __str__(self):
        output = super(PlaneWaveTest, self).__str__()
        output += "  Eigen Method = %s\n" % self.rundata.multilayer_data.eigen_method
        return output


# Defintion of tests defined here
tests = []
tests.append(PlaneWaveTest())
tests.append(PlaneWaveTest(numpy.pi / 4.0, numpy.pi / 4.0))
tests.append(PlaneWaveTest(numpy.pi / 4.0, 0.0))
tests.append(PlaneWaveTest(0.0, numpy.pi / 4.0))
tests.append(PlaneWaveTest(numpy.pi / 8.0, 0.0))
tests.append(PlaneWaveTest(numpy.pi / 4.0, numpy.pi / 8.0))
for method in [1,2,3,4]:
    tests.append(BubbleTest(eigen_method=method))


if __name__ == "__main__":

    if len(sys.argv) > 1:
        if sys.argv[1].lower() == 'all':
            tests_to_run = tests
        else:
            tests_to_run = []
            for test in sys.argv[1:]:
                tests_to_run.append(tests[int(test)])

        controller = batch.BatchController(tests_to_run)
        print(controller)
        controller.run()

    else:
        controller = batch.BatchController(tests)
        print(controller)


