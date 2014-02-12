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

import sys
import numpy

import batch

import clawpack.geoclaw.topotools as tt

import setrun

class PlaneWaveJob(batch.Job):

    r"""Base plane wave test case.

    """

    def __init__(self, angle=0.0, bathy_angle=0.0, location=[-0.1, 0.0], 
                       bathy_location=0.15):

        super(PlaneWaveJob, self).__init__()

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

    def __str__(self):
        output = super(PlaneWaveJob, self).__str__()
        output += "  Angle = %s\n" % self.rundata.qinit_data.angle
        output += "  Bathy Angle = %s\n" % self.bathy_angle
        return output

    def write_data_objects(self):
        super(PlaneWaveJob, self).write_data_objects()

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

if __name__ == "__main__":

    tests = []
    tests.append(PlaneWaveJob())
    tests.append(PlaneWaveJob(numpy.pi / 4.0, numpy.pi / 4.0))
    tests.append(PlaneWaveJob(numpy.pi / 4.0, 0.0))
    tests.append(PlaneWaveJob(0.0, numpy.pi / 4.0))
    tests.append(PlaneWaveJob(numpy.pi / 8.0, 0.0))
    tests.append(PlaneWaveJob(numpy.pi / 4.0, numpy.pi / 8.0))

    if len(sys.argv) > 1:
        if sys.argv[1].lower() == 'all':
            tests_to_run = tests
        else:
            tests_to_run = []
            for test in sys.argv[1:]:
                tests_to_run.append(tests[int(test)])

        controller = batch.BatchController(tests)
        controller.run()

    else:
        controller = batch.BatchController(tests)
        print controller


