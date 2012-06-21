#!/usr/bin/env python
# encoding: utf-8
r"""
Contains all test suites for the 1d multi-layer code

Can run a particular test by giving the number of the test at the comand line

:Authors:
    Kyle Mandli (2011-05-11) Initial version
"""
# ============================================================================
#      Copyright (C) 2011 Kyle Mandli <mandli@amath.washington.edu>
#
#  Distributed under the terms of the Berkeley Software Distribution (BSD) 
#  license
#                     http://www.opensource.org/licenses/
# ============================================================================

import os
import sys
import numpy as np

import topotools as tt
import test_runs

tests = []

def transform_c2p(x,y,x0,y0,theta):
    return ((x+x0)*np.cos(theta) - (y+y0)*np.sin(theta),
            (x+x0)*np.sin(theta) + (y+y0)*np.cos(theta))

def transform_p2c(x,y,x0,y0,theta):
    return ( x*np.cos(theta) + y*np.sin(theta) - x0,
            -x*np.sin(theta) + y*np.cos(theta) - y0)

    
# Base Plane Wave Test
class PlaneWaveBaseTest(test_runs.TestML2D):
    
    # Domain centered at (0.0,0.0), rotate around this point
    
    def __init__(self,angle=0.0,bathy_angle=0.0,location=[-0.1,0.0],bathy_location=0.15):
        super(PlaneWaveBaseTest,self).__init__()

        self.name = "plane_wave"
        
        # self.run_data.clawdata.nout = 1
        # self.run_data.clawdata.tfinal = 0.1
        self.run_data.clawdata.mx = 3*100
        self.run_data.clawdata.my = 3*100
        
        self.ml_data.eigen_method = 2
        self.ml_data.inundation_method = 2
        
        self.ml_data.init_type = 2
        self.ml_data.angle = angle
        self.ml_data.init_location = location
        
        self.ml_data.bathy_type = 1
        self.ml_data.bathy_location = bathy_location
        self.ml_data.bathy_angle = bathy_angle
        self.ml_data.bathy_left = -1.0
        self.ml_data.bathy_right = -0.2
        
        # Add gauges down perpendicular
        gauge_locations = [-0.1,0.0,0.1,0.2,0.3]
        for (i,x_c) in enumerate(gauge_locations):
            # y0 = (self.run_data.clawdata.yupper - self.run_data.clawdata.ylower) / 2.0
            # x_p,y_p = transform_c2p(x_c,0.0,location[0],location[1],angle)
            x_p = x_c*np.cos(angle)
            y_p = x_c*np.sin(angle)
            print "+=====+"
            print x_c,0.0
            print x_p,y_p
            if (self.run_data.clawdata.xlower < x_p < self.run_data.clawdata.xupper and
                    self.run_data.clawdata.ylower < y_p < self.run_data.clawdata.yupper):
                self.run_data.geodata.gauges.append([i, x_p, y_p, 0.0, 1e10])
                print "Gauge %s: (%s,%s)" % (i,x_p,y_p)
        print "+=====+"
        # Convert angle to degrees for the label
        self.prefix = "ml_2d_ia%s_ba%s" % (int(angle * 180.0 / np.pi),
            int(bathy_angle * 180.0 / np.pi))
        
    def write_data_objects(self):
        super(PlaneWaveBaseTest,self).write_data_objects()
        
        # Write out bathy file
        mx = self.run_data.clawdata.mx
        my = self.run_data.clawdata.my
        xlower = self.run_data.clawdata.xlower
        xupper = self.run_data.clawdata.xupper
        ylower = self.run_data.clawdata.ylower
        yupper = self.run_data.clawdata.yupper
        dx = (xupper - xlower) / mx
        dy = (yupper - ylower) / my
        d = min(dx,dy)
        mx = int((xupper - xlower) / d) + 8
        my = int((yupper - ylower) / d) + 8
        
        xlower = xlower - d*4.0
        ylower = ylower - d*4.0
        xupper = xupper + d*4.0
        yupper = yupper + d*4.0
        
        def step(x,y):
            x_c,y_c = transform_p2c(x,y,self.ml_data.bathy_location,0.0,
                                    self.ml_data.bathy_angle)
            return ((x_c <= 0.0) * self.ml_data.bathy_left 
                  + (x_c >  0.0) * self.ml_data.bathy_right)
        
        tt.topo2writer('./topo.data',step,
                xlower,xupper,ylower,yupper,mx,my,nodata_value=-99999)
                
               
# Angle tests
# tests.append(PlaneWaveBaseTest(0.0,location=[0.1,0.0]))
# tests.append(PlaneWaveBaseTest(np.pi/4.0,np.pi/4.0,location=[0.5,0.0],
#                                                    bathy_location=0.9))
# 
# tests.append(PlaneWaveBaseTest(np.pi/4.0,0.0,location=[0.4,0.0],
#                                              bathy_location=0.6))
# tests.append(PlaneWaveBaseTest(0.0,np.pi/4.0,location=[0.1,0.5],
#                                              bathy_location=0.8))
# 
# tests.append(PlaneWaveBaseTest(np.pi/8.0,0.0,location=[0.3,0.3]))
# tests.append(PlaneWaveBaseTest(np.pi/4.0,np.pi/8.0,location=[0.4,0.0],
#                                                    bathy_location=0.7))
tests.append(PlaneWaveBaseTest())
tests.append(PlaneWaveBaseTest(np.pi/4.0,np.pi/4.0))

tests.append(PlaneWaveBaseTest(np.pi/4.0,0.0))
tests.append(PlaneWaveBaseTest(0.0,np.pi/4.0))

tests.append(PlaneWaveBaseTest(np.pi/8.0,0.0))
tests.append(PlaneWaveBaseTest(np.pi/4.0,np.pi/8.0))


if __name__ == '__main__':
    if len(sys.argv) > 1:
        if sys.argv[1].lower() == 'all':
            tests_to_be_run = tests
        else:
            tests_to_be_run = []
            for test in sys.argv[1:]:
                tests_to_be_run.append(tests[int(test)])
        
        test_runs.run_tests(tests_to_be_run,parallel=True)
    
    else:
        test_runs.print_tests(tests)
