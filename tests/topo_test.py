#!/usr/bin/env python
# encoding: utf-8

from __future__ import absolute_import
from __future__ import print_function
import os
import glob
import tempfile


import matplotlib.pyplot as plt
from clawpack.geoclaw import topotools
import clawpack.geoclaw.topo as topo
import numpy as np


def test1():
    """
    Make two topo files and then read them in and plot them.
    The second file is a finer grid over a smaller region.
    """

    fname = 'bowl.tt2'
    maketopo1a(fname)

    tpd = topotools.TopoPlotData(fname)
    tpd.imshow = True
    tpd.cmin = -1000.
    tpd.cmax = 2000.
    tpd.addcolorbar = True
    tpd.plot()

    fname = 'hill.tt2'
    maketopo1b(fname)

    tpd = topotools.TopoPlotData(fname)
    tpd.imshow = True
    tpd.cmin = -1000.
    tpd.cmax = 2000.
    tpd.addcolorbar = False
    tpd.plot()

    plt.title('Bathymetry / topography')
    fname = 'topotest1.png'
    plt.savefig(fname)
    # print 'Created ',fname


def topo1(x,y):
    """
    Sample topography
    """
    # Parabolic bowl
    z = 1000.*(x**2 + y**2 - 1.)
    # Add a Gaussian hill
    z = z + 1000.*np.exp(-100*((x-0.7)**2 + (y-0.8)**2))
    return z


def maketopo1a(path):
    """
    Output topography file for the entire domain
    """

    nxpoints = 101
    nypoints = 76
    xlower = -1.5
    xupper = 2.5
    ylower = -1.
    yupper = 2.
    topotools.topo2writer(path,topo1,xlower,xupper,ylower,yupper,\
                          nxpoints,nypoints)

def maketopo1b(path):
    """
    Output topography file for the entire domain
    """
    nxpoints = 101
    nypoints = 71
    xlower = 0.0
    xupper = 1.0
    ylower = 0.5
    yupper = 1.2
    topotools.topo2writer(path,topo1,xlower,xupper,ylower,yupper,\
                          nxpoints,nypoints)


def test_topography_object(plot=False):
    """
    Test the Topography object's functionality
    """

    try:
        base_path = tempfile.mkdtemp()

        # Create initial test bathymetry
        maketopo1a(os.path.join(base_path, 'bowl.tt2'))
        maketopo1b(os.path.join(base_path, 'hill.tt2'))

        hill_topo = []
        bowl_topo = []
        topo_types = [2,3,1,2]
        for (n, topo_type) in enumerate(topo_types):
            bowl_path = os.path.join(base_path, 'bowl.tt%s' % topo_type)
            hill_path = os.path.join(base_path, 'hill.tt%s' % topo_type)

            bowl_topo.append(topo.Topography(bowl_path))
            hill_topo.append(topo.Topography(hill_path))

            if plot:
                fig = plt.figure()
                axes = fig.add_subplot(1,1,1)
                hill_topo[-1].plot(axes=axes, limits=[-2000,0])
                bowl_topo[-1].plot(axes=axes, limits=[-2000,0])
                fig.suptitle('Bathymetry / topography, topo type = %s' % topo_type)
                plt.savefig('topotest%s.png' % (n + 2))

            print(n, topo_type)
            if n + 1 != len(topo_types):
                bowl_path = os.path.join(base_path, 'bowl.tt%s' % topo_types[n+1])
                hill_path = os.path.join(base_path, 'hill.tt%s' % topo_types[n+1])
                bowl_topo[-1].write(bowl_path)
                hill_topo[-1].write(hill_path)

        # Check data
        for (n,topo_type) in enumerate(topo_types):
            for (m,topo_type) in enumerate(topo_types):
                assert np.all(bowl_topo[n].X == bowl_topo[m].X), \
                       "bowl[%s].X != bowl[%s].X" % (n,m)
                assert np.all(bowl_topo[n].Y == bowl_topo[m].Y), \
                       "bowl[%s].Y != bowl[%s].Y" % (n,m)
                assert np.all(bowl_topo[n].Z == bowl_topo[m].Z), \
                       "bowl[%s].Z != bowl[%s].Z" % (n,m)
                assert np.all(hill_topo[n].X == hill_topo[m].X), \
                       "hill[%s].X != hill[%s].X" % (n,m)
                assert np.all(hill_topo[n].Y == hill_topo[m].Y), \
                       "hill[%s].Y != hill[%s].Y" % (n,m)
                assert np.all(hill_topo[n].Z == hill_topo[m].Z), \
                       "hill[%s].Z != hill[%s].Z" % (n,m)

    finally:
        paths = glob.glob(os.path.join(base_path,"*"))
        for path in paths:
            os.remove(path)
        os.rmdir(base_path)


if __name__=='__main__':
    print("Starting procedural test...")
    test1()
    print("Done performing procedural test.")

    print("Starting object test...")
    test_topography_object(plot=True)
    print("Done performing object test...")