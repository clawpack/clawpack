#!/usr/bin/env python
# encoding: utf-8

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
    print 'Created ',fname


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


def test_topography_object():
    """
    Test the Topography object's functionality
    """

    try:
        base_path = tempfile.mkdtemp()

        # Test basic reading of topo type 2
        topo_files = []
        paths = [os.path.join(base_path,'bowl.tt2'), 
                 os.path.join(base_path,'hill.tt2')]
        maketopo1a(paths[0])
        topo_files.append(topo.Topography(paths[0]))
        maketopo1b(paths[1])
        topo_files.append(topo.Topography(paths[1]))

        fig = plt.figure()

        axes = fig.add_subplot(1,1,1)
        topo_files[1].plot(axes=axes, limits=[-2000,0])
        topo_files[0].plot(axes=axes, limits=[-2000,0])
        fig.suptitle('Bathymetry / topography')
        plt.savefig('topotest2.png')

        paths[0] = os.path.join(base_path,'bowl.tt3')
        paths[1] = os.path.join(base_path,'hill.tt3')
        topo_files[0].write(paths[0])
        topo_files[1].write(paths[1])

        del topo_files

        # Test topo type 3
        topo_files = []
        topo_files.append(topo.Topography(paths[0]))
        topo_files.append(topo.Topography(paths[1]))


        fig = plt.figure()

        axes = fig.add_subplot(1,1,1)
        topo_files[1].plot(axes=axes, limits=[-2000,0])
        topo_files[0].plot(axes=axes, limits=[-2000,0])
        fig.suptitle('Bathymetry / topography')
        plt.savefig('topotest3.png')

        paths[0] = os.path.join(base_path,'bowl.tt1')
        paths[1] = os.path.join(base_path,'hill.tt1')
        topo_files[0].write(paths[0])
        topo_files[1].write(paths[1])

        del topo_files

        # Test topo type 1
        topo_files = []
        topo_files.append(topo.Topography(paths[0]))
        topo_files.append(topo.Topography(paths[1]))

        fig = plt.figure()

        axes = fig.add_subplot(1,1,1)
        topo_files[1].plot(axes=axes, limits=[-2000,0])
        topo_files[0].plot(axes=axes, limits=[-2000,0])
        fig.suptitle('Bathymetry / topography')
        plt.savefig('topotest4.png')

        del topo_files

    finally:
        paths = glob.glob(os.path.join(base_path,"*"))
        for path in paths:
            os.remove(path)
        os.rmdir(base_path)


if __name__=='__main__':
    print "Starting procedural test..."
    test1()
    print "Done performing procedural test."

    print "Starting object test..."
    test_topography_object()
    print "Done performing object test..."