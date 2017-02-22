
"""
Module to create topo and qinit data files for this example.
"""

from __future__ import absolute_import
from __future__ import print_function
from clawpack.geoclaw.topotools import Topography

from clawpack.geoclaw import dtopotools
import numpy
import matplotlib.pyplot as plt
import os


def maketopo():
    """
    Output topography file for the entire domain
    """
    nxpoints=201
    nypoints=201
    xlower = -10e0
    xupper= 10e0
    ylower = -10e0
    yupper= 10e0
    outfile= "topo1.topotype2"
    topography = Topography(topo_func=topo)
    topography.x = numpy.linspace(xlower,xupper,nxpoints)
    topography.y = numpy.linspace(ylower,yupper,nypoints)
    topography.write(outfile, topo_type=2, Z_format="%22.15e")
    dx = (xupper-xlower)/(nxpoints-1)
    dy = (yupper-ylower)/(nypoints-1)
    print("==> topo1 has dx = %g, dy = %g" % (dx,dy))

def maketopo2():
    """
    Another topo file for testing
    """
    nxpoints=21
    nypoints=51
    xlower = -0.5
    xupper= -0.3
    ylower = -0.1
    yupper= 0.4
    outfile= "topo2.topotype2"
    topography = Topography(topo_func=topo2)
    topography.x = numpy.linspace(xlower,xupper,nxpoints)
    topography.y = numpy.linspace(ylower,yupper,nypoints)
    topography.write(outfile, topo_type=2, Z_format="%22.15e")
    dx = (xupper-xlower)/(nxpoints-1)
    dy = (yupper-ylower)/(nypoints-1)
    print("==> topo2 has dx = %g, dy = %g" % (dx,dy))

h0 = 1000.

def topo(x,y):
    #z = -10*numpy.ones(x.shape)
    z = -h0*(1 + 0.5*numpy.cos(x-y))
    return z

def topo2(x,y):
    #z = -15*numpy.ones(x.shape)
    z = -h0*(1. + numpy.exp(x+y))
    return z


def read_fault(fname_subfaults, plotfig=None):
    """
    Test data
    """

    input_units = {'slip': 'm', 'depth': 'km', 'length': 'km', 'width': 'km'}

    fault = dtopotools.CSVFault()
    fault.read(fname_subfaults, input_units=input_units, 
                coordinate_specification="top center")
    fault.rupture_type = 'dynamic'
                    
    if plotfig:
        plt.figure(plotfig)
        fault.plot_subfaults(slip_color=True,plot_rake=True)
        
    return fault


def make_dtopo1(plotfig=None):
    """
    Test data.
    """

    fname_subfaults = 'dtopo1.csv'
    fault = read_fault(fname_subfaults)

    dtopo_fname = fname_subfaults.split('.')[0] + '.tt3'

    print("Using Okada model to create %s " % dtopo_fname)

    # Needed for extent of dtopo file:
    xlower = -0.4
    xupper = 0.6
    ylower = -0.4
    yupper = 0.4

    # number of grid points in dtopo file
    mx = 151
    my = 121

    dx = (xupper-xlower)/(mx-1)
    dy = (yupper-ylower)/(my-1)
    print("==> dtopo1 has dx = %g, dy = %g" % (dx,dy))
    x = numpy.linspace(xlower,xupper,mx)
    y = numpy.linspace(ylower,yupper,my)

    tfinal = 1.
    times = numpy.linspace(0.,tfinal,25)
    dtopo = fault.create_dtopography(x,y,times=times)
    dtopo.write(dtopo_fname,3)


    if plotfig:
        plt.figure(plotfig)
        # plot final deformation:
        dtopo.plot_dz_colors(tfinal, cmax_dz=16,dz_interval=1)
        # make animation as html:
        dtopo.animate_dz_colors(times, cmax_dz=16,dz_interval=1,style='html')

    return dtopo


def make_dtopo2(plotfig=None):
    """
    Test data.
    """

    fname_subfaults = 'dtopo2.csv'
    fault = read_fault(fname_subfaults)

    dtopo_fname = fname_subfaults.split('.')[0] + '.tt3'

    print("Using Okada model to create %s " % dtopo_fname)

    # Needed for extent of dtopo file:
    xlower = -0.9
    xupper = 0.1
    ylower = -0.4
    yupper = 0.4

    # number of grid points in dtopo file
    mx = 201
    my = 161

    dx = (xupper-xlower)/(mx-1)
    dy = (yupper-ylower)/(my-1)
    print("==> dtopo2 has dx = %g, dy = %g" % (dx,dy))


    x = numpy.linspace(xlower,xupper,mx)
    y = numpy.linspace(ylower,yupper,my)

    tfinal = 1.2
    times = numpy.linspace(0.5,tfinal,25)
    dtopo = fault.create_dtopography(x,y,times=times)
    dtopo.write(dtopo_fname,3)


    if plotfig:
        plt.figure(plotfig)
        # plot final deformation:
        dtopo.plot_dz_colors(tfinal, cmax_dz=16,dz_interval=1)
        # make animation as html:
        dtopo.animate_dz_colors(times, cmax_dz=16,dz_interval=1,style='html')

        return dtopo

if __name__=='__main__':
    maketopo()
    maketopo2()
    make_dtopo1()
    make_dtopo2()
