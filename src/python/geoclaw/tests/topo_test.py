
import pylab
from clawpack.geoclaw import topotools
import numpy as np


def test1():
    """
    Make two topo files and then read them in and plot them.
    The second file is a finer grid over a smaller region.
    """

    maketopo1a()

    fname = 'bowl.tt2'
    tpd = topotools.TopoPlotData(fname)
    tpd.imshow = True
    tpd.cmin = -1000.
    tpd.cmax = 2000.
    tpd.addcolorbar = True
    tpd.plot()

    maketopo1b()

    fname = 'hill.tt2'
    tpd = topotools.TopoPlotData(fname)
    tpd.imshow = True
    tpd.cmin = -1000.
    tpd.cmax = 2000.
    tpd.addcolorbar = False
    tpd.plot()

    pylab.title('Bathymetry / topography')
    fname = 'topotest1.png'
    pylab.savefig(fname)
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


def maketopo1a():
    """
    Output topography file for the entire domain
    """
    nxpoints = 101
    nypoints = 76
    xlower = -1.5
    xupper = 2.5
    ylower = -1.
    yupper = 2.
    outfile = "bowl.tt2"
    topotools.topo2writer(outfile,topo1,xlower,xupper,ylower,yupper,\
                          nxpoints,nypoints)

def maketopo1b():
    """
    Output topography file for the entire domain
    """
    nxpoints = 101
    nypoints = 71
    xlower = 0.0
    xupper = 1.0
    ylower = 0.5
    yupper = 1.2
    outfile = "hill.tt2"
    topotools.topo2writer(outfile,topo1,xlower,xupper,ylower,yupper,\
                          nxpoints,nypoints)



if __name__=='__main__':
    test1()
