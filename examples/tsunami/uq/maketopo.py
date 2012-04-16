"""
Create topo files needed for this example:
    etopo10min120W60W60S0S.asc        download from GeoClaw topo repository
    usgs100227.tt1                    create using Okada model and .cfg file
    
"""

from geoclaw import topotools
import os,sys
from numpy import where,sqrt,cos,pi

def maketopo():
    """
    Output topography file for the entire domain and for a coastal feature.
    """
    nxpoints = 101
    nypoints = 151
    xlower = -0.95
    xupper = 1.05
    ylower = -1.
    yupper = 2.
    outfile= "domain.topotype2"     
    topotools.topo2writer(outfile,topo,xlower,xupper,ylower,yupper,nxpoints,nypoints)

    nxpoints = 81
    nypoints = 101
    xlower = 0.99
    xupper = 1.03
    ylower = 0.97
    yupper = 1.02
    outfile= "coast.topotype2"     
    topotools.topo2writer(outfile,topo,xlower,xupper,ylower,yupper,nxpoints,nypoints)

def cosine_bump(x,y, x0, y0, Radius, Height):
    """
    A cosine bump with Radius and Height centered at (x0,y0).
    """
    r = sqrt((x-x0)**2 + (y-y0)**2)
    z = where(r > Radius, 0., 0.5*Height*(1. + cos(pi*r/Radius)))
    return z

def topo(x,y):
    """
    Shelf and linear beach, with bay and hill.
    """
    # shelf and linear beach:
    z = where(x<0.8, -100., 500.*(x-1.))
    # add bay:
    z = z + cosine_bump(x,y, x0=1.01, y0=1., Radius=.02, Height=-5)
    # add hill:
    z = z + cosine_bump(x,y, x0=1.0, y0=0.985, Radius=.01, Height=5)
    return z
    
def makedtopo():
    """
    Create dtopo data file for deformation of sea floor due to earthquake.
    Uses the Okada model with fault parameters and mesh specified in the
    .cfg file.
    """
    from geoclaw import okada
    dtopo_fname = 'fault1.tt1'
    dtopo_cfg = 'fault1.cfg'
    print "Using Okada model to create %s " % dtopo_fname
    okada.builddynamicdeffile(dtopo_cfg, dtopo_cfg, dtopo_fname)

    dtopo_fname = 'fault2.tt1'
    dtopo_cfg = 'fault2.cfg'
    print "Using Okada model to create %s " % dtopo_fname
    okada.builddynamicdeffile(dtopo_cfg, dtopo_cfg, dtopo_fname)


if __name__=='__main__':
    maketopo()
    makedtopo()
