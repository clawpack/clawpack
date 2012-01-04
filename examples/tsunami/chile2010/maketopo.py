"""
Create topo files needed for this example:
    etopo10min120W60W60S0S.asc        download from GeoClaw topo repository
    usgs100227.tt1                    create using Okada model and .cfg file
    
"""

from pyclaw.geotools import topotools
import os,sys

def gettopo():
    """
    Retrieve the topo file from the GeoClaw repository.
    """
    remote_directory = 'http://kingkong.amath.washington.edu/topo/etopo2'
    topo_fname = 'etopo10min120W60W60S0S.asc'
    topotools.get_topo(topo_fname, remote_directory)

    
def makedtopo():
    """
    Create dtopo data file for deformation of sea floor due to earthquake.
    Uses the Okada model with fault parameters and mesh specified in the
    .cfg file.
    """
    from pyclaw.geotools import okada
    dtopo_fname = 'usgs100227.tt1'
    dtopo_cfg = 'usgs100227.cfg'
    if os.path.exists(dtopo_fname):
        print "*** Not regenerating dtopo file (already exists): %s" % dtopo_fname
    else:
        print "Using Okada model to create %s " % dtopo_fname
        okada.builddynamicdeffile(dtopo_cfg, dtopo_cfg, dtopo_fname)



def makeqinit():
    """
    Create qinit data file
    """
    nxpoints=100
    nypoints=100
    xlower=-84.e0
    xupper=-80.e0
    ylower=-49.e0
    yupper=-45.e0
    outfile= "hump.xyz"
    topotools.topo1writer(outfile,qinit,xlower,xupper,ylower,yupper,nxpoints,nypoints)


def qinit(x,y):
    """
    Gaussian hump:
    """
    from numpy import where,exp
    x0 = -82.0
    y0 = -47.0
    d = topotools.gcdist(x,y,x0,y0)  # great circle distance on earth
    ze = -d**2/2.e9
    z = where(ze>-100., 20.e0*exp(ze), 0.)
    return z

if __name__=='__main__':
    gettopo()
    #makeqinit()
    makedtopo()
