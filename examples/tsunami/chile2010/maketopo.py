"""
Create topo files needed for this example:
    etopo10min120W60W60S0S.asc        download from GeoClaw topo repository
    usgs100227.tt1                    create using Okada model and .cfg file
    
"""

from clawpack.geoclaw import topotools
import os,sys

def gettopo():
    """
    Retrieve the topo file from the GeoClaw repository.
    """
    remote_directory = 'http://www.clawpack.org/geoclaw/topo/etopo'
    topo_fname = 'etopo10min120W60W60S0S.asc'
    topotools.get_topo(topo_fname, remote_directory)

    
def makedtopo():
    """
    Create dtopo data file for deformation of sea floor due to earthquake.
    Uses the Okada model with fault parameters and mesh specified in the
    .cfg file.
    """
    from clawpack.geoclaw import okada
    dtopo_fname = 'usgs100227.tt1'
    dtopo_cfg = 'usgs100227.cfg'
    if os.path.exists(dtopo_fname):
        print "*** Not regenerating dtopo file (already exists): %s" % dtopo_fname
    else:
        print "Using Okada model to create %s " % dtopo_fname
        okada.builddynamicdeffile(dtopo_cfg, dtopo_cfg, dtopo_fname)


if __name__=='__main__':
    gettopo()
    makedtopo()
