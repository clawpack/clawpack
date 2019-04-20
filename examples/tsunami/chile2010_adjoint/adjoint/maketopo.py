"""
    Download topo and dtopo files needed for this example.
    
    Call functions with makeplots==True to create plots of topo, slip, and dtopo.
"""

from __future__ import print_function
import os,sys
import clawpack.clawutil.data
from clawpack.geoclaw import topotools
from numpy import *

try:
    CLAW = os.environ['CLAW']
except:
    raise Exception("*** Must first set CLAW enviornment variable")

# Scratch directory for storing topo and dtopo files:
scratch_dir = os.path.join(CLAW, 'geoclaw', 'scratch')


# Initial data for adjoint is Gaussian hump around this location:
# DART 32412 location:
xcenter = -86.392
ycenter = -17.975


def get_topo(makeplots=False):
    """
    Retrieve the topo file from the GeoClaw repository.
    """
    from clawpack.geoclaw import topotools
    topo_fname = 'etopo10min120W60W60S0S.asc'
    url = 'http://depts.washington.edu/clawpack/geoclaw/topo/etopo/' + topo_fname
    clawpack.clawutil.data.get_remote_file(url, output_dir=scratch_dir, 
            file_name=topo_fname, verbose=True)

    if makeplots:
        from matplotlib import pyplot as plt
        topo = topotools.Topography(os.path.join(scratch_dir,topo_fname), 
                                    topo_type=2)
        topo.plot()
        fname = os.path.splitext(topo_fname)[0] + '.png'
        plt.savefig(fname)
        print("Created ",fname)


def makeqinit():
    """
        Create qinit data file
    """
    
    nxpoints = 201
    nypoints = 201
    
    xlower = xcenter - 1.5
    xupper = xcenter + 1.5
    ylower = ycenter - 1.5
    yupper = ycenter + 1.5
    
    outfile= "hump.xyz"     
    topotools.topo1writer(outfile,qinit,xlower,xupper,ylower,yupper,nxpoints,nypoints)

def qinit(x,y):
    from numpy import where
    from clawpack.geoclaw.util import haversine

    #radius = 1.0e0
    #ze = sqrt((x-xcenter)**2 + (y-ycenter)**2)
    #z = where(ze<radius, 1.0e0, 0.e0)

    # Gaussian using distance in meters:
    r = haversine(x,y,xcenter,ycenter)
    z = exp(-(r/20e3)**2)
    return z

if __name__=='__main__':
    get_topo(False)
    makeqinit()
