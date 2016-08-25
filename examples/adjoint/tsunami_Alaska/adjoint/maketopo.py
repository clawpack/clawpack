"""
    Download topo and dtopo files needed for this example.
    
    Call functions with makeplots==True to create plots of topo, slip, and dtopo.
"""

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
xcenter = 235.80917e0
radius = 1.0e0
ycenter =41.74111e0

def get_topo(makeplots=False):
    """
        Retrieve the topo file from online.
    """
    from clawpack.geoclaw import topotools
    topo_fname = 'etopo1min170E124W40N61N.asc'
    url = 'http://students.washington.edu/bndavis/misc/topo/' + topo_fname
    clawpack.clawutil.data.get_remote_file(url, output_dir=scratch_dir,
            file_name=topo_fname, verbose=True)
        
    topo_fname = 'etopo4min120E110W0N62N.asc'
    url = 'http://students.washington.edu/bndavis/misc/topo/' + topo_fname
    clawpack.clawutil.data.get_remote_file(url, output_dir=scratch_dir,
            file_name=topo_fname, verbose=True)
                                           
    if makeplots:
        from matplotlib import pyplot as plt
        topo = topotools.Topography(topo_fname, topo_type=2)
        topo.plot()
        fname = os.path.splitext(topo_fname)[0] + '.png'
        plt.savefig(fname)
        print "Created ",fname

def makeqinit():
    """
        Create qinit data file
    """
    
    nxpoints = 201
    nypoints = 201
    
    xlower = 140.0e0
    xupper = 250.0e0
    ylower = 10.0e0
    yupper = 62.0e0
    
    outfile= "hump.xyz"     
    topotools.topo1writer(outfile,qinit,xlower,xupper,ylower,yupper,nxpoints,nypoints)

def qinit(x,y):
    from numpy import where
    ze = sqrt((x-xcenter)**2 + (y-ycenter)**2)
    z = where(ze<radius, 1.0e0, 0.e0)
    return z

if __name__=='__main__':
    get_topo(False)
    makeqinit()
