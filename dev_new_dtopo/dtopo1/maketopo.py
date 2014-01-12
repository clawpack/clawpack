
"""
Module to create topo and qinit data files for this example.
"""

from clawpack.geoclaw.topotools import topo1writer, topo2writer
from numpy import *


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
    topo2writer(outfile,topo,xlower,xupper,ylower,yupper,nxpoints,nypoints)

def maketopo2():
    """
    Another topo file for testing
    """
    nxpoints=101
    nypoints=101
    xlower = -10e0
    xupper= 10e0
    ylower = -10e0
    yupper= 10e0
    outfile= "topo2.topotype2"
    topo2writer(outfile,topo2,xlower,xupper,ylower,yupper,nxpoints,nypoints)


def topo(x,y):
    z = -10*ones(x.shape)
    return z

def topo2(x,y):
    z = -2000*ones(x.shape)
    return z



if __name__=='__main__':
    maketopo()
    #maketopo2()
    #makedtopo()
