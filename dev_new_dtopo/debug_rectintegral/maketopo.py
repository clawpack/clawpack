
"""
Module to create topo and qinit data files for this example.
"""

from clawpack.geoclaw.topotools import topo1writer, topo2writer
from numpy import *

#from pyclaw.data import Data
#probdata = Data('setprob.data')

a = 1.
sigma = 0.5
h0 = 0.1
grav = 9.81
omega = sqrt(2.*grav*h0) / a

def maketopo():
    """
    Output topography file for the entire domain
    """
    nxpoints=201
    nypoints=201
    xupper=2.5e0
    yupper=2.5e0
    xlower = -2.5e0
    ylower = -2.5e0
    outfile= "bowl.topotype2"
    topo2writer(outfile,topo,xlower,xupper,ylower,yupper,nxpoints,nypoints)

def maketopo2():
    """
    Another topo file for testing
    """
    nxpoints=101
    nypoints=101
    xupper=1.e0
    yupper=1.e0
    xlower = 0.5e0
    ylower = 0.5e0
    outfile= "bowl2.topotype2"
    topo2writer(outfile,topo2,xlower,xupper,ylower,yupper,nxpoints,nypoints)

def maketopo3():
    """
    Another topo file for testing
    """
    nxpoints=21
    nypoints=101
    xlower = 0.6e0
    xupper= 0.7e0
    ylower = 0.7e0
    yupper=1.2e0
    outfile= "bowl3.topotype2"
    topo2writer(outfile,topo3,xlower,xupper,ylower,yupper,nxpoints,nypoints)

def maketopo4():
    """
    Another topo file for testing
    """
    nxpoints=301
    nypoints=21
    xlower = 0.0e0
    xupper= 1.5e0
    ylower = 0.8e0
    yupper= 0.9e0
    outfile= "bowl4.topotype2"
    topo2writer(outfile,topo4,xlower,xupper,ylower,yupper,nxpoints,nypoints)

def maketopo5():
    """
    Another topo file for testing
    """
    nxpoints=151
    nypoints=251
    xlower = 0.65e0
    xupper= 0.8e0
    ylower = 0.85e0
    yupper= 1.1e0
    outfile= "bowl5.topotype2"
    topo2writer(outfile,topo5,xlower,xupper,ylower,yupper,nxpoints,nypoints)

def maketopo6():
    """
    Another topo file for testing
    """
    nxpoints=31
    nypoints=46 
    xlower = 0.3e0
    xupper= 0.9e0
    ylower = 0.6e0
    yupper= 1.5e0
    outfile= "bowl6.topotype2"
    topo2writer(outfile,topo6,xlower,xupper,ylower,yupper,nxpoints,nypoints)

def topo(x,y):
    """
    Parabolic bowl
    """
    #z = h0*(x**2 + y**2)/a**2 - h0
    z = ones(x.shape)
    return z

def topo2(x,y):
    z = 2*ones(x.shape)
    return z

def topo3(x,y):
    z = 3*ones(x.shape)
    return z

def topo4(x,y):
    z = 4*ones(x.shape)
    return z

def topo5(x,y):
    z = 5*ones(x.shape)
    return z

def topo6(x,y):
    z = 6*ones(x.shape)
    return z


if __name__=='__main__':
    maketopo()
    maketopo2()
    maketopo3()
    maketopo4()
    maketopo5()
    maketopo6()
