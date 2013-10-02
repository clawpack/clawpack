
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
    nxpoints=200
    nypoints=200
    xupper=2.e0
    yupper=2.e0
    xlower = -2.e0
    ylower = -2.e0
    outfile= "bowl.topotype2"
    topo2writer(outfile,topo,xlower,xupper,ylower,yupper,nxpoints,nypoints)

def topo(x,y):
    """
    Parabolic bowl
    """
    z = h0*(x**2 + y**2)/a**2 - h0
    return z


if __name__=='__main__':
    maketopo()
