
"""
Module to create topo and qinit data files for this example.
"""

import clawpack.geoclaw.topotools as topotools
import numpy

a = 1.
sigma = 0.5
h0 = 0.1
grav = 9.81
omega = numpy.sqrt(2.0 * grav * h0) / a

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
    outfile= "bowl.nc"

    parabolic_bowl = lambda x,y: h0 * (x**2 + y**2) / a**2 - h0

    topo = topotools.Topography(topo_func=parabolic_bowl)
    topo.x = numpy.linspace(xlower, xupper, nxpoints)
    topo.y = numpy.linspace(ylower, yupper, nypoints)
    topo.write(outfile)


if __name__=='__main__':
    maketopo()
