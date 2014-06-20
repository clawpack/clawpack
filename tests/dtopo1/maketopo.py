
"""
Module to create topo and qinit data files for this example.
"""

from clawpack.geoclaw.topotools import topo1writer, topo2writer
from clawpack.geoclaw import dtopotools
from numpy import *
import os


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
    dx = (xupper-xlower)/(nxpoints-1)
    dy = (yupper-ylower)/(nypoints-1)
    print "==> topo1 has dx = %g, dy = %g" % (dx,dy)

def maketopo2():
    """
    Another topo file for testing
    """
    nxpoints=21
    nypoints=51
    xlower = -0.5
    xupper= -0.3
    ylower = -0.1
    yupper= 0.4
    outfile= "topo2.topotype2"
    topo2writer(outfile,topo2,xlower,xupper,ylower,yupper,nxpoints,nypoints)
    dx = (xupper-xlower)/(nxpoints-1)
    dy = (yupper-ylower)/(nypoints-1)
    print "==> topo2 has dx = %g, dy = %g" % (dx,dy)

h0 = 1000.

def topo(x,y):
    #z = -10*ones(x.shape)
    z = -h0*(1 + 0.5*cos(x-y))
    return z

def topo2(x,y):
    #z = -15*ones(x.shape)
    z = -h0*(1. + exp(x+y))
    return z


def read_subfaults(fname_subfaults, plotfig=None):
    """
    Test data
    """

    # Format of subfault file:
    columns = """longitude latitude depth length width strike dip rake slip
       """.split()
    defaults = {'latlong_location': 'top center'}
    units = {'slip': 'm', 'depth': 'km', 'length': 'km', 'width': 'km'}

    subfaults = dtopotools.read_subfault_model(fname_subfaults, \
                        columns=columns, units=units, \
                        defaults = defaults, skiprows=1, delimiter=',')
                    
    if plotfig:
        figure(plotfig)
        dtopotools.plot_subfaults(subfaults,slip_color=True,plot_rake=True)
        
    return subfaults


def make_dtopo1(plotfig=None):
    """
    Test data.
    """

    fname_subfaults = 'dtopo1.csv'
    subfaults = read_subfaults(fname_subfaults)

    dtopo_fname = fname_subfaults.split('.')[0] + '.tt3'

    #if os.path.exists(dtopo_fname):
    #    print "*** Not regenerating dtopo file (already exists): %s" \
    #            % dtopo_fname
    #else:
    if 1:
        print "Using Okada model to create %s " % dtopo_fname

        # Needed for extent of dtopo file:
        xlower = -0.4
        xupper = 0.6
        ylower = -0.4
        yupper = 0.4

        # number of grid points in dtopo file
        mx = 151
        my = 121

        dx = (xupper-xlower)/(mx-1)
        dy = (yupper-ylower)/(my-1)
        print "==> dtopo1 has dx = %g, dy = %g" % (dx,dy)

        # Create dtopo_params dictionary with parameters for dtopo file: 
        dtopo_params = {}
        dtopo_params['fname'] = dtopo_fname
        dtopo_params['faulttype'] = 'static'
        dtopo_params['dtopotype'] = 3
        dtopo_params['mx'] = mx
        dtopo_params['my'] = my
        dtopo_params['xlower'] = xlower
        dtopo_params['xupper'] = xupper
        dtopo_params['ylower'] = ylower
        dtopo_params['yupper'] = yupper
        dtopo_params['t0'] = 0.
        dtopo_params['tfinal'] = 1.
        dtopo_params['ntimes'] = 25

        dtopo = dtopotools.make_dtopo_from_subfaults(subfaults, dtopo_params)

        if plotfig:
            figure(plotfig)
            x = dtopo.x
            y = dtopo.y
            dz_final = dtopo.dz_list[-1]
            dtopotools.plot_dz_colors(x,y,dz_final,cmax_dz=16,dz_interval=1)

        return dtopo


def make_dtopo2(plotfig=None):
    """
    Test data.
    """

    fname_subfaults = 'dtopo2.csv'
    subfaults = read_subfaults(fname_subfaults)

    dtopo_fname = fname_subfaults.split('.')[0] + '.tt3'

    #if os.path.exists(dtopo_fname):
    #    print "*** Not regenerating dtopo file (already exists): %s" \
    #            % dtopo_fname
    #else:
    if 1:
        print "Using Okada model to create %s " % dtopo_fname

        # Needed for extent of dtopo file:
        xlower = -0.9
        xupper = 0.1
        ylower = -0.4
        yupper = 0.4

        # number of grid points in dtopo file
        mx = 201
        my = 161

        dx = (xupper-xlower)/(mx-1)
        dy = (yupper-ylower)/(my-1)
        print "==> dtopo2 has dx = %g, dy = %g" % (dx,dy)

        # Create dtopo_params dictionary with parameters for dtopo file: 
        dtopo_params = {}
        dtopo_params['fname'] = dtopo_fname
        dtopo_params['faulttype'] = 'static'
        dtopo_params['dtopotype'] = 3
        dtopo_params['mx'] = mx
        dtopo_params['my'] = my
        dtopo_params['xlower'] = xlower
        dtopo_params['xupper'] = xupper
        dtopo_params['ylower'] = ylower
        dtopo_params['yupper'] = yupper
        dtopo_params['t0'] = 0.5
        dtopo_params['tfinal'] = 1.2
        dtopo_params['ntimes'] = 25

        dtopo = dtopotools.make_dtopo_from_subfaults(subfaults, dtopo_params)

        if plotfig:
            figure(plotfig)
            x = dtopo.x
            y = dtopo.y
            dz_final = dtopo.dz_list[-1]
            dtopotools.plot_dz_colors(x,y,dz_final,cmax_dz=16,dz_interval=1)

        return dtopo

if __name__=='__main__':
    maketopo()
    maketopo2()
    make_dtopo1()
    make_dtopo2()
