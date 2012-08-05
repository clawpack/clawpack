"""
Create initial perturbation file
"""

import os
import sys

import numpy as np
from numpy import ma
import scipy.interpolate as interp
import matplotlib as mpl
import matplotlib.pyplot as plt

from geoclaw import topotools
from clawpack.visclaw import geoplot
from clawpack.visclaw import colormaps

def plot_bathy(bathy_file,coarse_factor=5,topo_type=3):
    r"""docstring for plot_bathy"""

    # Geometry
    num_cells = [0,0]
    x_range = [0.0,0.0]
    y_range = [0.0,0.0]
    delta = 0.0
    no_data_value = 0
    
    if abs(topo_type) == 1:
        # Load in data
        file_data = np.loadtxt(bathy_file)
        x = file_data[:,0]
        y = file_data[:,1]
        bathy = file_data[:,2]
        
        # Figure out extents of domain
        lower_corner = [np.min(x),np.min(y)]
    elif abs(topo_type) == 3:
        # Parse header
        bathy_file_handle = open(bathy_file,'r')
        num_cells[0] = int(bathy_file_handle.readline().split()[0])
        num_cells[1] = int(bathy_file_handle.readline().split()[0])
        x_range[0] = float(bathy_file_handle.readline().split()[0])
        y_range[0] = float(bathy_file_handle.readline().split()[0])
        delta = float(bathy_file_handle.readline().split()[0])
        no_data_value = int(bathy_file_handle.readline().split()[0])
            
        bathy_file_handle.close()
            
        # Calculate additional geometry
        x_range[1] = x_range[0] + num_cells[0] * delta
        y_range[1] = y_range[0] + num_cells[1] * delta
            
        # Load in bathymetry
        bathy = np.flipud(np.loadtxt(bathy_file,skiprows=6))
        X,Y = np.meshgrid(np.linspace(x_range[0],x_range[1],num_cells[0]),
                          np.linspace(y_range[0],y_range[1],num_cells[1]))

    else:
        raise NotImplemented("Topography type %s has not been implemented." % topo_type)
    if topo_type < 0:
        bathy = -bathy
    
    # Plot bathymetry
    fig = plt.figure()
    axes = fig.add_subplot(111)
    nan_array = np.nan * np.empty(bathy.shape)
    land = ma.masked_where(bathy > 0.0,nan_array)
    sea_floor = ma.masked_where(bathy <= 0.0,nan_array)
    max_elevation = np.max(bathy)
    max_depth = np.min(bathy)
    land_cmap = colormaps.make_colormap({0:[.5,.7,0],max_elevation:[.2,.5,.2]})
    sea_cmap = colormaps.make_colormap({max_depth:[0.3,0.2,0.1],0:[0.95,0.9,0.7]})
    
    plot = axes.pcolor(X[0:-1:coarse_factor,0:-1:coarse_factor],
                       Y[0:-1:coarse_factor,0:-1:coarse_factor],
                       land[0:-1:coarse_factor,0:-1:coarse_factor],
                       cmap=land_cmap,vmin=0.0,vmax=max_elevation)
    plot = axes.pcolor(X[0:-1:coarse_factor,0:-1:coarse_factor],
                       Y[0:-1:coarse_factor,0:-1:coarse_factor],
                       sea_floor[0:-1:coarse_factor,0:-1:coarse_factor],
                       cmap=sea_cmap,vmin=max_depth,vmax=0.0)
        
    axes.set_xlim(x_range)
    axes.set_ylim(y_range)
    axes.set_xlabel('Latitude')
    axes.set_ylabel('Longitude')
    
    fig.colorbar(plot)
    
    plt.show()
    
def plot_bathy_2():
    from clawpack.visclaw import geoplot

    T = geoplot.TopoPlotData('bathy/gulf_coarse_bathy.tt3')
    T.topotype = 3
    T.cmin = -20
    T.cmax = 50.0
    T.plot()

    T.topotype = 3
    T.fname = 'bathy/galveston_bay.tt3'
    T.plot()

    plt.xlim(-99.0,-80.0)
    plt.ylim(17.0,32.0)
    
    plt.colorbar()
    
    plt.show()

    import pdb; pdb.set_trace()


def convert_adcirc_bathy(topo_file,output_file="converted.xyz",topo_type=1):
    
    raise NotImplemented("This function does not work yet.")
    
    # Read header
    topo_file_handle = open(topo_file,'r')
    data_name = topo_file_handle.readline().strip()
    print "Data name = %s" % data_name
    data_size = topo_file_handle.readline().split()
    num_cells = int(data_size[0])
    num_nodes = int(data_size[1])
    print "Number of cells = %s" % num_cells
    print "Number of nodes = %s" % num_nodes
    
    # Extract bathymetry tuples
    print "Reading in data from %s..." % topo_file
    topo_data = np.genfromtxt(topo_file,usecols=(1,2,3),skip_header=2,
                                                        skip_footer=num_nodes)
    print "Done reading."
    topo_data[:,2] = -topo_data[:,2]
    
    # Convert scatter data to gridded data
    print "Converting scattered data to gridded data..."
    lower = [-99.0,17.0]
    upper = [-80.0,24.0]

    refinement_factor = int(np.exp(0.5 * np.log(float(num_cells) / 
                                                float(abs(upper[0] - lower[0]) * 
                                                      abs(upper[1]-lower[1])) )))
    N = (int(abs(upper[0] - lower[0]) * refinement_factor),
         int(abs(upper[1] - lower[1]) * refinement_factor))
    print "Calculated refinement factor = %s" % refinement_factor
    grid_x, grid_y = np.mgrid[lower[0]:upper[0]:200j, lower[1]:upper[1]:100j]
    interp.griddata(topo_data[:,0:2],topo_data[:,2],(grid_x,grid_y),method='linear')
    print "Done converting."
    
    # Write out new file
    print "Writing out new bathymetry data."
    if topo_type == 1:
        np.savetxt(output_file,topo_data)
    else:
        raise Exception("Unsupported topography type %s." % topo_type)
    print "Wrote out data file to %s." % output_file
        
        
def make_qinit(outfile="hump.xyz"):
    """
    Create qinit data file
    """
    nxpoints = 200
    nypoints = 200
    xlower = -94.0
    xupper = -92.0
    ylower = 27.0
    yupper = 29.0
    topotools.topo1writer(outfile,qinit,xlower,xupper,ylower,yupper,nxpoints,nypoints)


def qinit(x,y):
    """
    Gaussian hump:
    """
    x0 = -93.0
    y0 = 28.0
    d = topotools.gcdist(x,y,x0,y0)  # great circle distance on earth
    ze = -d**2 / 2.e9
    amplitude = 1.0
    z = np.where(ze > -100., amplitude * np.exp(ze), 0.0)
    return z

def plot_qinit(qinit_file="hump.xyz"):
    r"""Plot the given qinit file"""
    
    data = np.loadtxt(qinit_file)
    x = data[:,0].reshape((100,100))
    y = data[:,1].reshape((100,100))
    perturbation = data[:,2].reshape((100,100))
    
    fig = plt.figure()
    axes = fig.add_subplot(111)
    
    hump_plot = axes.pcolor(x,y,perturbation)
    fig.colorbar(hump_plot)
    
    plt.show()

if __name__=='__main__':
    # convert_adcirc_bathy('gulf_bathy.xyz')
    # convert_adcirc_bathy('test_bathy.xyz')
    # plot_bathy_2()
    make_qinit()
    plot_qinit()
    # plot_bathy("./gulf_coarse_bathy.tt3",topo_type=-3)
    