#!/usr/bin/env python
# encoding: utf-8
""" 
Functions for creating topography
    
"""

import os
import sys

import numpy as np
from scipy.interpolate import griddata

import matplotlib.pyplot as plt
import matplotlib.colors as colors 

import clawpack.visclaw.colormaps as colormaps
import bathy

# Degree to meter conversion function
R_earth = 6378.1 * 1000.0
deg2meters = lambda theta,lat:R_earth * theta * np.pi / 180.0 * np.cos(lat * np.pi / 180.0)
meters2deg = lambda d,lat:d / (R_earth * np.pi / 180.0 * np.cos(lat * np.pi / 180.0))

# Read in fort.14 file
def extract_bathy(path, extent, no_data_value=999999, method='nearest', 
                                delta_limit=20.0, TOLERANCE=1e-3):
    r"""Extract sub-section of bathymetry from file at path

    Function to extract a sub-section given by extent of the bathymetry file at
    path assumed to be in a x,y,z format which can be unstructured.  Returns the
    data interpolated onto a grid determined by the resolution of the original
    file or the limiting resolution delta_limit.

    :Input:
     *path* (string) - Path to the bathymetry file which the data is being 
                       pulled from.
     *extent* (tuple) - A tuple defining the rectangle of the sub-section.  Must
                        be in the form (x lower,x upper,y lower, y upper).
     *no_data_value* (float) - Value to use if no data was found to fill in a 
                               missing value, ignored if `method = 'nearest'`.
                               Default is `999999`.
     *method* (string) - Method for interpolation, valid methods are found in
                         the scipy module scipy.interpolate.  Default is 
                         `nearest`.
     *delta_limit* (float) - Limit of finest horizontal resolution, default is
                             20 meters.
     *tolerance* (float) -  Tolerance allowed for extent matching.  Since the 
                            requested extents and the eventual output may not 
                            match due to round off, this parameter is used to 
                            check if they are within acceptable tolerances.
                            Default is `1e-3`.

    :Output:
     *Z* (ndarray) - Interpolated 2D array of bathymetry depths starting in the
                     upper right corner of the sub-section specified by extent.
     *delta* (float) - Final choice used for the horizontal resolution.
    """

    # Extract data
    print "Loading data from file %s" % path
    points = []
    values = []
    dx = np.infty
    dy = np.infty

    data_file = open(path,'r')
    data_file.readline()
    num_lines = int(data_file.readline().split()[1])

    data = np.empty((3,num_lines))

    for n in xrange(num_lines):
        line = data_file.readline().split()
        coordinates = [float(value) for value in line[1:]]
        if extent[0] <= coordinates[0] <= extent[1]:
            if extent[2] <= coordinates[1] <= extent[3]:
                points.append(coordinates[0:2])
                values.append(coordinates[2])

                # Try to determine smallest dx and dy
                if len(points) > 1:
                    if np.abs(points[-1][0] - points[-2][0]) < dx:
                        dx = np.abs(points[-1][0] - points[-2][0])
                    if np.abs(points[-1][1] - points[-2][1]) < dy:
                        dy = np.abs(points[-1][1] - points[-2][1])

    data_file.close()

    if len(points) == 0:
        raise Exception("No points were found inside requested extent.")

    # Cast lists as ndarrays
    points = np.array(points)
    values = np.array(values)    

    # Create regularized grid
    print "Computing grid data"
    delta = max(min(dx,dy),delta_limit) # Limit to size of delta
    N = (np.ceil((extent[1] - extent[0]) / delta),
         np.ceil((extent[3] - extent[2]) / delta))
    print "  delta = %s, N = %s" % (delta,N)
    if N[0] > 2000 or N[1] > 2000:
        raise Exception("Calculated resolution %s too high!")
    x = np.linspace(extent[0],extent[1],N[0])
    y = np.linspace(extent[2],extent[3],N[1])
    X,Y = np.meshgrid(x,y)

    # Check extents
    if abs(x[0]  - extent[0]) > TOLERANCE or \
       abs(x[-1] - extent[1]) > TOLERANCE or \
       abs(y[0]  - extent[2]) > TOLERANCE or \
       abs(y[-1] - extent[3]) > TOLERANCE:

       raise Exception("Calculated grid out of extent tolerance.")

    # Interpolate known points onto regularized grid
    print "Creating interpolating function..."
    Z = griddata(points,values,(X,Y), method=method, fill_value=no_data_value)

    return Z,delta


def write_bathy(path,Z,lower,delta,no_data_value=999999,topotype=3):
    r"""Write out a topography file to path of type topotype

    Writes out a bathymetry file of type 3 to path from data in Z.  The rest of
    the arguments are used to write the header data.
    """

    outfile = open(path,'w')

    if topotype == 3:
        # Write out header
        outfile.write('%s ncols\n' % Z.shape[1])
        outfile.write('%s nrows\n' % Z.shape[0])
        outfile.write('%s xll\n' % lower[0])
        outfile.write('%s yll\n' % lower[1])
        outfile.write('%s cellsize\n' % delta)
        outfile.write('%s nodata_value\n' % no_data_value)

        # Write out bathy data
        # We flip the output data here since we write from the upper left corner
        # to lower right and the data is ordered from lower left to upper right
        Z_flipped = np.flipud(Z)
        for i in xrange(Z.shape[0]):
            for j in xrange(Z.shape[1]):
                outfile.write("%s   " % (Z_flipped[i,j]))
            outfile.write("\n")

    else:
        raise NotImplemented("Output type %s not implemented." % topotype)

    outfile.close()


def plot_bathy(paths,region_path,patch_edges=True,patch_names=True,names=None,
               plot_coastline=True):
    r"""Plot the bathymetry files specified in paths and region_path."""
    
    # Setup region figure
    region_fig = plt.figure(1)
    region_axes = region_fig.add_subplot(111)

    # Setup patch figure
    patch_fig = plt.figure(2)
    columns = 3
    rows = np.ceil(len(paths) / float(columns))
    patch_axes = [patch_fig.add_subplot(rows,columns,i) for i in xrange(len(paths))]

    # Read in region bathymetry
    X,Y,Z = bathy.read_topo(region_path)
    region_extent = (np.min(X),np.max(X),np.min(Y),np.max(Y))
    depth_extent = (np.min(Z),np.max(Z))

    # Create color map
    cmap = colormaps.make_colormap({-1:[0.3,0.2,0.1],
                                       -0.00001:[0.95,0.9,0.7],
                                       0.00001:[.5,.7,0],
                                       1:[.2,.5,.2]})
    color_norm = colors.Normalize(depth_extent[0],depth_extent[1],clip=True)

    # Plot region data
    region_plot = region_axes.imshow(Z,vmin=depth_extent[0],
                                       vmax=depth_extent[1],
                                       extent=region_extent,
                                       cmap=cmap,norm=color_norm)
    
    if plot_coastline:
        region_axes.contour(X,Y,Z,levels=[0.0],colors='r')


    # Read in and plot each patch
    for (i,patch_path) in enumerate(paths):
        X,Y,Z = bathy.read_topo(patch_path)
        extent = (np.min(X),np.max(X),np.min(Y),np.max(Y))

        # Plot on region figure
        region_axes.imshow(Z,vmin=depth_extent[0],vmax=depth_extent[1],
                             extent=extent,cmap=cmap,norm=color_norm)

        # Plot boundaries of local bathy on region plot
        if patch_edges:
            # Bottom boundary
            region_axes.plot((extent[0],extent[1]),(extent[2],extent[2]),'k')
            # Upper boundary
            region_axes.plot((extent[0],extent[1]),(extent[3],extent[3]),'k')
            # Left boundary
            region_axes.plot((extent[0],extent[0]),(extent[2],extent[3]),'k')
            # Right boundary
            region_axes.plot((extent[1],extent[1]),(extent[2],extent[3]),'k')

        # Write name near edge
        if names is None:
            file_name = os.path.splitext(patch_path)[0]
        else:
            file_name = names[i]
        if patch_names:
            delta = X[0,0] - X[1,0]
            region_axes.text(extent[0]+delta,extent[2]+delta,file_name,color='m')

        # Plot on local bathy
        patch_axes[i].imshow(Z,vmin=depth_extent[0],vmax=depth_extent[1],
                             extent=extent,cmap=cmap,norm=color_norm)
        patch_axes[i].set_title(file_name)
        patch_axes[i].set_xlim(extent[0:2])
        patch_axes[i].set_ylim(extent[2:])
        patch_axes[i].set_xlabel('longitude')
        patch_axes[i].set_ylabel('latitude')

    # Fix up figures
    region_axes.set_xlim(region_extent[0:2])
    region_axes.set_ylim(region_extent[2:])
    region_axes.set_title('Region')
    region_fig.colorbar(region_plot)
    
    # patch_fig.colorbar(region_plot)

    return region_fig,patch_fig

        
if __name__ == "__main__":
    bathy_file = './fort.14'
    regions = {'full_bathy':[-2250.0, 2250.0, -1500.0, 1500.0]}
    for region_name,extent in regions.iteritems():
        out_path = "%s.tt3" % region_name
        Z,delta = extract_bathy(bathy_file,extent)
        print np.max(Z),np.min(Z)
        write_bathy(out_path,-Z,(extent[0],extent[2]),delta)

    # Plot results
    plot_bathy(['.'.join((name,'tt3')) for name in regions.keys()],
               'full_bathy.tt3',names=regions.keys(),patch_names=False)
    plt.show()


