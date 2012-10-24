#!/usr/bin/env python

import sys

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.colors as colors 

import clawpack.visclaw.colormaps as colormaps

def read_topo_header(path,topo_type=3):
    r"""Read in header of topography file at path.

    If a value returns numpy.nan then the value was not retrievable.
    """

    # Default values to track errors
    num_cells = [np.nan,np.nan]
    extent = [np.nan,np.nan,np.nan,np.nan]
    delta = np.nan
    no_data_value = np.nan

    bathy_file = open(path,'r')

    if topo_type == 3:
        num_cells[0] = int(bathy_file.readline().split()[0])
        num_cells[1] = int(bathy_file.readline().split()[0])
        extent[0] = float(bathy_file.readline().split()[0])
        extent[2] = float(bathy_file.readline().split()[0])
        delta = float(bathy_file.readline().split()[0])
        no_data_value = float(bathy_file.readline().split()[0])
        
        extent[1] = extent[0] + num_cells[0] * delta
        extent[3] = extent[2] + num_cells[1] * delta
    else:
        raise NotImplemented("Topo type header reading not implemented.")

    bathy_file.close()

    return num_cells,extent,delta,no_data_value


def read_topo(path,topo_type=3):
    r"""Read in topography data

    Depending on the topography type, returns:
     1) 1D arrays x,y,z
     3) 2D arrays X,Y,Z
    """

    if topo_type == 3:
        N,extent,delta,no_data_value = read_topo_header(path)
        x = np.linspace(extent[0],extent[1],N[0])
        y = np.linspace(extent[3],extent[2],N[1])
        X,Y = np.meshgrid(x,y)
        # Data is read in starting at the top right corner
        Z = np.loadtxt(path,skiprows=6)

    else:
        raise NotImplemented('Topo type reading not implemented.')

    return X,Y,Z

def plot(path,coastlines=True,axes=None):
    r"""Plot the bathymetry file at path.

    Returns an axes instance.
    """

    # Create an axes instance if not provided
    if axes is None:
        fig = plt.figure()
        axes = fig.add_subplot(111)

    # Read in bathy
    X,Y,Z = read_topo(path)
    region_extent = (np.min(X),np.max(X),np.min(Y),np.max(Y))
    depth_extent = (np.min(Z),np.max(Z))

    # Create color map
    # cmap = colormaps.make_colormap({-1:[0.3,0.2,0.1],
    #                                    -0.00001:[0.95,0.9,0.7],
    #                                    0.00001:[.5,.7,0],
    #                                    1:[.2,.5,.2]})
    # cmap = plt.get_cmap('terrain')
    # color_norm = colors.Normalize(depth_extent[0],depth_extent[1],clip=True)

    # Plot data
    # plot = axes.imshow(Z,vmin=depth_extent[0], vmax=depth_extent[1],
    #                      extent=region_extent, cmap=cmap, norm=color_norm)
    levels = range(0,int(-np.min(Z)),500)
    axes.contour(X,Y,-Z,levels=levels,colors='gray')

    # Plot coastlines
    if coastlines:
        axes.contour(X,Y,Z,levels=[0.0],colors='r')

    axes.set_xlim(region_extent[0:2])
    axes.set_ylim(region_extent[2:])
    # axes.set_title('Region')

    return axes

if __name__ == '__main__':
    if len(sys.argv) > 1:
        path = sys.argv[1]
    else:
        path = './gulf_coarse_bathy.tt3'
    plot(path)

    plt.show()

