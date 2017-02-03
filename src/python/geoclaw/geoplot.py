"""
Useful things for plotting GeoClaw results.
"""

from __future__ import absolute_import
from clawpack.visclaw import colormaps
from matplotlib.colors import Normalize 
from clawpack.geoclaw import topotools
from numpy import ma
from six.moves import zip


# Colormaps from geoclaw
# Color attributes, single instance per run
# Colors
black = [0.0,0.0,0.0]
white = [1.0,1.0,1.0]
red = [1.0,0.0,0.0]
green = [0.0,1.0,0.0];
dark_green = [0.1,0.4,0.0];
light_green = [0.8,1.0,0.5];
blue = [0.0,0.0,1.0];
dark_blue = [0.2,0.2,0.7];
light_blue = [0.5,0.5,1.0];
blue_green = [0.0,1.0,1.0];
tan = [0.9,0.8,0.2];
tan = [0.8,0.5,0.2];
brown = [0.9,0.8,0.2];
gray8 = [0.8,0.8,0.8];
purple = [0.8,0.3,0.8];

# Colormaps
TSUNAMI_MAX_AMPLITUDE = 0.6
tsunami_colormap = colormaps.make_colormap({-TSUNAMI_MAX_AMPLITUDE:blue,
                                            0.0:blue_green,
                                            TSUNAMI_MAX_AMPLITUDE:red})
                                            
land1_colormap = colormaps.make_colormap({0.0:dark_green,
                                          1000.0:green,
                                          2000.0:light_green,
                                          4000.0:tan})
                                         
land2_colormap = colormaps.make_colormap({0:dark_green,
                                          50:green,
                                          100:light_green,
                                          200:tan})
                                          
water_land_colormap = colormaps.make_colormap({-1000:dark_blue,
                                               -500:blue,
                                               0:light_blue,
                                               .1:tan,
                                               5:tan,
                                               6:dark_green,
                                               1000:green,
                                               2000:light_green,
                                               4000:tan})
                                               
bathy1_colormap = colormaps.make_colormap({-1000:brown,
                                           0:tan,
                                           .1:dark_green,
                                           1000:green,
                                           2000:light_green})
       
bathy2_colormap = colormaps.make_colormap({-1000:brown,
                                           -100:tan,
                                           0:dark_green,
                                           .1:dark_green,
                                           1000:green,
                                           2000:light_green})

bathy3_colormap = colormaps.make_colormap({-1:[0.3,0.2,0.1],
                                           -0.01:[0.95,0.9,0.7],
                                           .01:[.5,.7,0],
                                           1:[.2,.5,.2]})

seafloor_colormap = colormaps.make_colormap({-1:[0.3,0.2,0.1],
                                              0:[0.95,0.9,0.7]})

land_colormap = colormaps.make_colormap({ 0:[0.95,0.9,0.7],
                                          1:[.2,.5,.2]})


                                           
colormaps_list = {"tsunami":tsunami_colormap,"land1":land1_colormap,
             "land2":land2_colormap,"water_land":water_land_colormap,
             "bathy1":bathy1_colormap,"bathy2":bathy2_colormap}
        
def plot_colormaps():
    r"""Plots all colormaps avaiable or the ones specified"""
        
    import numpy as np
    import matplotlib.pyplot as plt
    
    a = np.linspace(0, 1, 256).reshape(1,-1)
    a = np.vstack((a,a))
    
    nmaps = len(colormaps_list) + 1

    fig = plt.figure(figsize=(5,10))
    fig.subplots_adjust(top=0.99, bottom=0.01, left=0.2, right=0.99)
    
    for i,name in enumerate(colormaps_list):
        ax = plt.subplot(nmaps,1,i+1)
        plt.axis("off")
        plt.imshow(a, aspect='auto', cmap=colormaps_list[name], origin='lower')
        pos = list(ax.get_position().bounds)
        fig.text(pos[0] - 0.01, pos[1], name, fontsize=10, horizontalalignment='right')

    #plt.show()

land_colors = colormaps.make_colormap({0:[.5,.7,0], 1:[.2,.5,.2]})
# water_colors = colormaps.make_colormap({-1.:'r', 0.:[0, .8, .8], 1.: 'b'})
# land_colors = land2_colormap
water_colors = tsunami_colormap

# Plotting functions

# The drytol parameter is used in masking land and water and
# affects what color map is used for cells with small water depth h.
# The best value to use often depends on the application and can
# be set for an application by setting current_data.user.drytol in
# a beforeframe function, for example.  If it's not set by the user,
# the following default value is used (in meters):

drytol_default = 1.e-3


def topo(current_data):
   """ 
   Return topography = eta - h. 
   Surface eta is assumed to be output as 4th column of fort.q files.
   """
   q = current_data.q
   h = q[0,:,:]
   eta = q[3,:,:]
   topo = eta - h
   return topo

def land(current_data):
   """
   Return a masked array containing the surface elevation only in dry cells.
   """
   from numpy import ma
   drytol = getattr(current_data.user, 'drytol', drytol_default)
   q = current_data.q
   h = q[0,:,:]
   eta = q[3,:,:]
   land = ma.masked_where(h>drytol, eta)
   return land

def water(current_data):
   """Deprecated: use surface instead."""
   raise DeprecationWarning("Deprecated function, use surface instead.")
   from numpy import ma
   drytol = getattr(current_data.user, 'drytol', drytol_default)
   q = current_data.q
   h = q[0,:,:]
   eta = q[3,:,:]
   water = ma.masked_where(h<=drytol, eta)
   return water

def depth(current_data):
   """
   Return a masked array containing the depth of fluid only in wet cells.
   """
   from numpy import ma
   drytol = getattr(current_data.user, 'drytol', drytol_default)
   q = current_data.q
   h = q[0,:,:]
   depth = ma.masked_where(h<=drytol, h)
   return depth

def surface(current_data):
    """
    Return a masked array containing the surface elevation only in wet cells.
    Surface is eta = h+topo, assumed to be output as 4th column of fort.q
    files.
    """
    from numpy import ma
    drytol = getattr(current_data.user, 'drytol', drytol_default)
    q = current_data.q
    h = q[0,:,:]
    eta = q[3,:,:]
    water = ma.masked_where(h<=drytol, eta)
    return water

def surface_or_depth(current_data):
    """
    Return a masked array containing the surface elevation where the topo is 
    below sea level or the water depth where the topo is above sea level.
    Mask out dry cells.  Assumes sea level is at topo=0.
    Surface is eta = h+topo, assumed to be output as 4th column of fort.q
    files.
    """
    from numpy import ma, where
    drytol = getattr(current_data.user, 'drytol', drytol_default)
    q = current_data.q
    h = q[0,:,:]
    eta = q[3,:,:]
    topo = eta - h
    surface = ma.masked_where(h<=drytol, eta)
    depth = ma.masked_where(h<=drytol, h)
    surface_or_depth = where(topo<0, surface, depth)
    return surface_or_depth

# Some discrete color maps useful for contourf plots of fgmax results:

def discrete_cmap_1(clines):
    """
    Construct a discrete color map for the regions between the contour lines
    given in clines. Colors go from turqouise through yellow to red.
    Good for zeta.
    """
    from numpy import floor, linspace, hstack, ones, zeros
    nlines = len(clines)
    n1 = int(floor((nlines-1)/2.))
    n2 = nlines - 1 - n1
    Green = hstack([linspace(1,1,n1),linspace(1,0,n2)])
    Red = hstack([linspace(0,0.8,n1), ones(n2)])
    Blue = hstack([linspace(1,0.2,n1), zeros(n2)])
    colors = list(zip(Red,Green,Blue))
    return colors

def discrete_cmap_2(clines):
    """
    Construct a discrete color map for the regions between the contour lines
    given in clines.  Colors go from red to turquoise.
    Good for arrival times.
    """
    from numpy import floor, linspace, hstack, ones, zeros, flipud
    nlines = len(clines)
    n1 = int(floor((nlines-1)/2.))
    n2 = nlines - 1 - n1
    Green = flipud(hstack([linspace(1,1,n1),linspace(1,0,n2)]))
    Red = flipud(hstack([linspace(0,0.8,n1), ones(n2)]))
    Blue = flipud(hstack([linspace(1,0.2,n1), zeros(n2)]))
    colors = list(zip(Red,Green,Blue))
    return colors

