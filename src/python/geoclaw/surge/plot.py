#!/usr/bin/env python
# encoding: utf-8
r"""
Plotting routines for storm surge simulations with GeoClaw

:Authors:
    Kyle Mandli (2012-10-23) Initial version
"""
# ==============================================================================
#  Distributed under the terms of the Berkeley Software Distribution (BSD) 
#  license
#                     http://www.opensource.org/licenses/
# ==============================================================================

from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

import clawpack.visclaw.colormaps as colormaps
import clawpack.visclaw.gaugetools as gaugetools

import clawpack.geoclaw.geoplot as geoplot

# Default location of storm fields in aux array output
# Note that selective output can move these around, when importing the module
# make sure to change these values before plotting to use
bathy_index = 0
friction_field = 3
wind_field = 4
pressure_field = 6
# bathy_index = 0
# friction_field = 1
# wind_field = 2
# pressure_field = 4

class figure_counter(object):

    def __init__(self,initial_value=1):
        self._counter = 1

    def get_counter(self):
        self._counter += 1
        return self._counter - 1


class track_data(object):

    def __init__(self,path):
        try:
            self._path = path
            self._data = np.loadtxt(self._path)
        except:
            self._data = None

    def get_track(self,frame):
        # If data was not load successfully return None
        if self._data is None or len(self._data.shape) < 2:
            return None, None, None
        # If it appears that our data is not long enough, try reloading file
        # print self._data.shape,frame
        if self._data.shape[0] < frame + 1:
            self._data = np.loadtxt(self._path)
            # print "new data",self._data.shape

            # Check to make sure that this fixed the problem
            if self._data.shape[0] < frame + 1:
                print(" *** WARNING *** Could not find track data for frame %s." % frame)
                return None, None, None
                # raise Exception("Could not find data for frame %s." % frame)

        return self._data[frame,1:]

# ==========================================================================
# Gauge functions
# ==========================================================================
def gauge_locations(current_data,gaugenos='all'):
    plt.hold(True)
    gaugetools.plot_gauge_locations(current_data.plotdata, \
                                    gaugenos=gaugenos, format_string='kx', 
                                    add_labels=True, xoffset=0.02, yoffset=0.02)
    plt.hold(False)

def gaugetopo(current_data):
    q = current_data.q
    h = q[0,:]
    eta = q[3,:]
    topo = eta - h
    return topo
    
def gauge_afteraxes(current_data):
    # Add sea level line
    add_zeroline(current_data)

    # Change time to days
    plt.xlabel('t (days)')
    plt.ylabel('m')
    locs,labels = plt.xticks()
    labels = [r"$%s$" % str(np.trunc(value/(3600.0 * 24))) for value in locs]
    # locs = np.linspace(-12.0,40,52)
    # labels = range(-12,41)
    plt.xticks(locs,labels)

def addgauges(current_data):
    gaugetools.plot_gauge_locations(current_data.plotdata, \
         gaugenos='all', format_string='ko', add_labels=True)

def add_zeroline(current_data):
    # from pylab import plot, legend, xticks, floor
    t = current_data.t
    #legend(('surface','topography'),loc='lower left')
    plt.plot(t, 0*t, 'k')
    # n = int(floor(t.max()/3600.) + 2)
    # xticks([3600*i for i in range(n)])


# ==========================================================================
#  Generic helper functions
# ==========================================================================
def bathy_ref_lines(current_data):
    plt.hold(True)
    y = [amrdata.ylower,amrdata.yupper]
    for ref_line in ref_lines:
        plt.plot([ref_line,ref_line],y,'y--')
    plt.hold(False)


# ========================================================================
#  Surge related helper functions
# ========================================================================
# Surge eye location
def eye_location(cd,track):
    return track.get_track(cd.frameno)
    
def days_figure_title(current_data, land_fall=0.0):
    t = (current_data.t - land_fall) / (60**2 * 24) 
    days = int(t)
    hours = (t - int(t)) * 24.0

    title = current_data.plotaxes.title
    plt.title('%s at day %3i, hour %2.1f' % (title,days,hours))

def surge_afteraxes(current_data, track, land_fall=0.0, plot_direction=False):
    x,y,theta = eye_location(current_data,track)
    if x is not None and y is not None:
        plt.hold(True)
        plt.plot(x,y,'rD',markersize=2)
        if plot_direction:
            plt.quiver(x, y, np.cos(theta), np.sin(theta))
        plt.hold(False)
    days_figure_title(current_data,land_fall)

def friction(cd):
    return cd.aux[friction_field,:,:]
    
def storm_wind(current_data):
    if current_data.level == 1:
        t = current_data.t
        u = wind_x(cd)
        v = wind_y(cd)
        plt.hold(True)
        Q = plt.quiver(current_data.x[::3,::3],current_data.y[::3,::3],
                    u[::3,::3],v[::3,::3])
        # plt.quiverkey(Q,0.5,0.5,50,r'$50 \frac{m}{s}$',labelpos='W',
        #                 fontproperties={'weight':'bold'})
        plt.hold(False)

def wind_x(cd):
    return cd.aux[wind_field,:,:]
def wind_y(cd):
    return cd.aux[wind_field+1,:,:]
def wind_speed(cd):
    return np.sqrt(wind_x(cd)**2 + wind_y(cd)**2)

def pressure(cd):
    # The division by 100.0 is to convert from Pa to millibars
    return cd.aux[pressure_field,:,:] / 100.0

def pressure_gradient_x(cd):
    return cd.aux[pressure_field+1,:,:]

def pressure_gradient_y(cd):
    return cd.aux[pressure_field+2,:,:]
    
def wind_contours(current_data):
    plt.hold(True)
    w = wind_speed(current_data)
    max_w = np.max(np.max(w))
    levels = [0.0,0.25*max_w,0.5*max_w,0.75*max_w,max_w*0.999]
    C = plt.contour(current_data.x,current_data.y,w,levels)
    plt.clabel(C,inline=1)
    plt.hold(False)


# ========================================================================
#  Water helper functions
# ========================================================================
def b(cd):
    return cd.aux[bathy_index,:,:]
    
def extract_eta(h,eta,DRY_TOL=10**-3):
    index = np.nonzero((np.abs(h) < DRY_TOL) + (h == np.nan))
    eta[index[0],index[1]] = np.nan
    return eta

def extract_velocity(h,hu,DRY_TOL=10**-8):
    u = np.zeros(hu.shape)
    index = np.nonzero((np.abs(h) > DRY_TOL) * (h != np.nan))
    u[index[0],index[1]] = hu[index[0],index[1]] / h[index[0],index[1]]
    return u

def eta(cd):
    return extract_eta(cd.q[0,:,:],cd.q[3,:,:])
    
def water_u(cd):
    return extract_velocity(cd.q[0,:,:],cd.q[1,:,:])
    
def water_v(cd):
    return extract_velocity(cd.q[0,:,:],cd.q[2,:,:])
    
def water_speed(current_data):
    u = water_u(current_data)
    v = water_v(current_data)
        
    return np.sqrt(u**2+v**2)
    
def water_quiver(current_data):
    u = water_u(current_data)
    v = water_v(current_data)
        
    plt.hold(True)
    Q = plt.quiver(current_data.x[::2,::2],current_data.y[::2,::2],
                    u[::2,::2],v[::2,::2])
    max_speed = np.max(np.sqrt(u**2+v**2))
    label = r"%s m/s" % str(np.ceil(0.5*max_speed))
    plt.quiverkey(Q,0.15,0.95,0.5*max_speed,label,labelpos='W')
    plt.hold(False)


# ========================================================================
#  Profile functions
# ========================================================================
class PlotProfile(object):

    def __init__(self,slice_value = 0.0):
        self.slice_value = slice_value

    def slice_index(self,cd):
        if cd.grid.y.lower < self.slice_value < cd.grid.y.upper:
            return int((self.slice_value - cd.grid.y.lower) / cd.dy - 0.5)
        else:
            return None

    def bathy_profile(self,current_data):
        index = self.slice_index(current_data)
        if index:
            return current_data.x[:,index], b(current_data)[:,index]
        else:
            return None, None
    
    def surface_profile(self,current_data):
        index = self.slice_index(current_data)
        if index:
            return current_data.x[:,index], eta(current_data)[:,index]
        else:
            return None, None


# ========================================================================
#  Plot items
# ========================================================================
def add_surface_elevation(plotaxes, plot_type='pcolor', bounds=None, 
                                    contours=None, shrink=1.0):
    if plot_type == 'pcolor' or plot_type == 'imshow':            
        plotitem = plotaxes.new_plotitem(name='surface',plot_type='2d_pcolor')
        plotitem.plot_var = geoplot.surface_or_depth

        if bounds is not None:
            if bounds[0] == 0.0:
                plotitem.pcolor_cmap = plt.get_cmap('OrRd')
            else:
                plotitem.pcolor_cmap = \
                              colormaps.make_colormap({1.0:'r',0.5:'w',0.0:'b'})
            plotitem.pcolor_cmin = bounds[0]
            plotitem.pcolor_cmax = bounds[1]
        plotitem.add_colorbar = True
        plotitem.colorbar_shrink = shrink
        plotitem.colorbar_label = "Surface Height (m)"
        plotitem.amr_celledges_show = [0,0,0,0,0,0,0]
        plotitem.amr_patchedges_show = [1,1,1,0,0,0,0]

    elif plot_type == 'contour':
        plotitem = plotaxes.new_plotitem(name='surface', plot_type='2d_contour')
        if bounds is None:
            plotitem.contour_levels = [-2.5,-1.5,-0.5,0.5,1.5,2.5]

        plotitem.plot_var = geoplot.surface_or_depth
        # plotitem.contour_nlevels = 21
        # plotitem.contour_min = -2.0
        # plotitem.contour_max = 2.0
        # plotitem.kwargs = {''}
        plotitem.amr_contour_show = [1,1,1,1,1,1,1]
        plotitem.amr_celledges_show = [0,0,0,0,0,0,0]
        plotitem.amr_patchedges_show = [1,1,1,1,0,0,0]
        plotitem.amr_contour_colors = 'k'
        # plotitem.amr_contour_colors = ['r','k','b']  # color on each level
        # plotitem.amr_grid_bgcolor = ['#ffeeee', '#eeeeff', '#eeffee']

    elif plot_type == 'contourf':
        plotitem = plotaxes.new_plotitem(name='surface', plot_type='2d_contourf')
        plotitem.plot_var = geoplot.surface_or_depth
        if bounds is not None:
            contours = numpy.linspace(bounds[0],bounds[1],11)
            plotitem.contour_levels = contours
            plotitem.fill_cmin = bounds[0]
            plotitem.fill_cmax = bounds[1]
        elif contours is not None:
            plotitem.contour_levels = contours
            plotitem.fill_cmin = min(contours)
            plotitem.fill_cmax = max(contours)

        plotitem.add_colorbar = True
        plotitem.fill_cmap = geoplot.tsunami_colormap
        plotitem.colorbar_shrink = shrink
        plotitem.colorbar_label = "Surface Height (m)"
        plotitem.fill_cmap = plt.get_cmap('OrRd')
        if any((value < 0 for value in plotitem.contour_levels)):
            plotitem.fill_cmap = \
                            colormaps.make_colormap({1.0:'r',0.5:'w',0.0:'b'})

        plotitem.amr_contour_show = [1,1,1,1,1,1,1]
        plotitem.amr_celledges_show = [0,0,0,0,0,0,0]
        plotitem.amr_patchedges_show = [1,1,1,1,0,0,0]
        plotitem.amr_contour_colors = 'k'


def add_speed(plotaxes, plot_type='pcolor', bounds=None,  contours=None,  
                        shrink=1.0):
    if plot_type == 'pcolor' or plot_type == 'imshow':
        plotitem = plotaxes.new_plotitem(name='speed',plot_type='2d_pcolor')
        plotitem.plot_var = water_speed
        # plotitem.plot_var = 1
        plotitem.pcolor_cmap = plt.get_cmap('PuBu')
        if bounds is not None:
            plotitem.pcolor_cmin = bounds[0]
            plotitem.pcolor_cmax = bounds[1]
        plotitem.add_colorbar = True
        plotitem.colorbar_shrink = shrink
        plotitem.colorbar_label = "Current (m/s)"
        plotitem.amr_celledges_show = [0,0,0,0,0,0,0]
        plotitem.amr_patchedges_show = [1,1,1,0,0,0,0]

    elif plot_type == 'quiver':
        plotitem = plotaxes.new_plotitem(plot_type='2d_quiver')
        plotitem.quiver_var_x = water_u
        plotitem.quiver_var_y = water_v
        plotitem.amr_quiver_show = [4,10,10]
        plotitem.amr_show_key = [True,True,False]
        plotitem.key_units = 'm/s'
        
    elif plot_type == 'contour':
        plotitem = plotaxes.new_plotitem(name='speed', plot_type='2d_contour')
        if bounds is None:
            plotitem.contour_levels = [0.5,1.5,3,4.5,6.0]
        plotitem.kwargs = {'linewidths':1}
        
        plotitem.plot_var = water_speed
        plotitem.amr_contour_show = [1,1,1,1,1,1,1]
        plotitem.amr_celledges_show = [0,0,0]
        plotitem.amr_patchedges_show = [1,1,1,1,1,0,0]
        plotitem.amr_contour_colors = 'k'

    elif plot_type == 'contourf':

        plotitem = plotaxes.new_plotitem(name='speed', plot_type='2d_contourf')

        plotitem.add_colorbar = True
        plotitem.colorbar_label = "Current (m/s)"
        plotitem.colorbar_shrink = shrink
        plotitem.fill_cmap = plt.get_cmap('PuBu')
        if bounds is not None:
            plotitem.contour_levels = numpy.linspace(bounds[0],bounds[1],11)
            plotitem.fill_cmin = bounds[0]
            plotitem.fill_cmap = bounds[1]
        elif contours is not None:
            plotitem.contour_levels = contours
            plotitem.fill_cmin = min(contours)
            plotitem.fill_cmax = max(contours)

        # Modify the 'extends' plot attribute as we don't want this to extend
        # below 0
        plotitem.kwargs['extend'] = 'max'
        
        plotitem.plot_var = water_speed
        plotitem.amr_contour_show = [1,1,1,1,1,1,1]
        plotitem.amr_celledges_show = [0,0,0]
        plotitem.amr_patchedges_show = [1,1,1,1,1,0,0]
        plotitem.amr_contour_colors = 'k'


def add_friction(plotaxes,bounds=None,plot_type='pcolor',shrink=1.0):
    if plot_type == 'pcolor' or plot_type == 'imshow':
        plotitem = plotaxes.new_plotitem(name='friction',plot_type='2d_pcolor')
        plotitem.plot_var = friction
        plotitem.pcolor_cmap = plt.get_cmap('YlOrRd')
        plotitem.colorbar_shrink = shrink
        if bounds is not None:
            plotitem.pcolor_cmin = bounds[0]
            plotitem.pcolor_cmax = bounds[1]
        plotitem.add_colorbar = True
        plotitem.colorbar_shrink = shrink
        plotitem.colorbar_label = "Manning's-$n$ Coefficient"
        plotitem.amr_celledges_show = [0,0,0]
        plotitem.amr_patchedges_show = [1,1,1,1,1,0,0]

def add_wind(plotaxes,bounds=None,plot_type='pcolor',shrink=1.0):
    if plot_type == 'pcolor' or plot_type == 'imshow':
        plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
        plotitem.plot_var = wind_speed
        plotitem.pcolor_cmap = plt.get_cmap('PuBu')
        if bounds is not None:
            plotitem.pcolor_cmin = bounds[0]
            plotitem.pcolor_cmax = bounds[1]
        plotitem.add_colorbar = True
        plotitem.colorbar_shrink = shrink
        plotitem.colorbar_label = "Wind Speed (m/s)"
        plotitem.amr_celledges_show = [0,0,0]
        plotitem.amr_patchedges_show = [1,1,1,1,1,0,0]
    elif plot_type == 'contour':
        plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
        plotitem.plot_var = wind_speed
        plotitem.contour_nlevels = len(surge_data.wind_refine)
        plotitem.countour_min = surge_data.wind_refine[0]
        plotitem.patchedges_show = 1
    elif plot_type == 'quiver':
        plotitem = plotaxes.new_plotitem(plot_type='2d_quiver')
        plotitem.quiver_var_x = wind_x
        plotitem.quiver_var_y = wind_y
        plotitem.amr_quiver_show = [0,0,1]
        plotitem.amr_quiver_key_show = [True,False,False]
        plotitem.amr_quiver_key_units = 'm/s'
        
def add_pressure(plotaxes, bounds=None, plot_type='pcolor', shrink=1.0):
    if plot_type == 'pcolor' or plot_type == 'imshow':
        plotitem = plotaxes.new_plotitem(name="pressure", plot_type='2d_pcolor')
        plotitem.plot_var = pressure
        plotitem.colorbar_shrink = shrink
        plotitem.pcolor_cmap = plt.get_cmap('PuBu')
        if bounds is not None:
            plotitem.pcolor_cmin = bounds[0]
            plotitem.pcolor_cmax = bounds[1]
        plotitem.add_colorbar = True
        plotitem.colorbar_shrink = shrink
        plotitem.colorbar_label = "Pressure (mbar)"
        plotitem.amr_celledges_show = [0,0,0]
        plotitem.amr_patchedges_show = [1,1,1,1,1,0,0]
    elif plot_type == 'contour':
        pass
        
def add_vorticity(plotaxes,bounds=None,plot_type="pcolor"):
    if plot_type == 'pcolor' or plot_type == 'imshow':            
        plotitem = plotaxes.new_plotitem(name="vorticity", plot_type='2d_pcolor')
        plotitem.plot_var = 9
        plotitem.pcolor_cmap = plt.get_cmap('PRGn')
        if bounds is not None:
            plotitem.pcolor_cmin = bounds[0]
            plotitem.pcolor_cmax = bounds[1]
        plotitem.add_colorbar = True
        plotitem.amr_celledges_show = [0,0,0]
        plotitem.amr_patchedges_show = [1]
        
def add_land(plotaxes,plot_type='pcolor',topo_min=-10,topo_max=10.0):
    if plot_type == 'pcolor':
        plotitem = plotaxes.new_plotitem(name='land',plot_type='2d_pcolor')
        plotitem.show = True
        cmap = colormaps.make_colormap({-1:[0.3,0.2,0.1],
                                       -0.00001:[0.95,0.9,0.7],
                                       0.00001:[.5,.7,0],
                                       1:[.2,.5,.2]})
        # color_norm = colors.Normalize(topo_min,topo_max,clip=True)
        plotitem.plot_var = geoplot.land
        # plotitem.pcolor_cmap = geoplot.land_colors
        plotitem.pcolor_cmap = cmap
        plotitem.pcolor_cmin = topo_min
        plotitem.pcolor_cmax = topo_max
        plotitem.add_colorbar = False
        plotitem.amr_celledges_show = [0,0,0,0,0,0,0]
        plotitem.amr_patchedges_show = [1,1,1,1,1,0,0]
    elif plot_type == 'contour':            
        plotitem = plotaxes.new_plotitem(name="land", plot_type='2d_contour')
        plotitem.plot_var = geoplot.land
        plotitem.contour_nlevels = 40
        plotitem.contour_min = 0.0
        plotitem.contour_max = 100.0
        plotitem.amr_contour_colors = ['g']  # color on each level
        plotitem.amr_patch_bgcolor = ['#ffeeee', '#eeeeff', '#eeffee']
        plotitem.celledges_show = 0
        plotitem.patchedges_show = 0   

def add_bathy_contours(plotaxes,contour_levels=None,color='k'):
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.plot_var = geoplot.topo
    if contour_levels is None:
        contour_levels = [0.0]
    plotitem.contour_levels = contour_levels
    plotitem.amr_contour_colors = [color]  # color on each level
    plotitem.kwargs = {'linestyles':'solid','linewidths':2}
    plotitem.amr_contour_show = [0,0,0,0,0,0,1]
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0


# ===== Storm related plotting =======
sec2days = lambda seconds:seconds / (60.0**2 * 24.0)

def plot_track(t,x,y,wind_radius,wind_speed,Pc):
    r"""Plot hurricane track given a storm.data file"""
    colors = ['r','b']

    divide = (np.max(Pc) + np.min(Pc)) / 2.0

    fig = plt.figure(1)
    axes = fig.add_subplot(111)
    indices = Pc < divide
    axes.scatter(x[indices],y[indices],color='r',marker='o')
    indices = Pc >= divide
    axes.scatter(x[indices],y[indices],color='b',marker='o')
    axes.set_title("Track - Hurricane Ike")

    fig = plt.figure(2,figsize=(8*3,6))
    axes = fig.add_subplot(131)
    axes.plot(sec2days(t),wind_speed)
    axes.set_title("Maximum Wind Speed - Hurricane Ike")

    axes = fig.add_subplot(132)
    axes.plot(sec2days(t),wind_radius)
    axes.set_title("Maximum Wind Radius - Hurricane Ike")

    axes = fig.add_subplot(133)
    axes.plot(sec2days(t),Pc)
    axes.plot(sec2days(t),np.ones(t.shape) * divide,'k--')
    axes.set_title("Central Pressure - Hurricane Ike")

    # data = np.loadtxt('_output/storm.data')
    # plot_track(data[:,0],data[:,1],data[:,2],data[:,-3],data[:,-2],data[:,-1])

    # plt.show()
