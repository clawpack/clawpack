
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
"""

article = False

import os

import numpy

# Plot customization
import matplotlib

# Use LaTeX for all text
matplotlib.rcParams['text.usetex'] = True

# Markers and line widths
matplotlib.rcParams['lines.linewidth'] = 2.0
matplotlib.rcParams['lines.markersize'] = 6
matplotlib.rcParams['lines.markersize'] = 8

# Font Sizes
matplotlib.rcParams['font.size'] = 16
matplotlib.rcParams['axes.labelsize'] = 16
matplotlib.rcParams['legend.fontsize'] = 12
matplotlib.rcParams['xtick.labelsize'] = 16
matplotlib.rcParams['ytick.labelsize'] = 16

# DPI of output images
matplotlib.rcParams['savefig.dpi'] = 100

import matplotlib.pyplot as plt
import datetime

from clawpack.visclaw import colormaps
import clawpack.clawutil.data as clawutil
import clawpack.amrclaw.data as amrclaw
import clawpack.geoclaw.data as geodata

import clawpack.geoclaw.surge.plot as surgeplot

try:
    from setplotfg import setplotfg
except:
    setplotfg = None

# Gauge support
days2seconds = lambda days: days * 60.0**2 * 24.0
date2seconds = lambda date: days2seconds(date.days) + date.seconds
seconds2days = lambda secs: secs / (24.0 * 60.0**2)
min2deg = lambda minutes: minutes / 60.0
ft2m = lambda x:0.3048 * x

# Gauge name translation
gauge_name_trans = {1:"W", 2:"X", 3:"Y", 4:"Z"}
gauge_surface_offset = [0.0, 0.0]
gauge_landfall = []
gauge_landfall.append(datetime.datetime(2008,9,13 + 1,7) 
                                            - datetime.datetime(2008,1,1,0))
gauge_landfall.append(datetime.datetime(2008,9,13 - 1,7) 
                                            - datetime.datetime(2008,1,1,0))
gauge_landfall.append(days2seconds(4.25))

def setplot(plotdata):
    r"""Setplot function for surge plotting"""
    

    plotdata.clearfigures()  # clear any old figures,axes,items data
    plotdata.format = 'binary'

    fig_num_counter = surgeplot.figure_counter()

    # Load data from output
    clawdata = clawutil.ClawInputData(2)
    clawdata.read(os.path.join(plotdata.outdir,'claw.data'))
    amrdata = amrclaw.AmrclawInputData(clawdata)
    amrdata.read(os.path.join(plotdata.outdir,'amr.data'))
    physics = geodata.GeoClawData()
    physics.read(os.path.join(plotdata.outdir,'geoclaw.data'))
    surge_data = geodata.SurgeData()
    surge_data.read(os.path.join(plotdata.outdir,'surge.data'))
    friction_data = geodata.FrictionData()
    friction_data.read(os.path.join(plotdata.outdir,'friction.data'))

    # Load storm track
    track = surgeplot.track_data(os.path.join(plotdata.outdir,'fort.track'))

    # Calculate landfall time, off by a day, maybe leap year issue?
    landfall_dt = datetime.datetime(2008,9,13,7) - datetime.datetime(2008,1,1,0)
    landfall = (landfall_dt.days - 1.0) * 24.0 * 60**2 + landfall_dt.seconds

    # Set afteraxes function
    surge_afteraxes = lambda cd: surgeplot.surge_afteraxes(cd,
                                        track, landfall, plot_direction=False)

    # Color limits
    surface_range = 5.0
    speed_range = 3.0
    eta = physics.sea_level
    if not isinstance(eta,list):
        eta = [eta]
    surface_limits = [eta[0]-surface_range,eta[0]+surface_range]
    # surface_contours = numpy.linspace(-surface_range, surface_range,11)
    surface_contours = [-5,-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5]
    surface_ticks = [-5,-4,-3,-2,-1,0,1,2,3,4,5]
    surface_labels = [str(value) for value in surface_ticks]
    speed_limits = [0.0,speed_range]
    speed_contours = numpy.linspace(0.0,speed_range,13)
    speed_ticks = [0,1,2,3]
    speed_labels = [str(value) for value in speed_ticks]
    
    wind_limits = [0,64]
    # wind_limits = [-0.002,0.002]
    pressure_limits = [935,1013]
    friction_bounds = [0.01,0.04]
    # vorticity_limits = [-1.e-2,1.e-2]

    # def pcolor_afteraxes(current_data):
    #     surge_afteraxes(current_data)
    #     surge.plot.gauge_locations(current_data,gaugenos=[6])
    
    def contour_afteraxes(current_data):
        surge_afteraxes(current_data)

    def add_custom_colorbar_ticks_to_axes(axes, item_name, ticks, tick_labels=None):
        axes.plotitem_dict[item_name].colorbar_ticks = ticks
        axes.plotitem_dict[item_name].colorbar_tick_labels = tick_labels

    # ==========================================================================
    # ==========================================================================
    #   Plot specifications
    # ==========================================================================
    # ==========================================================================

    # ========================================================================
    #  Entire Gulf
    # ========================================================================
    gulf_xlimits = [clawdata.lower[0],clawdata.upper[0]]
    gulf_ylimits = [clawdata.lower[1],clawdata.upper[1]]
    gulf_shrink = 0.9
    def gulf_after_axes(cd):
        plt.subplots_adjust(left=0.08, bottom=0.04, right=0.97, top=0.96)
        surge_afteraxes(cd)
    #
    #  Surface
    #
    plotfigure = plotdata.new_plotfigure(name='Surface - Entire Domain', 
                                         figno=fig_num_counter.get_counter())
    plotfigure.show = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Surface'
    plotaxes.scaled = True
    plotaxes.xlimits = gulf_xlimits
    plotaxes.ylimits = gulf_ylimits
    plotaxes.afteraxes = gulf_after_axes

    surgeplot.add_surface_elevation(plotaxes, plot_type='contourf',
                                               contours=surface_contours,
                                               shrink=gulf_shrink)
    surgeplot.add_land(plotaxes,topo_min=-10.0,topo_max=5.0)
    # surge.plot.add_bathy_contours(plotaxes)
    if article:
        plotaxes.plotitem_dict['surface'].add_colorbar = False
    else:
        add_custom_colorbar_ticks_to_axes(plotaxes, 'surface', surface_ticks, surface_labels)
    plotaxes.plotitem_dict['surface'].amr_patchedges_show = [1,1,1,1,1,1,1,1]

    #
    #  Water Speed
    #
    plotfigure = plotdata.new_plotfigure(name='Currents - Entire Domain',  
                                         figno=fig_num_counter.get_counter())
    plotfigure.show = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Currents'
    plotaxes.scaled = True
    plotaxes.xlimits = gulf_xlimits
    plotaxes.ylimits = gulf_ylimits
    plotaxes.afteraxes = gulf_after_axes

    # Speed
    surgeplot.add_speed(plotaxes, plot_type='contourf',
                                   contours=speed_contours,
                                   shrink=gulf_shrink)
    if article:
        plotaxes.plotitem_dict['speed'].add_colorbar = False
    else:
        add_custom_colorbar_ticks_to_axes(plotaxes, 'speed', speed_ticks, speed_labels)

    # Land
    surgeplot.add_land(plotaxes)
    surgeplot.add_bathy_contours(plotaxes)

    #
    # Friction field
    #
    plotfigure = plotdata.new_plotfigure(name='Friction',
                                         figno=fig_num_counter.get_counter())
    plotfigure.show = friction_data.variable_friction and True

    def friction_after_axes(cd):
        plt.subplots_adjust(left=0.08, bottom=0.04, right=0.97, top=0.96)
        plt.title(r"Manning's $n$ Coefficient")
        # surge_afteraxes(cd)

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = gulf_xlimits
    plotaxes.ylimits = gulf_ylimits
    # plotaxes.title = "Manning's N Coefficient"
    plotaxes.afteraxes = friction_after_axes
    plotaxes.scaled = True

    surgeplot.add_friction(plotaxes,bounds=friction_bounds,shrink=0.9)
    plotaxes.plotitem_dict['friction'].amr_patchedges_show = [0,0,0,0,0,0,0]
    plotaxes.plotitem_dict['friction'].colorbar_label = "$n$"


    # ========================================================================
    #  LaTex Shelf
    # ========================================================================
    latex_xlimits = [-97.5,-88.5]
    latex_ylimits = [27.5,30.5]
    latex_shrink = 1.0
    def latex_after_axes(cd):
        if article:
            plt.subplots_adjust(left=0.07, bottom=0.14, right=1.0, top=0.86)
        else:
            plt.subplots_adjust(right=1.0)
        surge_afteraxes(cd)

    #
    # Surface
    #
    plotfigure = plotdata.new_plotfigure(name='Surface - LaTex Shelf', 
                                         figno=fig_num_counter.get_counter())
    plotfigure.show = True
    if article:
        plotfigure.kwargs = {'figsize':(8,2.7), 'facecolor':'none'}
    else:
        plotfigure.kwargs = {'figsize':(9,2.7), 'facecolor':'none'}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Surface'
    plotaxes.scaled = True
    plotaxes.xlimits = latex_xlimits
    plotaxes.ylimits = latex_ylimits
    plotaxes.afteraxes = latex_after_axes
    surgeplot.add_surface_elevation(plotaxes, plot_type='contourf',
                                               contours=surface_contours,
                                               shrink=latex_shrink)

    if article:
        plotaxes.plotitem_dict['surface'].add_colorbar = False
        # plotaxes.afteraxes = lambda cd: article_latex_after_axes(cd, landfall)
    else:
        add_custom_colorbar_ticks_to_axes(plotaxes, 'surface', [-5,-2.5,0,2.5,5.0], 
                                    ["-5.0","-2.5"," 0"," 2.5"," 5.0"])
    # plotaxes.plotitem_dict['surface'].contour_cmap = plt.get_cmap('OrRd')
    # surge.plot.add_surface_elevation(plotaxes,plot_type='contour')
    surgeplot.add_land(plotaxes)
    # plotaxes.plotitem_dict['surface'].amr_patchedges_show = [1,1,1,0,0,0,0]
    plotaxes.plotitem_dict['surface'].amr_patchedges_show = [0,0,0,0,0,0,0]
    # plotaxes.plotitem_dict['land'].amr_patchedges_show = [1,1,1,0,0,0,0]
    plotaxes.plotitem_dict['land'].amr_patchedges_show = [0,0,0,0,0,0,0]

    # Plot using jet and 0.0 to 5.0 to match figgen generated ADCIRC results
    # plotaxes.plotitem_dict['surface'].pcolor_cmin = 0.0
    # plotaxes.plotitem_dict['surface'].pcolor_cmax = 5.0
    # plotaxes.plotitem_dict['surface'].pcolor_cmap = plt.get_cmap('jet')

    #
    # Water Speed
    #
    plotfigure = plotdata.new_plotfigure(name='Currents - LaTex Shelf',  
                                         figno=fig_num_counter.get_counter())
    plotfigure.show = True
    if article:
        plotfigure.kwargs = {'figsize':(8,2.7), 'facecolor':'none'}
    else:
        plotfigure.kwargs = {'figsize':(9,2.7), 'facecolor':'none'}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Currents'
    plotaxes.scaled = True
    plotaxes.xlimits = latex_xlimits
    plotaxes.ylimits = latex_ylimits
    plotaxes.afteraxes = latex_after_axes
    surgeplot.add_speed(plotaxes, plot_type='contourf',
                                   contours=speed_contours,
                                   shrink=latex_shrink)

    if article:
        plotaxes.plotitem_dict['speed'].add_colorbar = False
    else:
        add_custom_colorbar_ticks_to_axes(plotaxes, 'speed', speed_ticks, speed_labels)
    # surge.plot.add_surface_elevation(plotaxes,plot_type='contour')
    surgeplot.add_land(plotaxes)
    # plotaxes.plotitem_dict['speed'].amr_patchedges_show = [1,1,0,0,0,0,0]
    # plotaxes.plotitem_dict['land'].amr_patchedges_show = [1,1,1,0,0,0,0]
    plotaxes.plotitem_dict['speed'].amr_patchedges_show = [0,0,0,0,0,0,0]
    plotaxes.plotitem_dict['land'].amr_patchedges_show = [0,0,0,0,0,0,0]



    # ========================================================================
    #  Houston/Galveston
    # ========================================================================
    houston_xlimits = [-(95.0 + 26.0 / 60.0), -(94.0 + 25.0 / 60.0)]
    houston_ylimits = [29.1, 29.0 + 55.0 / 60.0]
    houston_shrink = 0.9
    def houston_after_axes(cd):
        if article:
            plt.subplots_adjust(left=0.05, bottom=0.07, right=0.99, top=0.92)
        else:
            plt.subplots_adjust(left=0.12, bottom=0.06, right=0.97, top=0.97)
        surge_afteraxes(cd)
        # surge.plot.gauge_locations(cd)
    
    #
    # Surface Elevations
    #
    plotfigure = plotdata.new_plotfigure(name='Surface - Houston/Galveston',  
                                         figno=fig_num_counter.get_counter())
    plotfigure.show = True
    # if article:
    #     plotfigure.kwargs['figsize'] = 

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Surface'
    plotaxes.scaled = True
    plotaxes.xlimits = houston_xlimits
    plotaxes.ylimits = houston_ylimits
    plotaxes.afteraxes = houston_after_axes
    surgeplot.add_surface_elevation(plotaxes, plot_type='contourf',
                                               contours=surface_contours,
                                               shrink=houston_shrink)
    
    if article:
        plotaxes.plotitem_dict['surface'].add_colorbar = False
    else:
        add_custom_colorbar_ticks_to_axes(plotaxes, 'surface', surface_ticks, surface_labels)
    surgeplot.add_land(plotaxes)
    plotaxes.plotitem_dict['surface'].amr_patchedges_show = [0,0,0,0,0,0,0]
    plotaxes.plotitem_dict['land'].amr_patchedges_show = [0,0,0,0,0,0,0]
    # surge.plot.add_bathy_contours(plotaxes)

    # Plot using jet and 0.0 to 5.0 to match figgen generated ADCIRC results
    # plotaxes.plotitem_dict['surface'].pcolor_cmin = 0.0
    # plotaxes.plotitem_dict['surface'].pcolor_cmax = 5.0
    # plotaxes.plotitem_dict['surface'].pcolor_cmap = plt.get_cmap('jet')

    #
    # Water Speed
    #
    plotfigure = plotdata.new_plotfigure(name='Currents - Houston/Galveston',  
                                         figno=fig_num_counter.get_counter())
    plotfigure.show = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Currents'
    plotaxes.scaled = True
    plotaxes.xlimits = houston_xlimits
    plotaxes.ylimits = houston_ylimits
    plotaxes.afteraxes = houston_after_axes
    surgeplot.add_speed(plotaxes, plot_type='contourf',
                                   contours=speed_contours,
                                   shrink=houston_shrink)
    
    if article:
        plotaxes.plotitem_dict['speed'].add_colorbar = False
    else:
        add_custom_colorbar_ticks_to_axes(plotaxes, 'speed', speed_ticks, speed_labels)
    surgeplot.add_land(plotaxes)
    # surge.plot.add_bathy_contours(plotaxes)
    # plotaxes.plotitem_dict['speed'].amr_patchedges_show = [1,1,1,1,1,1,1,1]
    # plotaxes.plotitem_dict['land'].amr_patchedges_show = [1,1,1,1,1,1,1,1]
    plotaxes.plotitem_dict['speed'].amr_patchedges_show = [0,0,0,0,0,0,0]
    plotaxes.plotitem_dict['land'].amr_patchedges_show = [0,0,0,0,0,0,0]

    # ==========================
    #  Hurricane Forcing fields
    # ==========================
    
    # Pressure field
    plotfigure = plotdata.new_plotfigure(name='Pressure',  
                                         figno=fig_num_counter.get_counter())
    plotfigure.show = surge_data.pressure_forcing and True
    
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = gulf_xlimits
    plotaxes.ylimits = gulf_ylimits
    plotaxes.title = "Pressure Field"
    plotaxes.afteraxes = gulf_after_axes
    plotaxes.scaled = True
    surgeplot.add_pressure(plotaxes, bounds=pressure_limits, shrink=gulf_shrink)
    surgeplot.add_land(plotaxes)

    # Wind field
    plotfigure = plotdata.new_plotfigure(name='Wind Speed', 
                                         figno=fig_num_counter.get_counter())
    plotfigure.show = surge_data.wind_forcing and True
    
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = gulf_xlimits
    plotaxes.ylimits = gulf_ylimits
    plotaxes.title = "Wind Field"
    plotaxes.afteraxes = gulf_after_axes
    plotaxes.scaled = True
    surgeplot.add_wind(plotaxes, bounds=wind_limits, plot_type='pcolor',
                                  shrink=gulf_shrink)
    surgeplot.add_land(plotaxes)

    # ========================================================================
    #  Figures for gauges
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='Surface & topo', figno=300, \
                    type='each_gauge')
    plotfigure.show = True
    plotfigure.clf_each_gauge = True
    # plotfigure.kwargs['figsize'] = (16,10)

    def gauge_after_axes(cd):

        if cd.gaugeno in [1,2,3,4]:
            axes = plt.gca()
            # # Add Kennedy gauge data
            # kennedy_gauge = kennedy_gauges[gauge_name_trans[cd.gaugeno]]
            # axes.plot(kennedy_gauge['t'] - seconds2days(date2seconds(gauge_landfall[0])), 
            #          kennedy_gauge['mean_water'] + kennedy_gauge['depth'], 'k-', 
            #          label='Gauge Data')

            # Add GeoClaw gauge data
            geoclaw_gauge = cd.gaugesoln
            axes.plot(seconds2days(geoclaw_gauge.t - date2seconds(gauge_landfall[1])),
                  geoclaw_gauge.q[3,:] + gauge_surface_offset[0], 'b--', 
                  label="GeoClaw")

            # Add ADCIRC gauge data
            # ADCIRC_gauge = ADCIRC_gauges[kennedy_gauge['gauge_no']]
            # axes.plot(seconds2days(ADCIRC_gauge[:,0] - gauge_landfall[2]), 
            #          ADCIRC_gauge[:,1] + gauge_surface_offset[1], 'r-.', label="ADCIRC")

            # Fix up plot
            axes.set_title('Station %s' % cd.gaugeno)
            axes.set_xlabel('Days relative to landfall')
            axes.set_ylabel('Surface (m)')
            axes.set_xlim([-2,1])
            axes.set_ylim([-1,5])
            axes.set_xticks([-2,-1,0,1])
            axes.set_xticklabels([r"$-2$",r"$-1$",r"$0$",r"$1$"])
            axes.grid(True)
            axes.legend()

            plt.hold(False)

        # surge.plot.gauge_afteraxes(cd)


    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [-2,1]
    # plotaxes.xlabel = "Days from landfall"
    # plotaxes.ylabel = "Surface (m)"
    plotaxes.ylimits = [-1,5]
    plotaxes.title = 'Surface'
    plotaxes.afteraxes = gauge_after_axes

    # Plot surface as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 3
    plotitem.plotstyle = 'b-'

    # =====================
    #  Gauge Location Plot
    # =====================
    gauge_xlimits = [-95.5, -94]
    gauge_ylimits = [29.0, 30.0]
    gauge_location_shrink = 0.75
    def gauge_after_axes(cd):
        plt.subplots_adjust(left=0.12, bottom=0.06, right=0.97, top=0.97)
        surge_afteraxes(cd)
        surgeplot.gauge_locations(cd, gaugenos=[1, 2, 3, 4])
        plt.title("Gauge Locations")

    plotfigure = plotdata.new_plotfigure(name='Gauge Locations',  
                                         figno=fig_num_counter.get_counter())
    plotfigure.show = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Surface'
    plotaxes.scaled = True
    plotaxes.xlimits = gauge_xlimits
    plotaxes.ylimits = gauge_ylimits
    plotaxes.afteraxes = gauge_after_axes
    surgeplot.add_surface_elevation(plotaxes, plot_type='contourf',
                                               contours=surface_contours,
                                               shrink=gauge_location_shrink)
    # surge.plot.add_surface_elevation(plotaxes, plot_type="contourf")
    add_custom_colorbar_ticks_to_axes(plotaxes, 'surface', surface_ticks, surface_labels)
    surgeplot.add_land(plotaxes)
    # plotaxes.plotitem_dict['surface'].amr_patchedges_show = [0,0,0,0,0,0,0]
    # plotaxes.plotitem_dict['surface'].add_colorbar = False
    # plotaxes.plotitem_dict['surface'].pcolor_cmap = plt.get_cmap('jet')
    # plotaxes.plotitem_dict['surface'].pcolor_cmap = plt.get_cmap('gist_yarg')
    # plotaxes.plotitem_dict['surface'].pcolor_cmin = 0.0
    # plotaxes.plotitem_dict['surface'].pcolor_cmax = 5.0
    plotaxes.plotitem_dict['surface'].amr_patchedges_show = [0,0,0,0,0,0,0]
    plotaxes.plotitem_dict['land'].amr_patchedges_show = [0,0,0,0,0,0,0]
    
    # ==============================================================
    #  Debugging Plots, only really work if using interactive plots
    # ==============================================================
    #
    # Water Velocity Components
    #
    plotfigure = plotdata.new_plotfigure(name='Velocity Components - Entire Domain',  
                                         figno=fig_num_counter.get_counter())
    plotfigure.show = False

    # X-Component
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = "subplot(121)"
    plotaxes.title = 'Velocity, X-Component'
    plotaxes.scaled = True
    plotaxes.xlimits = gulf_xlimits
    plotaxes.ylimits = gulf_ylimits
    plotaxes.afteraxes = gulf_after_axes

    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = surgeplot.water_u
    plotitem.pcolor_cmap = colormaps.make_colormap({1.0:'r',0.5:'w',0.0:'b'})
    plotitem.pcolor_cmin = -speed_limits[1]
    plotitem.pcolor_cmax = speed_limits[1]
    plotitem.colorbar_shrink = gulf_shrink
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [1,1,1]

    surgeplot.add_land(plotaxes)

    # Y-Component
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = "subplot(122)"
    plotaxes.title = 'Velocity, Y-Component'
    plotaxes.scaled = True
    plotaxes.xlimits = gulf_xlimits
    plotaxes.ylimits = gulf_ylimits
    plotaxes.afteraxes = gulf_after_axes

    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = surgeplot.water_v
    plotitem.pcolor_cmap = colormaps.make_colormap({1.0:'r',0.5:'w',0.0:'b'})
    plotitem.pcolor_cmin = -speed_limits[1]
    plotitem.pcolor_cmax = speed_limits[1]
    plotitem.colorbar_shrink = gulf_shrink
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [1,1,1]
    surgeplot.add_land(plotaxes)

    # 
    # Depth
    # 
    plotfigure = plotdata.new_plotfigure(name='Depth - Entire Domain', 
                                         figno=fig_num_counter.get_counter())
    plotfigure.show = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'depth'
    plotaxes.scaled = True
    plotaxes.xlimits = gulf_xlimits
    plotaxes.ylimits = gulf_ylimits
    plotaxes.afteraxes = gulf_after_axes

    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = 0
    plotitem.imshow_cmap = colormaps.make_colormap({1.0:'r',0.5:'w',0.0:'b'})
    plotitem.imshow_cmin = 0
    plotitem.imshow_cmax = 100
    plotitem.colorbar_shrink = gulf_shrink
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [1,1,1,1,1,1,1,1,1]
    
    # Surge field
    plotfigure = plotdata.new_plotfigure(name='Surge Field', 
                                         figno=fig_num_counter.get_counter())
    plotfigure.show = ((surge_data.wind_forcing or surge_data.pressure_forcing) 
                        and False)
    
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = gulf_xlimits
    plotaxes.ylimits = gulf_ylimits
    plotaxes.title = "Storm Surge Source Term S"
    plotaxes.afteraxes = gulf_after_axes
    plotaxes.scaled = True
    
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = surgeplot.pressure_field + 1
    plotitem.pcolor_cmap = plt.get_cmap('PuBu')
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 1e-3
    plotitem.add_colorbar = True
    plotitem.colorbar_shrink = gulf_shrink
    plotitem.colorbar_label = "Source Strength"
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [1,1,1,1,1,0,0]
    surgeplot.add_land(plotaxes)

    plotfigure = plotdata.new_plotfigure(name='Friction/Coriolis Source', 
                                         figno=fig_num_counter.get_counter())
    plotfigure.show = False
    
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = gulf_xlimits
    plotaxes.ylimits = gulf_ylimits
    plotaxes.title = "Friction/Coriolis Source"
    plotaxes.afteraxes = surge_afteraxes
    plotaxes.scaled = True
    
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = surgeplot.pressure_field + 2
    plotitem.pcolor_cmap = plt.get_cmap('PuBu')
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 1e-3
    plotitem.add_colorbar = True
    plotitem.colorbar_shrink = gulf_shrink
    plotitem.colorbar_label = "Source Strength"
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [1,1,1,1,1,0,0]
    surgeplot.add_land(plotaxes)

    #-----------------------------------------

    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    if article:
        plotdata.printfigs = True                # print figures
        plotdata.print_format = 'png'            # file format
        plotdata.print_framenos = [54,60,66,72,78,84]            # list of frames to print
        plotdata.print_gaugenos = [1,2,3,4]          # list of gauges to print
        plotdata.print_fignos = [4,5,6,7,10,3,300]            # list of figures to print
        plotdata.html = True                     # create html files of plots?
        plotdata.html_homelink = '../README.html'   # pointer for top of index
        plotdata.latex = False                    # create latex file of plots?
        plotdata.latex_figsperline = 2           # layout of plots
        plotdata.latex_framesperline = 1         # layout of plots
        plotdata.latex_makepdf = False           # also run pdflatex?

    else:
        plotdata.printfigs = True                # print figures
        plotdata.print_format = 'png'            # file format
        plotdata.print_framenos = 'all'            # list of frames to print
        plotdata.print_gaugenos = [1,2,3,4]          # list of gauges to print
        plotdata.print_fignos = 'all'            # list of figures to print
        plotdata.html = True                     # create html files of plots?
        plotdata.html_homelink = '../README.html'   # pointer for top of index
        plotdata.latex = True                    # create latex file of plots?
        plotdata.latex_figsperline = 2           # layout of plots
        plotdata.latex_framesperline = 1         # layout of plots
        plotdata.latex_makepdf = False           # also run pdflatex?

    return plotdata

