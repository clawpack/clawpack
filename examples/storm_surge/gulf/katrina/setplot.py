
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
"""

import os

import numpy as np
import matplotlib

import matplotlib.pyplot as plt

from geoclaw import topotools
from clawpack.visclaw import colormaps, geoplot, gaugetools
from clawpack.clawutil.oldclawdata import Data

import geoclaw.surge as surge

try:
    from setplotfg import setplotfg
except:
    setplotfg = None

def setplot(plotdata):
    r"""Setplot function for surge plotting"""
    

    plotdata.clearfigures()  # clear any old figures,axes,items data

    fig_num_counter = surge.plot.figure_counter()

    # Load data from output
    amrdata = Data(os.path.join(plotdata.outdir,'amr2ez.data'))
    physics = Data(os.path.join(plotdata.outdir,'physics.data'))
    surge_data = Data(os.path.join(plotdata.outdir,'surge.data'))

    # Load storm track
    track = surge.plot.track_data(os.path.join(plotdata.outdir,'fort.track'))
    surge_afteraxes = lambda cd: surge.plot.surge_afteraxes(cd,track)

    # Limits for plots
    full_xlimits = [-99.0,-50.0]
    full_ylimits = [8.0,32.0]

    # Color limits
    surface_range = 6.0
    speed_range = 6.0

    xlimits = full_xlimits
    ylimits = full_ylimits
    eta = physics.eta_init
    if not isinstance(eta,list):
        eta = [eta]
    surface_limits = [eta[0]-surface_range,eta[0]+surface_range]
    speed_limits = [0.0,speed_range]
    
    wind_limits = [0,55]
    pressure_limits = [966,1013]
    vorticity_limits = [-1.e-2,1.e-2]

    def pcolor_afteraxes(current_data):
        surge_afteraxes(current_data)
        surge.plot.gauge_locations(current_data)
    
    def contour_afteraxes(current_data):
        surge_afteraxes(current_data)

    # ==========================================================================
    # ==========================================================================
    #   Plot specifications
    # ==========================================================================
    # ==========================================================================

    # ========================================================================
    #  Surface Elevations - Entire Gulf
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='Surface - Entire Domain', 
                                         figno=fig_num_counter.get_counter())
    plotfigure.show = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Surface'
    plotaxes.scaled = True
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.afteraxes = pcolor_afteraxes
    
    surge.plot.add_surface_elevation(plotaxes,bounds=surface_limits)
    surge.plot.add_land(plotaxes)


    # ========================================================================
    #  Water Speed - Entire Gulf
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='Currents - Entire Domain',  
                                         figno=fig_num_counter.get_counter())
    plotfigure.show = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Currents'
    plotaxes.scaled = True
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.afteraxes = pcolor_afteraxes

    # Speed
    surge.plot.add_speed(plotaxes,bounds=speed_limits)

    # Land
    surge.plot.add_land(plotaxes)

    # # ========================================================================
    # #  Surface Elevations - Florida
    # # ========================================================================
    # plotfigure = plotdata.new_plotfigure(name='Surface - Houston/Galveston',  
    #                                      figno=fig_num_counter.get_counter())
    # plotfigure.show = False

    # # Set up for axes in this figure:
    # plotaxes = plotfigure.new_plotaxes()
    # plotaxes.title = 'Surface'
    # plotaxes.scaled = True
    # plotaxes.xlimits = houston_xlimits
    # plotaxes.ylimits = houston_ylimits
    # plotaxes.afteraxes = pcolor_afteraxes
    
    # surge.plot.add_surface_elevation(plotaxes,bounds=surface_limits)
    # surge.plot.add_land(plotaxes)


    # # ========================================================================
    # #  Water Speed - Florida
    # # ========================================================================
    # plotfigure = plotdata.new_plotfigure(name='Currents - Houston/Galveston',  
    #                                      figno=fig_num_counter.get_counter())
    # plotfigure.show = False

    # # Set up for axes in this figure:
    # plotaxes = plotfigure.new_plotaxes()
    # plotaxes.title = 'Currents'
    # plotaxes.scaled = True
    # plotaxes.xlimits = houston_xlimits
    # plotaxes.ylimits = houston_ylimits
    # plotaxes.afteraxes = pcolor_afteraxes

    # # Speed
    # surge.plot.add_speed(plotaxes,bounds=speed_limits)

    # # Land
    # surge.plot.add_land(plotaxes)

    # ========================================================================
    # Hurricane forcing - Entire gulf
    # ========================================================================
    # Pressure field
    plotfigure = plotdata.new_plotfigure(name='Pressure',  
                                         figno=fig_num_counter.get_counter())
    plotfigure.show = surge_data.pressure_forcing
    
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = full_xlimits
    plotaxes.ylimits = full_ylimits
    plotaxes.title = "Pressure Field"
    plotaxes.afteraxes = surge_afteraxes
    plotaxes.scaled = True
    
    surge.plot.add_pressure(plotaxes,bounds=pressure_limits)
    # add_pressure(plotaxes)
    surge.plot.add_land(plotaxes)

    # Pressure gradient
    dp = 1.e3
    plotfigure = plotdata.new_plotfigure(name='Pressure Gradient', 
                                         figno=fig_num_counter.get_counter())
    plotfigure.show = surge_data.pressure_forcing and False
    plotfigure.kwargs = {'figsize':(16,6)}
    
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = "subplot(121)"
    plotaxes.xlimits = full_xlimits
    plotaxes.ylimits = full_ylimits
    plotaxes.title = "X-Component of Pressure Gradient"
    plotaxes.afteraxes = surge_afteraxes
    plotaxes.scaled = True

    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = surge.plot.pressure_gradient_x
    plotitem.imshow_cmap = colormaps.make_colormap({1.0:'r',0.5:'w',0.0:'b'})
    plotitem.imshow_cmin = -dp
    plotitem.imshow_cmax = dp
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [1,1,1]
    
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = "subplot(122)"
    plotaxes.xlimits = full_xlimits
    plotaxes.ylimits = full_ylimits
    plotaxes.title = "Y-Component of Pressure Gradient"
    plotaxes.afteraxes = surge_afteraxes
    plotaxes.scaled = True

    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = surge.plot.pressure_gradient_y
    plotitem.imshow_cmap = colormaps.make_colormap({1.0:'r',0.5:'w',0.0:'b'})
    plotitem.imshow_cmin = -dp
    plotitem.imshow_cmax = dp
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [1,1,1]
    
    # Wind field
    plotfigure = plotdata.new_plotfigure(name='Wind Speed', 
                                         figno=fig_num_counter.get_counter())
    plotfigure.show = surge_data.wind_forcing
    
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = full_xlimits
    plotaxes.ylimits = full_ylimits
    plotaxes.title = "Wind Field"
    plotaxes.afteraxes = surge_afteraxes
    plotaxes.scaled = True
    
    surge.plot.add_wind(plotaxes,bounds=wind_limits,plot_type='imshow')
    # add_wind(plotaxes,bounds=wind_limits,plot_type='contour')
    # add_wind(plotaxes,bounds=wind_limits,plot_type='quiver')
    surge.plot.add_land(plotaxes)
    
    # Wind field components
    plotfigure = plotdata.new_plotfigure(name='Wind Components', 
                                         figno=fig_num_counter.get_counter())
    plotfigure.show = surge_data.wind_forcing and False
    
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = "subplot(121)"
    plotaxes.xlimits = full_xlimits
    plotaxes.ylimits = full_ylimits
    plotaxes.title = "X-Component of Wind Field"
    plotaxes.afteraxes = surge_afteraxes
    plotaxes.scaled = True

    surge.plot.add_land(plotaxes)

    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = surge.plot.wind_x
    plotitem.imshow_cmap = colormaps.make_colormap({1.0:'r',0.5:'w',0.0:'b'})
    plotitem.imshow_cmin = -wind_limits[1]
    plotitem.imshow_cmax = wind_limits[1]
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [1,1,1]
    
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = "subplot(122)"
    plotaxes.xlimits = full_xlimits
    plotaxes.ylimits = full_ylimits
    plotaxes.title = "Y-Component of Wind Field"
    plotaxes.afteraxes = surge_afteraxes
    plotaxes.scaled = True

    surge.plot.add_land(plotaxes)

    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = surge.plot.wind_y
    plotitem.imshow_cmap = colormaps.make_colormap({1.0:'r',0.5:'w',0.0:'b'})
    plotitem.imshow_cmin = -wind_limits[1]
    plotitem.imshow_cmax = wind_limits[1]
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [1,1,1]

    # ========================================================================
    #  Figures for gauges
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='Surface & topo', figno=300, \
                    type='each_gauge')
    plotfigure.show = True
    plotfigure.clf_each_gauge = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    try:
        plotaxes.xlimits = [amrdata.t0,amrdata.tfinal]
    except:
        pass
    # plotaxes.ylimits = [0,150.0]
    plotaxes.ylimits = surface_limits
    plotaxes.title = 'Surface'
    plotaxes.afteraxes = surge.plot.gauge_afteraxes

    # Plot surface as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 3
    plotitem.plotstyle = 'r-'

    # # ==========================================================================
    # # Galveston Channel
    # # ==========================================================================
    # plotfigure = plotdata.new_plotfigure(name='Galveston Bay Entry Channel', figno=2)

    # # Set up for axes in this figure:
    # plotaxes = plotfigure.new_plotaxes('pcolor')
    # plotaxes.title = 'Surface'
    # plotaxes.scaled = True
    # plotaxes.afteraxes = fixup

    # # Water
    # plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    # # plotitem.plot_var = geoplot.surface
    # plotitem.plot_var = geoplot.surface_or_depth
    # plotitem.pcolor_cmap = geoplot.tsunami_colormap
    # plotitem.pcolor_cmin = -0.2
    # plotitem.pcolor_cmax = 0.2
    # plotitem.add_colorbar = True
    # plotitem.amr_celledges_show = [0,0,0]
    # plotitem.patchedges_show = 1

    # # Land
    # plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    # plotitem.plot_var = geoplot.land
    # plotitem.pcolor_cmap = geoplot.land_colors
    # plotitem.pcolor_cmin = 0.0
    # plotitem.pcolor_cmax = 100.0
    # plotitem.add_colorbar = False
    # # plotitem.amr_celledges_show = [1,1,0]
    # plotitem.amr_celledges_show = [0,0,0]
    # plotitem.patchedges_show = 1
    # plotaxes.xlimits = galveston_xlimits
    # plotaxes.ylimits = galveston_ylimits

    # # ==========================================================================
    # #  City of Houston Area
    # # ==========================================================================
    # plotfigure = plotdata.new_plotfigure(name='City of Houston', figno=3)

    # # Set up for axes in this figure:
    # plotaxes = plotfigure.new_plotaxes('pcolor')
    # plotaxes.title = 'Surface'
    # plotaxes.scaled = True
    # plotaxes.afteraxes = fixup

    # # Water
    # plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    # # plotitem.plot_var = geoplot.surface
    # plotitem.plot_var = geoplot.surface_or_depth
    # plotitem.pcolor_cmap = geoplot.tsunami_colormap
    # plotitem.pcolor_cmin = -0.2
    # plotitem.pcolor_cmax = 0.2
    # plotitem.add_colorbar = True
    # plotitem.amr_celledges_show = [0,0,0]
    # plotitem.patchedges_show = 1

    # # Land
    # plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    # plotitem.plot_var = geoplot.land
    # plotitem.pcolor_cmap = geoplot.land_colors
    # plotitem.pcolor_cmin = 0.0
    # plotitem.pcolor_cmax = 100.0
    # plotitem.add_colorbar = False
    # # plotitem.amr_celledges_show = [1,1,0]
    # plotitem.amr_celledges_show = [0,0,0]
    # plotitem.patchedges_show = 1
    # plotaxes.xlimits = houston_harbor_xlimits
    # plotaxes.ylimits = houston_harbor_ylimits


    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    # plotfigure = plotdata.new_plotfigure(name='Surface & topo', figno=300, \
    #                 type='each_gauge')
    # plotfigure.clf_each_gauge = True
    # 
    # # Set up for axes in this figure:
    # plotaxes = plotfigure.new_plotaxes()
    # plotaxes.xlimits = 'auto'
    # plotaxes.ylimits = 'auto'
    # plotaxes.title = 'Surface'
    # 
    # # Plot surface as blue curve:
    # plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    # plotitem.plot_var = 3
    # plotitem.plotstyle = 'b-'
    # 
    # # Plot topo as green curve:
    # plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    # plotitem.show = False
    # 
    # def gaugetopo(current_data):
    #     q = current_data.q
    #     h = q[0,:]
    #     eta = q[3,:]
    #     topo = eta - h
    #     return topo
    #     
    # plotitem.plot_var = gaugetopo
    # plotitem.plotstyle = 'g-'
    # 
    # def add_zeroline(current_data):
    #     from pylab import plot, legend, xticks, floor
    #     t = current_data.t
    #     #legend(('surface','topography'),loc='lower left')
    #     plot(t, 0*t, 'k')
    #     n = int(floor(t.max()/3600.) + 2)
    #     xticks([3600*i for i in range(n)])
    # 
    # plotaxes.afteraxes = add_zeroline


    #-----------------------------------------
    
    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_gaugenos = 'all'          # list of gauges to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?

    return plotdata

