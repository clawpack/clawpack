
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt

from clawpack.geoclaw import topotools
from six.moves import range
import os


#--------------------------
def setplot(plotdata=None):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 

    import os

    import numpy as np
    import matplotlib.pyplot as plt

    from clawpack.visclaw import geoplot, gaugetools, colormaps

    import clawpack.clawutil.data as clawutil
    import clawpack.amrclaw.data as amrclaw
    import clawpack.geoclaw.data

    import clawpack.geoclaw.multilayer.plot as ml_plot

    if plotdata is None:
        from clawpack.visclaw.data import ClawPlotData
        plotdata = ClawPlotData()

    from numpy import linspace



    plotdata.clearfigures()  # clear any old figures,axes,items data
    plotdata.save_frames = False

    # Load data from output
    clawdata = clawutil.ClawInputData(2)
    clawdata.read(os.path.join(plotdata.outdir,'claw.data'))
    amrdata = amrclaw.AmrclawInputData(clawdata)
    amrdata.read(os.path.join(plotdata.outdir,'amr.data'))
    geodata = clawpack.geoclaw.data.GeoClawData()
    geodata.read(os.path.join(plotdata.outdir,'geoclaw.data'))
    multilayer_data = clawpack.geoclaw.data.MultilayerData()
    multilayer_data.read(os.path.join(plotdata.outdir,'multilayer.data'))

    # To plot gauge locations on pcolor or contour plot, use this as
    # an afteraxis function:

    def addgauges(current_data):
        from clawpack.visclaw import gaugetools
        gaugetools.plot_gauge_locations(current_data.plotdata, \
             gaugenos='all', format_string='ko', add_labels=True)

    # ========================================================================
    #  Generic helper functions
    # ========================================================================
    def pcolor_afteraxes(current_data):
        # bathy_ref_lines(current_data)
        gauge_locations(current_data)
        
    def contour_afteraxes(current_data):
        # gauge_locations(current_data)
        # m_to_km_labels()
        plt.hold(True)
        pos = -80.0 * (23e3 / 180) + 500e3 - 5e3
        plt.plot([pos,pos],[-300e3,300e3],'b',[pos-5e3,pos-5e3],[-300e3,300e3],'y')
        plt.hold(False)
        wind_contours(current_data)
        bathy_ref_lines(current_data)
        
    def profile_afteraxes(current_data):
        pass
        
    def gauge_locations(current_data,gaugenos='all'):
        plt.hold(True)
        gaugetools.plot_gauge_locations(current_data.plotdata, \
             gaugenos=gaugenos, format_string='kx', add_labels=True)
        plt.hold(False)
    
    # ========================================================================
    # Axis limits
    # xlimits = [amrdata.xlower,amrdata.xupper]
    xlimits = [-100.0, 100.0]

    # ylimits = [amrdata.ylower,amrdata.yupper]
    ylimits = [-100.0, 100.0]
    
    eta = [multilayer_data.eta[0],multilayer_data.eta[1]]

    top_surface_limits = [eta[0]-10,eta[0]+10]
    internal_surface_limits = [eta[1]-5,eta[1]+5]
    depth_limits = [0.0, 0.4]
    top_speed_limits = [0.0,0.1]
    internal_speed_limits = [0.0,0.03]
    

    # ========================================================================
    #  Surface Elevations
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='Surface', figno=0)
    plotfigure.show = True
    plotfigure.kwargs = {'figsize':(14,4)}
    
    # Top surface
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Top Surface'
    plotaxes.axescmd = 'subplot(1,2,1)'
    plotaxes.scaled = True
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    ml_plot.add_inundation(plotaxes, 1, bounds=depth_limits)
    ml_plot.add_surface_elevation(plotaxes,1,bounds=top_surface_limits)
    ml_plot.add_land(plotaxes, 1)
    
    # Bottom surface
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Internal Surface'
    plotaxes.axescmd = 'subplot(1,2,2)'
    plotaxes.scaled = True
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    ml_plot.add_inundation(plotaxes, 2, bounds=depth_limits)
    ml_plot.add_surface_elevation(plotaxes,2,bounds=internal_surface_limits)
    ml_plot.add_colorbar = True
    ml_plot.add_land(plotaxes, 2)


    # ========================================================================
    # Figure for cross section
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='cross-section', figno=4)
    plotfigure.show = True
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.title = 'Cross section at y=0'
    ml_plot.add_cross_section(plotaxes, 1)
    ml_plot.add_cross_section(plotaxes, 2)
    ml_plot.add_land_cross_section(plotaxes)


    # ========================================================================
    #  Water Speed
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='speed', figno=1)
    plotfigure.show = False
    plotfigure.kwargs = {'figsize':(14,4)}

    # Top layer speed
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Currents - Top Layer'
    plotaxes.scaled = True
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.axescmd = 'subplot(1,2,1)'
    ml_plot.add_speed(plotaxes,1,bounds=top_speed_limits)
    ml_plot.add_land(plotaxes, 1)
    
    # Bottom layer speed
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Currents - Bottom Layer'
    plotaxes.scaled = True
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.axescmd = 'subplot(1,2,2)'
    ml_plot.add_speed(plotaxes,2,bounds=internal_speed_limits)
    ml_plot.add_land(plotaxes, 2)
    
    # Individual components
    plotfigure = plotdata.new_plotfigure(name='speed_components',figno=401)
    plotfigure.show = False
    plotfigure.kwargs = {'figsize':(14,14)}
    
    # Top layer
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = "X-Velocity - Top Layer"
    plotaxes.scaled = True
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.axescmd = 'subplot(2,2,1)'
    ml_plot.add_x_velocity(plotaxes,1)
    ml_plot.add_land(plotaxes, 1)
    
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = "Y-Velocity - Top Layer"
    plotaxes.scaled = True
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.axescmd = 'subplot(2,2,2)'
    ml_plot.add_y_velocity(plotaxes,1)
    ml_plot.add_land(plotaxes, 1)
    
    # Bottom layer
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = "X-Velocity - Bottom Layer"
    plotaxes.scaled = True
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.axescmd = 'subplot(2,2,3)'
    ml_plot.add_x_velocity(plotaxes,2)
    ml_plot.add_land(plotaxes, 2)
    
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = "Y-Velocity - Bottom Layer"
    plotaxes.scaled = True
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.axescmd = 'subplot(2,2,4)'
    ml_plot.add_y_velocity(plotaxes,2)
    ml_plot.add_land(plotaxes, 2)


    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Surface at gauges', figno=300, \
                    type='each_gauge')
    plotfigure.clf_each_gauge = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Surface'

    # Plot surface as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 3
    plotitem.plotstyle = 'b-'

    # Plot topo as green curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.show = False

    def gaugetopo(current_data):
        q = current_data.q
        h = q[0,:]
        eta = q[3,:]
        topo = eta - h
        return topo
        
    plotitem.plot_var = gaugetopo
    plotitem.plotstyle = 'g-'

    def add_zeroline(current_data):
        from pylab import plot, legend, xticks, floor, axis, xlabel
        t = current_data.t 
        gaugeno = current_data.gaugeno

        if gaugeno == 32412:
            try:
                plot(TG32412[:,0], TG32412[:,1], 'r')
                legend(['GeoClaw','Obs'],loc='lower right')
            except: pass
            axis((0,t.max(),-0.3,0.3))

        plot(t, 0*t, 'k')
        n = int(floor(t.max()/3600.) + 2)
        xticks([3600*i for i in range(n)], ['%i' % i for i in range(n)])
        xlabel('time (hours)')


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
    plotdata.parallel = True                 # make multiple frame png's at once

    return plotdata

