
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

import numpy
from clawpack.geoclaw import dtopotools

dtopo1 = dtopotools.DTopography('dtopo1.tt3',3)
dtopo2 = dtopotools.DTopography('dtopo2.tt3',3)

#--------------------------
def setplot(plotdata):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 


    from clawpack.visclaw import colormaps, geoplot

    plotdata.clearfigures()  # clear any old figures,axes,items data

    # To plot gauge locations on pcolor or contour plot, use this as
    # an afteraxis function:

    def addgauges(current_data):
        from clawpack.visclaw import gaugetools
        gaugetools.plot_gauge_locations(current_data.plotdata, \
             gaugenos='all', format_string='ko', add_labels=True)


    #-----------------------------------------
    # Figure for pcolor plot
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='pcolor', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Surface'
    plotaxes.scaled = True

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.show = True
    plotitem.plot_var = geoplot.surface
    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmin = -5.
    plotitem.pcolor_cmax = 5.
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 1

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.land
    #plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmap = colormaps.yellow_red_blue
    plotitem.pcolor_cmin = 1.0
    plotitem.pcolor_cmax = 6.
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 1
    plotaxes.xlimits = [-1,1]
    plotaxes.ylimits = [-1,1]

    # Add contour lines of bathymetry:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.show = True
    plotitem.plot_var = geoplot.topo
    from numpy import arange, linspace
    #plotitem.contour_levels = [-14,-13,-12,-11]

    # need to adjust if h0 is changed in maketopo.py:
    plotitem.contour_levels = linspace(-3000,-1100,21)

    plotitem.amr_contour_colors = ['g']  # color on each level
    plotitem.kwargs = {'linestyles':'solid'}
    plotitem.amr_contour_show = [1,0]  
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0

    plotaxes.afteraxes = addgauges

    #-----------------------------------------
    # Figure for cross section
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='cross-section', figno=1)
    #plotfigure.show = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [-1,1]
    plotaxes.ylimits = [-5,7]
    y0 = 0.0
    plotaxes.title = 'Cross section at y=%s' % y0

    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')

    def xsec(current_data):
        # Return x value and surface eta at this point, along y=0
        from pylab import find,ravel
        x = current_data.x
        y = current_data.y
        dy = current_data.dy
        q = current_data.q

        ij = find((y <= y0+dy/2.+1e-8) & (y > y0-dy/2.))
        x_slice = ravel(x)[ij]
        eta_slice = ravel(q[3,:,:])[ij]
        return x_slice, eta_slice

    plotitem.map_2d_to_1d = xsec
    plotitem.plotstyle = 'o'
    plotitem.kwargs = {'markersize':5}
    plotitem.amr_show = [1]  # plot on all levels

    def add_dtopo_plot(current_data):
        from pylab import find, plot, legend
        for dtopo in [dtopo1,dtopo2]:
            j = max(find(dtopo.y <= y0))
            j1 = min(j+1, len(dtopo.y)-1)
            beta = (y0 - dtopo.y[j])/(dtopo.y[1] - dtopo.y[0])
    
            t = current_data.t
            itlist = find(dtopo.times <= t)
            if len(itlist) > 0:
                it = max(find(dtopo.times <= t))
                #print "+++ t, it, dtopo.times[it]: ",t, it, dtopo.times[it]
                dz = dtopo.dZ[it,:,:]
                dzj = dz[j,:]
                if it < len(dtopo.times)-1:
                    dz1 = dtopo.dZ[it+1,:,:]
                    dz1j = dz1[j,:]
                    alpha = (t - dtopo.times[it])/(dtopo.times[it+1] - dtopo.times[it])
                    #print "+++ alpha = ",alpha
                    dzj = (1.-alpha)*dzj + alpha*dz1j
                plot(dtopo.x, dzj, label="dtopo")
        legend()
    plotaxes.afteraxes = add_dtopo_plot


    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Surface at gauges', figno=300, \
                    type='each_gauge')
    plotfigure.clf_each_gauge = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(211)'
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Surface'

    # Plot surface as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 3
    plotitem.plotstyle = 'bo-'

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(212)'
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Topography'

    # Plot topo as green curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    #plotitem.show = False

    def gaugetopo(current_data):
        q = current_data.q
        h = q[0,:]
        eta = q[3,:]
        topo = eta - h
        return topo

    plotitem.plot_var = gaugetopo
    plotitem.plotstyle = 'g-'


    #-----------------------------------------
    
    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_gaugenos = []             # list of gauges to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?

    return plotdata

    
