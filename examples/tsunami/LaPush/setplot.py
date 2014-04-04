
"""
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.

"""


import pylab
import glob, os
from numpy import loadtxt
from matplotlib import image


# Specific to WA coast project:
# try:
#     WAcoast = os.environ['WAcoast']
# except:
#     raise Exception("Need to set WAcoast environment variable")
#
# # Google Earth image for plotting on top of...
# # need to change for specific location on coast:
#
# plot_zeta_map = True  # set to false if no image available
# GEmap = image.imread(WAcoast + '/maps/LaPush.png')
# GEextent = (-124.645,-124.53,47.9,47.946)  # for LaPush
#

# --------------------------
def setplot(plotdata):
# --------------------------

    """
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of clawpack.visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.

    """


    from clawpack.visclaw import colormaps, geoplot

    plotdata.clearfigures()  # clear any old figures,axes,items dat
    plotdata.format = 'binary'

    clim_ocean = 8.0
    clim_coast = 8.0

    sealevel = 0.  # Level of tide in run relative to MHW
    cmax_ocean = clim_ocean + sealevel
    cmin_ocean = -clim_ocean + sealevel
    cmax_coast = clim_coast + sealevel
    cmin_coast = -clim_coast + sealevel


    # To plot gauge locations on pcolor or contour plot, use this as
    # an afteraxis function:

    def addgauges(current_data):
        from clawpack.visclaw import gaugetools
        gaugetools.plot_gauge_locations(current_data.plotdata, \
             gaugenos=[0,1,2,4,5,6,7], format_string='ko', add_labels=True)

    def timeformat(t):
        from numpy import mod
        hours = int(t/3600.)
        tmin = mod(t,3600.)
        min = int(tmin/60.)
        sec = int(mod(tmin,60.))
        timestr = '%s:%s:%s' % (hours,str(min).zfill(2),str(sec).zfill(2))
        return timestr

    def title_hours(current_data):
        from pylab import title
        t = current_data.t
        timestr = timeformat(t)
        title('%s after earthquake' % timestr)


    #-----------------------------------------
    # Figure for big area
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Pacific', figno=0)
    plotfigure.kwargs = {'figsize': (8,10)}
    #plotfigure.show = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Pacific'
    plotaxes.scaled = False
    plotaxes.xlimits = [-126, -124]
    plotaxes.ylimits = [47,49]

    def aa(current_data):
        from pylab import ticklabel_format, xticks, gca, cos, pi, savefig
        title_hours(current_data)
        ticklabel_format(format='plain',useOffset=False)
        xticks(rotation=20)
        a = gca()
        a.set_aspect(1./cos(48*pi/180.))
    plotaxes.afteraxes = aa

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = geoplot.surface_or_depth
    my_cmap = colormaps.make_colormap({-1.0: [0.0,0.0,1.0], \
                                     -0.5: [0.5,0.5,1.0], \
                                      0.0: [1.0,1.0,1.0], \
                                      0.5: [1.0,0.5,0.5], \
                                      1.0: [1.0,0.0,0.0]})
    plotitem.imshow_cmap = my_cmap
    plotitem.imshow_cmin = cmin_ocean
    plotitem.imshow_cmax = cmax_ocean
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [1]

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = geoplot.land
    plotitem.imshow_cmap = geoplot.land_colors
    plotitem.imshow_cmin = 0.0
    plotitem.imshow_cmax = 100.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [1]
    #plotaxes.afteraxes = addgauges

    # Add contour lines of bathymetry:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.show = False
    plotitem.plot_var = geoplot.topo
    from numpy import arange, linspace
    plotitem.contour_levels = linspace(-6000,0,7)
    plotitem.amr_contour_colors = ['g']  # color on each level
    plotitem.kwargs = {'linestyles':'solid'}
    plotitem.amr_contour_show = [0,0,1,0]  # show contours only on finest level
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0

    # Add contour lines of topography:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.show = False
    plotitem.plot_var = geoplot.topo
    from numpy import arange, linspace
    plotitem.contour_levels = arange(0., 11., 1.)
    plotitem.amr_contour_colors = ['g']  # color on each level
    plotitem.kwargs = {'linestyles':'solid'}
    plotitem.amr_contour_show = [0,0,0,1]  # show contours only on finest level
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0


    #-----------------------------------------
    # Figure for zoom1
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name="Olympics", figno=11)
    #plotfigure.show = False
    plotfigure.kwargs = {'figsize': (10,9)}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.scaled = False
    plotaxes.xlimits = [-125, -124.3]
    plotaxes.ylimits = [47.5, 48.5]
    plotaxes.afteraxes = aa

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.imshow_cmap = my_cmap
    plotitem.imshow_cmin = cmin_coast
    plotitem.imshow_cmax = cmax_coast
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [0]

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = geoplot.land
    plotitem.imshow_cmap = geoplot.land_colors
    plotitem.imshow_cmin = 0.0
    plotitem.imshow_cmax = 100.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [0]



    #-----------------------------------------
    # Figure for zoom2
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name="LaPush", figno=12)
    #plotfigure.show = True
    plotfigure.kwargs = {'figsize': (10,9)}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = "LaPush"
    plotaxes.scaled = False
    plotaxes.xlimits = [-124.68, -124.55]
    plotaxes.ylimits = [47.86, 47.96]
    plotaxes.afteraxes = aa

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.imshow_cmap = my_cmap
    plotitem.imshow_cmin = cmin_coast
    plotitem.imshow_cmax = cmax_coast
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [0]

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = geoplot.land
    plotitem.imshow_cmap = geoplot.land_colors
    plotitem.imshow_cmin = 0.0
    plotitem.imshow_cmax = 100.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [0]


    #-----------------------------------------
    # Figure for fixed grid area with GE image
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name="fgmax", figno=13)
    plotfigure.show = False
    plotfigure.kwargs = {'figsize': (12,8)}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = "LaPush"
    plotaxes.scaled = False
    plotaxes.xlimits = [-124.68, -124.55]
    plotaxes.ylimits = [47.86, 47.96]
    plotaxes.afteraxes = aa

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.imshow_cmap = my_cmap
    plotitem.imshow_cmin = cmin_coast
    plotitem.imshow_cmax = cmax_coast
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [0]

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = geoplot.land
    plotitem.imshow_cmap = geoplot.land_colors
    plotitem.imshow_cmin = 0.0
    plotitem.imshow_cmax = 100.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [0]

    #-------------------------------------------------------------------
    # Figure for KML files
    #--------------------------------------------------------------------
    plotfigure = plotdata.new_plotfigure(name='kml_figure',figno=1)
    plotfigure.show = True   # Don't show this file in the html version
    plotfigure.use_for_kml = True
    plotfigure.kml_dpi = 1600
    # plotfigure.kml_xlimits = [-124.68, -124.55]
    # plotfigure.kml_ylimits = [47.86, 47.96]
    plotfigure.kml_xlimits = [-126,-124]
    plotfigure.kml_ylimits = [47,49]


    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('kml')
    plotaxes.scaled = True

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.surface
    plotitem.pcolor_cmap = my_cmap
    plotitem.pcolor_cmin = cmin_coast
    plotitem.pcolor_cmax = cmax_coast
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 0

#     #-----------------------------------------
#     # Figures for gauges
#     #-----------------------------------------
#     plotfigure = plotdata.new_plotfigure(name='gauge plot', figno=300, \
#                     type='each_gauge')
#     #plotfigure.clf_each_gauge = False
#
#     # Set up for axes in this figure:
#     plotaxes = plotfigure.new_plotaxes()
#     plotaxes.ylimits = [-1,10]
#     plotaxes.title = 'Flow depth'
#
#     # Plot depth as blue curve:
#     plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
#     plotitem.plot_var = 0
#     plotitem.plotstyle = 'b-'
#


    #-----------------------------------------

    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via clawpack.visclaw.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.print_gaugenos = 'all'          # list of gauges to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?

    plotdata.format = 'binary'
    plotdata.kml = True

    return plotdata
