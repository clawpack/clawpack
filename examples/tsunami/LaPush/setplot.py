
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

    # This file has been modifid significantly to illustrate Google Earth plotting
    # tools.

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

    # ---------------------------------------
    # Some KML info
    # ---------------------------------------
    plotdata.kml_index_fname = "LaPush"     # Name for .kmz and .kml files.
    plotdata.kml_name = "LaPush"
    plotdata.kml_publish = 'http://math.boisestate.edu/~calhoun/visclaw/GoogleEarth/kmz/'


    light_green_a = [0.0,1.0,1.0,1.0];
    transparent = [0.0,0.0,0.0,0.0]
    red_a = [1.0,0.0,0.0,1.0]
    lightblue = "#4E6498"  # Matches region fill color for comp. domain

    lightblue_cmap = colormaps.make_colormap({-1:light_green_a,
                                              0.0:lightblue,
                                              1:red_a})

    transparent_cmap = colormaps.make_colormap({-1:light_green_a,
                                                0.0:transparent,
                                                1:red_a})

    def cbar_transparent(filename):
        geoplot.kml_build_colorbar(filename,
                                   transparent_cmap,
                                   cmin_coast,cmax_coast)

    def cbar_lightblue(filename):
        geoplot.kml_build_colorbar(filename,
                                   lightblue_cmap,
                                   cmin_coast,cmax_coast)



    #-----------------------------------------
    # Figure for big area
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Pacific', figno=0)
    plotfigure.show = True

    plotfigure.use_for_kml = True
    plotfigure.kml_use_for_initial_view = True

    plotfigure.kml_xlimits = [-126.0, -124.0]
    plotfigure.kml_ylimits = [47.0,49.0]

    # Resolution : Use num_cells and refinement factors to figure out
    # correct resolution.
    # Refinement factors : [10,6,16]; maxlevel = 4
    rcl = 1   # resolution of coarsest level in this zoom;  rcl*figsize = num_cells
    plotfigure.kml_dpi = rcl*10*6   # Resolve levels 2,3
    plotfigure.kml_figsize = [12.0,12.0]
    plotfigure.kml_tile_images = True

    plotfigure.kml_colorbar = cbar_transparent

    # Plot sea surface height
    plotaxes = plotfigure.new_plotaxes()
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.pcolor_cmap = transparent_cmap
    plotitem.pcolor_cmin = cmin_ocean
    plotitem.pcolor_cmax = cmax_ocean

    #-----------------------------------------
    # Figure for zoom1
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name="Olympics", figno=1)
    plotfigure.show = True

    plotfigure.use_for_kml = True
    plotfigure.kml_use_for_initial_view = False

    # Limits should align with coarse grid.
    plotfigure.kml_xlimits = [-125.0, -124.0]
    plotfigure.kml_ylimits = [47.5,48.5]

    # rcl : Resolution of the coarsest level in this zoom, where
    #       rcl*figsize == number of coarsest level cells in each
    #       direction.
    rcl = 10
    plotfigure.kml_dpi = rcl*6    # Resolve level 3
    plotfigure.kml_figsize = [6.0,6.0]
    plotfigure.kml_tile_images = True

    plotfigure.kml_colorbar = cbar_transparent

    # Water
    plotaxes = plotfigure.new_plotaxes()
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.pcolor_cmap = transparent_cmap
    plotitem.pcolor_cmin = cmin_coast
    plotitem.pcolor_cmax = cmax_coast

    #-----------------------------------------
    # Figure for zoom2
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name="LaPush", figno=2)
    plotfigure.show = True

    plotfigure.use_for_kml = True
    plotfigure.kml_use_for_initial_view = False

    # Should align with a [720x720] grid in larger computational domain.
    plotfigure.kml_xlimits = [-124.7, -124.55]
    plotfigure.kml_ylimits = [47.8, 47.95]

    # Resolution
    # rcl : Resolution of the coarsest level in this zoom, where
    #       rcl*figsize == number of coarsest level cells in each
    #       direction.
    rcl = 10   # figsize*rcl = number of coarsest level cells in this zoom
    plotfigure.kml_dpi = rcl*16  # Resolve finest level.
    plotfigure.kml_figsize = [5.4,5.4]
    plotfigure.kml_tile_images = True

    plotfigure.kml_colorbar = cbar_lightblue


    # Water
    plotaxes = plotfigure.new_plotaxes()
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.pcolor_cmap = lightblue_cmap
    plotitem.pcolor_cmin = cmin_coast
    plotitem.pcolor_cmax = cmax_coast

    #-----------------------------------------

    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via clawpack.visclaw.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_fignos = [0,1,2]           # list of figures to print
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
