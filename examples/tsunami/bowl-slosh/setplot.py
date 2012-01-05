
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

from numpy import sqrt
a = 1.
sigma = 0.5
h0 = 0.1
grav = 9.81
omega = sqrt(2.*grav*h0) / a 


#--------------------------
def setplot(plotdata):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 


    from pyclaw.plotters import colormaps, geoplot

    plotdata.clearfigures()  # clear any old figures,axes,items data

    def set_drytol(current_data):
        # The drytol parameter is used in masking land and water and
        # affects what color map is used for cells with small water depth h.
        # The cell will be plotted as dry if h < drytol.
        # The best value to use often depends on the application and can
        # be set here (measured in meters):
        current_data.user.drytol = 1.e-3

    plotdata.beforeframe = set_drytol

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
    plotitem.plot_var = geoplot.surface
    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmin = -0.1
    plotitem.pcolor_cmax = 0.1
    plotitem.add_colorbar = True
    plotitem.amr_gridlines_show = [0,0,0]
    plotitem.gridedges_show = 1

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 100.0
    plotitem.add_colorbar = False
    plotitem.amr_gridlines_show = [0,0,0]
    plotitem.gridedges_show = 1
    plotaxes.xlimits = [-2,2]
    plotaxes.ylimits = [-2,2]

    # Add contour lines of bathymetry:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.plot_var = geoplot.topo
    from numpy import arange, linspace
    plotitem.contour_levels = linspace(-.1, 0.5, 20)
    plotitem.amr_contour_colors = ['k']  # color on each level
    plotitem.kwargs = {'linestyles':'solid'}
    plotitem.amr_contour_show = [1]  
    plotitem.gridlines_show = 0
    plotitem.gridedges_show = 0
    plotitem.show = True

    #-----------------------------------------
    # Figure for cross section
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='cross-section', figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [-2,2]
    plotaxes.ylimits = [-0.15,0.3]
    plotaxes.title = 'Cross section at y=0'
    def plot_topo_xsec(current_data):
        from pylab import plot, hold, cos,sin,where,legend,nan
        x = current_data.x
        y = current_data.y
        t = current_data.t
        q = current_data.q
        my2 = q.shape[1] / 2.
        x = x[:,my2]
        y = y[:,my2]
        hold(True)
        #topo = q[:,my2,3] - q[:,my2,0]  # should equal B below
        #plot(x, topo, 'g')
        B = h0*(x**2 + y**2)/a**2 - h0
        eta1 = sigma*h0/a**2 * (2.*x*cos(omega*t) + 2.*y*sin(omega*t) -sigma)
        etatrue = where(eta1>B, eta1, nan)
        plot(x, etatrue, 'r')
        legend(['computed','true'])
        plot(x, B, 'g')
        hold(False)
    plotaxes.afteraxes = plot_topo_xsec

    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')

    def xsec(current_data):
        # Return x value and surface eta at this point, along y=0
        x = current_data.x
        q = current_data.q
        my2 = q.shape[1] / 2.
        eta_slice = q[:,my2,3]
        return x[:,my2], eta_slice

    plotitem.map_2d_to_1d = xsec
    plotitem.plotstyle = 'b-'
    plotitem.amr_plot_show = [0,1]  # only plot on level 2


    #-----------------------------------------
    # Figure for grids alone
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='grids', figno=2)
    plotfigure.show = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [-2,2]
    plotaxes.ylimits = [-2,2]
    plotaxes.title = 'grids'
    plotaxes.scaled = True

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_grid')
    plotitem.amr_grid_bgcolor = ['#ffeeee', '#eeeeff', '#eeffee']
    plotitem.amr_gridlines_show = [1,1,0]   
    plotitem.amr_gridedges_show = [1]     



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

    
