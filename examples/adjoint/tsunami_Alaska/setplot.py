
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.

This plotting function is for plotting when the adjoint-flagging method
is used.
The main difference is that the inner product plot is produced.
    
"""

from clawpack.geoclaw import topotools
import pylab
import glob
from numpy import loadtxt

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
    plotdata.format = 'binary'      # 'ascii', 'binary', 'netcdf'

    try:
        tsudata = open(plotdata.outdir+'/geoclaw.data').readlines()
        for line in tsudata:
            if 'sea_level' in line:
                sea_level = float(line.split()[0])
                print "sea_level = ",sea_level
    except:
        print "Could not read sea_level, setting to 0."
        sea_level = 0.

    clim_ocean = 0.3
    clim_CC = 0.5

    cmax_ocean = clim_ocean + sea_level
    cmin_ocean = -clim_ocean + sea_level
    cmax_CC = clim_CC + sea_level
    cmin_CC = -clim_CC + sea_level

    # To plot gauge locations on pcolor or contour plot, use this as
    # an afteraxis function:

    def addgauges(current_data):
        from clawpack.visclaw import gaugetools
        gaugetools.plot_gauge_locations(current_data.plotdata, \
             gaugenos='all', format_string='ko', add_labels=True)

    def timeformat(t):
        from numpy import mod
        hours = int(t/3600.)
        tmin = mod(t,3600.)
        min = int(tmin/60.)
        sec = int(mod(tmin,60.))
        timestr = '%s:%s:%s' % (hours,str(min).zfill(2),str(sec).zfill(2))
        return timestr
        
    def title(current_data):
        from pylab import title
        title('Surface Height', fontsize = 28)
    
    def title_innerproduct(current_data):
        from pylab import title
        title('Inner Product', fontsize = 28)

    def plotcc(current_data):
        from pylab import plot,text
        plot([235.8162], [41.745616],'wo')
        text(235.8,41.9,'Cr.City',color='w',fontsize=10)
    

    #-----------------------------------------
    # Figure for big area
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Pacific', figno=0)
    plotfigure.kwargs = {'figsize': (11,5)}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'axes([0.03,0.13,0.45,0.75])'
    plotaxes.title = 'Pacific'
    plotaxes.scaled = True

    def aa(current_data):
        from pylab import ticklabel_format, xticks, gca, cos, pi, yticks
        plotcc(current_data)
        title(current_data)
        ticklabel_format(format='plain',useOffset=False)
        xticks([180, 200, 220, 240], rotation=20, fontsize = 28)
        yticks(fontsize = 28)
        a = gca()
        a.set_aspect(1./cos(41.75*pi/180.))
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
    plotitem.add_colorbar = False
    plotitem.colorbar_shrink = 0.7
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
    plotitem.amr_patchedges_show = [0]
    plotaxes.xlimits = [180,240] 
    plotaxes.ylimits = [10,62]

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

    # Adding inner product plot

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Inner Product'
    plotaxes.axescmd = 'axes([0.45,0.13,0.6,0.75])'
    plotaxes.scaled = True
    def aa_innerprod(current_data):
        from pylab import ticklabel_format, xticks, gca, cos, pi, yticks
        plotcc(current_data)
        title_innerproduct(current_data)
        ticklabel_format(format='plain',useOffset=False)
        xticks([180, 200, 220, 240], rotation=20, fontsize = 28)
        pylab.tick_params(axis='y', labelleft='off')
        a = gca()
        a.set_aspect(1./cos(41.75*pi/180.))
    plotaxes.afteraxes = aa_innerprod
    
    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = 4
    plotitem.imshow_cmap = colormaps.white_red
    plotitem.imshow_cmin = 0.0
    plotitem.imshow_cmax = 0.05
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [0]
    plotitem.amr_data_show = [1,0,0,0]
    
    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = geoplot.land
    plotitem.imshow_cmap = geoplot.land_colors
    plotitem.imshow_cmin = 0.0
    plotitem.imshow_cmax = 100.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [0]
    plotaxes.xlimits = [180,240]
    plotaxes.ylimits = [10,62]


    #-----------------------------------------
    # Figure for zoom
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='California', figno=10)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'California'
    plotaxes.scaled = True
    plotaxes.afteraxes = aa

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.imshow_cmap = my_cmap
    plotitem.imshow_cmin = cmin_ocean
    plotitem.imshow_cmax = cmax_ocean
    plotitem.add_colorbar = False
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
    plotaxes.xlimits = [235.5,236]
    plotaxes.ylimits = [41.6,41.8]
    plotaxes.afteraxes = addgauges

    # add contour lines of bathy if desired:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.plot_var = geoplot.topo
    plotitem.contour_levels = [0.]
    plotitem.amr_contour_colors = ['k']  # color on each level
    plotitem.kwargs = {'linestyles':'solid','linewidths':2}
    plotitem.amr_contour_show = [0,0,0,0,0,1]
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0

    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    
    def fix_gauge(current_data):
        from pylab import plot, legend, xticks, floor, yticks,xlabel,savefig, title
        t = current_data.t
        gaugeno = current_data.gaugeno
        xticks([18000, 21600, 25200, 28800, 32400, 36000],\
               [str(180/36), str(216/36), str(252/36), str(288/36), \
                str(324/36), str(360/36)], fontsize=17)
        yticks(fontsize=17)
        title('Surface at gauge ' + str(gaugeno), fontsize=17)
        xlabel(' ')

    plotfigure = plotdata.new_plotfigure(name='gauge plot', figno=300, \
                    type='each_gauge')
    plotfigure.kwargs = {'figsize': (10.5,5)}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'axes([0.12,0.12,0.79,0.79])'
    plotaxes.title = 'Surface'
    plotaxes.xlimits = [15000, 39600]
    plotaxes.afteraxes = fix_gauge

    # Plot surface as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 3
    plotitem.plotstyle = 'b-'

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
    plotdata.html_movie = 'JSAnimation'      # new style, or "4.x" for old style


    return plotdata

    
