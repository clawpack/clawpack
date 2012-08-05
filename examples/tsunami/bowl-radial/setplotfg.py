#!/usr/bin/env python
# encoding: utf-8
"""
Set up data for plotting fixed grids.
This is used in setplot.py.

To use interactively:
    >>> import setplotfg
    >>> fgdata = setplotfg.setplotfg(fgno, outdir)  # for grid number fgno 
    >>> fgdata.fgloop()          # to loop through frames "
    >>> fgdata.plotfg(frameno)   # to plot one frame"
    >>> fgdata.fg2html('all')    # to make html files of all plots"

"""

def setplotfg(fgno=1, outdir='_output'):
    from clawpack.visclaw import plotfg, colormaps, geoplot
    from numpy import arange

    # max used for surface and inundation colormaps below:
    etamax = 0.5

    fgdata = plotfg.ClawPlotFGData()

    # Set attributes as desired:

    fgdata.outdir = outdir
    fgdata.plotdir = '.'

    # Fixed grid to display:
    fgdata.fgno = fgno

    if fgno>0:

        # Could set things differently for each fgno if desired...

        # Contour levels for all plots:
        fgdata.clines = arange(-20,32,4)

        # For plot of surface eta each frame:
        fgdata.eta_show = True
        my_cmap_surface = colormaps.make_colormap({-1.0: [0.0,0.0,1.0], \
                                              -0.02: [0.75,0.75,1.0], \
                                               0.0: [1.0,1.0,1.0], \
                                               0.02: [1.0,0.75,0.75], \
                                               1.0: [1.0,0.0,0.0]})
        fgdata.water_cmap =  my_cmap_surface
        fgdata.water_clim = (-etamax,etamax)

        fgdata.land_cmap =  geoplot.land_colors
        fgdata.land_clim = (0,40)

        # For plot of inundation region:
        fgdata.inundated_show = True
        fgdata.inundated_clim =(0,etamax)
	# colormap that scales up to etamax:
        fgdata.inundated_cmap =  colormaps.make_colormap({0:[0,0.3,1],\
                0.1*etamax:[0,1,1],1.01*etamax:[0,1,1], \
                0.4*etamax:[0,1,0], etamax:[1,0.2,0.2]})
        fgdata.inundated_add_colorbar = True

        # For plot of exposed seafloor:
        fgdata.seafloor_show = True
        fgdata.seafloor_cmap =  geoplot.seafloor_colormap
        fgdata.seafloor_clim = (-1,0)

        fgdata.save_png = False
        fgdata.seafloor_add_colorbar = True

        fgdata.drytol = 1.e-2
        fgdata.exposed_tol = 1.e-2

    return fgdata
