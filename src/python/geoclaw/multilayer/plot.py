#!/usr/bin/env python
# encoding: utf-8
r"""
Plotting routines for multilayer simulations with GeoClaw

:Authors:
    Kyle Mandli (2011-2-07) Initial version
"""

from __future__ import absolute_import
import numpy
import matplotlib.pyplot as plt

from clawpack.visclaw import colormaps, geoplot, gaugetools
from six.moves import range

# Definition of colormpas
surface_cmap = plt.get_cmap("bwr")
speed_cmap = plt.get_cmap('PuBu')
land_cmap = geoplot.land_colors
velocity_cmap = plt.get_cmap('PiYG')
vorticity_cmap = plt.get_cmap('PRGn')

# ========================================================================
#  Data extraction routines
#  This array contains the location of the start of each layer's data.  To
#  access a layer's depth for instance you want to use
#    q[layer_index[0] + 1, :, :]
#  Momenta in the x and y directions are `+ 2` and `+ 3` respectively.
#
#  Note also that all layer indices are 0-indexed
layer_index = [0, 3]
eta_index = 6


def extract_eta(h, eta, DRY_TOL=1e-3):
    masked_eta = numpy.ma.masked_where(numpy.abs(h) < DRY_TOL, eta)
    return masked_eta


def eta(cd, layer):
    return extract_eta(cd.q[layer_index[layer], :, :],
                       cd.q[eta_index + layer, :, :])


def eta1(cd):
    return eta(cd, 0)


def eta2(cd):
    return eta(cd, 1)


def b(current_data):
    h1 = current_data.q[layer_index[0], :, :]
    h2 = current_data.q[layer_index[1], :, :]

    return eta1(cd) - h1 - h2


def extract_velocity(h, hu, DRY_TOL=10**-8):
    u = numpy.ones(hu.shape) * numpy.nan
    index = numpy.nonzero((numpy.abs(h) > DRY_TOL) * (h != numpy.nan))
    u[index[0], index[1]] = hu[index[0], index[1]] / h[index[0], index[1]]
    return u


def water_u(cd, direction, layer):
    return extract_velocity(cd.q[layer_index[layer], :, :],
                            cd.q[layer_index[layer] + direction + 1, :, :])


def water_u1(cd):
    return water_u(cd, 0, 0)


def water_u2(cd):
    return water_u(cd, 0, 1)


def water_v1(cd):
    return water_u(cd, 1, 0)


def water_v2(cd):
    return water_u(cd, 1, 1)


def water_speed1(current_data):
    u = water_u1(current_data)
    v = water_v1(current_data)

    return numpy.sqrt(u**2 + v**2)


def water_speed2(current_data):
    u = water_u2(current_data)
    v = water_v2(current_data)

    return numpy.sqrt(u**2 + v**2)


def water_speed_depth_ave(current_data):
    h1 = current_data.q[layer_index[0], :, :]
    h2 = current_data.q[layer_index[1], :, :]
    u1 = water_speed1(current_data)
    u2 = water_speed1(current_data)

    return (h1 * u1 + h2 * u2) / (h1 + h2)


# ========================================================================
#  Plot items
# ========================================================================
def add_surface_elevation(plotaxes, surface, bounds=None, plot_type='pcolor'):
    if plot_type == 'pcolor' or plot_type == 'imshow':
        plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
        # plotitem.plot_var = geoplot.surface
        if surface == 1:
            plotitem.plot_var = eta1
        elif surface == 2:
            plotitem.plot_var = eta2
        if bounds is not None:
            plotitem.pcolor_cmin = bounds[0]
            plotitem.pcolor_cmax = bounds[1]
        plotitem.pcolor_cmap = surface_cmap
        plotitem.add_colorbar = True
        plotitem.amr_celledges_show = [0] * 10
        plotitem.amr_patchedges_show = [0] * 10

    elif plot_type == 'contour':
        plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
        plotitem.plot_var = surface + 5
        if bounds is not None:
            plotitem.contour_levels = bounds
        plotitem.amr_contour_show = [1] * 10
        plotitem.amr_celledges_show = [0] * 10
        plotitem.amr_patchedges_show = [1] * 10
        plotitem.amr_contour_colors = 'k'

    else:
        raise NotImplementedError("Plot type %s not implemented" % plot_type)


def add_layer_depth(plotaxes, layer, bounds=None, plot_type='pcolor'):

    if plot_type == 'pcolor' or plot_type == 'imshow':
        plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
        if layer == 1:
            plotitem.plot_var = 0
        elif layer == 2:
            plotitem.plot_var = 3
        if bounds is not None:
            plotitem.imshow_cmin = bounds[0]
            plotitem.imshow_cmax = bounds[1]
        plotitem.imshow_cmap = surface_cmap
        plotitem.add_colorbar = True
        plotitem.amr_celledges_show = [0] * 10
        plotitem.amr_patchedges_show = [1] * 10


def add_speed(plotaxes, layer, bounds=None, plot_type='pcolor'):

    if plot_type == 'pcolor' or plot_type == 'imshow':
        plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
        if layer == 1:
            plotitem.plot_var = water_speed1
        elif layer == 2:
            plotitem.plot_var = water_speed2
        plotitem.imshow_cmap = speed_cmap
        if bounds is not None:
            plotitem.imshow_cmin = bounds[0]
            plotitem.imshow_cmax = bounds[1]
        plotitem.add_colorbar = True
        plotitem.amr_celledges_show = [0] * 10
        plotitem.amr_patchedges_show = [1] * 10


def add_x_velocity(plotaxes, layer, plot_type='pcolor', bounds=None):

    if plot_type == 'pcolor' or plot_type == 'imshow':
        plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
        if layer == 1:
            plotitem.plot_var = water_u1
        if layer == 2:
            plotitem.plot_var = water_u2
        if bounds is not None:
            plotitem.imshow_cmin = bounds[0]
            plotitem.imshow_cmax = bounds[1]
        plotitem.add_colorbar = True
        plotitem.imshow_cmap = velocity_cmap
        plotitem.amr_celledges_show = [0] * 10
        plotitem.amr_patchedges_show = [1] * 10


def add_y_velocity(plotaxes, layer, plot_type='pcolor', bounds=None):

    if plot_type == 'pcolor' or plot_type == 'imshow':
        plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
        if layer == 1:
            plotitem.plot_var = water_v1
        if layer == 2:
            plotitem.plot_var = water_v2
        if bounds is not None:
            plotitem.imshow_cmin = bounds[0]
            plotitem.imshow_cmax = bounds[1]
        plotitem.imshow_cmap = velocity_cmap
        plotitem.add_colorbar = True
        plotitem.amr_celledges_show = [0] * 10
        plotitem.amr_patchedges_show = [1] * 10


# Land
def add_land(plotaxes, plot_type='pcolor', bounds=[-10, 10]):
    r"""Add plot item for land"""

    if plot_type == 'pcolor':
        plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
        plotitem.show = True
        plotitem.plot_var = geoplot.land
        plotitem.imshow_cmap = land_cmap
        plotitem.imshow_cmin = bounds[0]
        plotitem.imshow_cmax = bounds[1]
        plotitem.add_colorbar = False
        plotitem.amr_celledges_show = [0] * 10
        plotitem.amr_patchedges_show = [1] * 10

    elif plot_type == 'contour':
        plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
        plotitem.plot_var = geoplot.land
        plotitem.contour_nlevels = 40
        plotitem.contour_min = bounds[0]
        plotitem.contour_max = bounds[1]
        plotitem.amr_contour_colors = ['g']  # color on each level
        plotitem.amr_patch_bgcolor = ['#ffeeee', '#eeeeff', '#eeffee']
        plotitem.celledges_show = 0
        plotitem.patchedges_show = 0
        plotitem.show = True

    else:
        raise NotImplementedError("Plot type %s not implemented" % plot_type)


def add_combined_profile_plot(plot_data, slice_value, direction='x',
                              figno=120):

    def slice_index(cd):
        if direction == 'x':
            if cd.grid.y.lower < slice_value < cd.grid.y.upper:
                return int((slice_value - cd.grid.y.lower) / cd.dy - 0.5)
            else:
                return None
        elif direction == 'y':
            if cd.grid.x.lower < slice_value < cd.grid.x.upper:
                return int((slice_value - cd.grid.x.lower) / cd.dx - 0.5)
            else:
                return None

    def bathy_profile(current_data):
        index = slice_index(current_data)
        if direction == 'x':
            if index:
                return current_data.x[:, index], b(current_data)[:, index]
            else:
                return None
        elif direction == 'y':
            if index:
                return current_data.y[index, :], b(current_data)[index, :]
            else:
                return None

    def lower_surface(current_data):
        index = slice_index(current_data)
        if direction == 'x':
            if index:
                return current_data.x[:, index], eta2(current_data)[:, index]
            else:
                return None
        elif direction == 'y':
            if index:
                return current_data.y[index, :], eta2(current_data)[index, :]
            else:
                return None

    def upper_surface(current_data):
        index = slice_index(current_data)
        if direction == 'x':
            if index:
                return current_data.x[:, index], eta1(current_data)[:, index]
            else:
                return None
        elif direction == 'y':
            if index:
                return current_data.y[index, :], eta1(current_data)[index, :]
            else:
                return None

    # Surfaces
    plotfigure = plotdata.new_plotfigure(name='combined_surface_%s' % figno,
                                         figno=figno)
    plotfigure.show = True
    plotfigure.kwargs = {'figsize': (6, 6)}

    # Top surface
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(2,1,1)'
    plotaxes.title = 'Surfaces Profile %s at %s' % (direction, slice_value)
    if multilayer_data.init_type == 2:
        plotaxes.xlimits = xlimits
    elif multilayer_data.init_type == 6:
        plotaxes.xlimits = ylimits

    plotaxes.ylimits = top_surf_zoomed

    def top_surf_afteraxes(cd):
        axes = plt.gca()
        axes.set_xlabel('')
        locs, labels = plt.xticks()
        labels = ['' for i in range(len(locs))]
        plt.xticks(locs, labels)
        axes.plot([multilayer_data.bathy_location,
                  multilayer_data.bathy_location], top_surf_zoomed, '--k')
        axes.set_ylabel('m')

    plotaxes.afteraxes = top_surf_afteraxes
    plotitem = plotaxes.new_plotitem(plot_type="1d_from_2d_data")
    plotitem.map_2d_to_1d = upper_surface
    plotitem.amr_plotstyle = ['-', '+', 'x']
    plotitem.show = True

    # Internal surface
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(2,1,2)'
    plotaxes.title = ''
    if multilayer_data.init_type == 2:
        plotaxes.xlimits = xlimits
    elif multilayer_data.init_type == 6:
        plotaxes.xlimits = ylimits
    plotaxes.ylimits = bottom_surf_zoomed

    def internal_surf_afteraxes(cd):
        axes = plt.gca()
        axes.set_title('')
        axes.set_ylabel('m')
        axes.subplots_adjust(hspace=0.05)
        axes.plot([multilayer_data.bathy_location,
                  multilayer_data.bathy_location], bottom_surf_zoomed, '--k')
    plotaxes.afteraxes = internal_surf_afteraxes
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.map_2d_to_1d = lower_surface
    plotitem.amr_plotstyle = ['-', '+', 'x']
    plotitem.color = 'k'
    plotitem.show = True


def add_velocities_profile_plot(plot_data, slice_value, direction='x',
                                figno=130):

    def slice_index(cd):
        if cd.grid.y.lower < slice_value < cd.grid.y.upper:
            return int((slice_value - cd.grid.y.lower) / cd.dy - 0.5)
        else:
            return None

    def upper_surface(current_data):
        index = slice_index(current_data)
        if index:
            return current_data.x[:, index], eta1(current_data)[:, index]
        else:
            return None

    def top_speed(current_data):
        index = slice_index(current_data)
        if index:
            return current_data.x[:, index], water_u1(current_data)[:, index]
        else:
            return None, None

    def bottom_speed(current_data):
        index = slice_index(current_data)
        if index:
            return current_data.x[:, index], water_u2(current_data)[:, index]
        else:
            return None, None

    # Velocities
    plotfigure = plotdata.new_plotfigure(name='combined_velocities_%s'
                                         % figno, figno=figno)
    plotfigure.show = True

    # Top surface
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Velocities Profile %s at %s' % (direction, slice_value)
    if multilayer_data.init_type == 2:
        plotaxes.xlimits = xlimits
    elif multilayer_data.init_type == 6:
        plotaxes.xlimits = ylimits
    plotaxes.ylimits = velocities_zoomed

    def velocity_afteraxes(cd):
        axes = plt.gca()
        axes.set_xlabel('')
        locs, labels = plt.xticks()
        labels = ['' for i in range(len(locs))]
        plt.xticks(locs, labels)
        axes.plot([multilayer_data.bathy_location,
                  multilayer_data.bathy_location], velocities_zoomed, '--k')
        axes.set_ylabel('m/s')
    plotaxes.afteraxes = velocity_afteraxes

    plotitem = plotaxes.new_plotitem(plot_type="1d_from_2d_data")
    plotitem.map_2d_to_1d = top_speed
    plotitem.amr_plotstyle = ['-', '+', 'x']
    plotitem.show = True

    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.map_2d_to_1d = bottom_speed
    plotitem.amr_plotstyle = ['-', '+', 'x']
    plotitem.color = 'k'
    plotitem.show = True
