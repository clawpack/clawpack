# encoding: utf-8
r"""
Plotting routines for storm surge simulations with GeoClaw

:Authors:
    Kyle Mandli (2012-10-23) Initial version
"""
# ==============================================================================
#  Distributed under the terms of the Berkeley Software Distribution (BSD)
#  license
#                     http://www.opensource.org/licenses/
# ==============================================================================

from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

import clawpack.visclaw.colormaps as colormaps
import clawpack.visclaw.gaugetools as gaugetools
import clawpack.geoclaw.geoplot as geoplot

# TODO:  Assign these absed on data files
bathy_index = 0
friction_field = 3
wind_field = 4
pressure_field = 6

surface_cmap = plt.get_cmap("bwr")
speed_cmap = plt.get_cmap('PuBu')
friction_cmap = plt.get_cmap('YlOrRd')
velocity_cmap = plt.get_cmap('PiYG')
vorticity_cmap = plt.get_cmap('PRGn')
wind_cmap = plt.get_cmap('PuBu')
pressure_cmap = plt.get_cmap('PuBu')
land_cmap = geoplot.land_colors


class track_data(object):
    """Read in storm track data from run output"""

    def __init__(self, path=None):
        if path is None:
            path = "fort.track"

        try:
            self._path = path
            self._data = np.loadtxt(self._path)
        except:
            self._data = None

    def get_track(self, frame):
        """Return storm location for frame requested"""

        # If data was not load successfully return None
        if self._data is None or len(self._data.shape) < 2:
            return None, None, None

        # If it appears that our data is not long enough, try reloading file
        if self._data.shape[0] < frame + 1:
            self._data = np.loadtxt(self._path)

            # Check to make sure that this fixed the problem
            if self._data.shape[0] < frame + 1:
                print(" *** WARNING *** Could not find track data for ",
                      "frame %s." % frame)
                return None, None, None

        return self._data[frame, 1:]


# ==========================================================================
# Gauge functions
# ==========================================================================
def gauge_locations(current_data, gaugenos='all'):
    gaugetools.plot_gauge_locations(current_data.plotdata,
                                    gaugenos=gaugenos, format_string='kx',
                                    add_labels=True, xoffset=0.02,
                                    yoffset=0.02)


def gaugetopo(current_data):
    q = current_data.q
    h = q[0, :]
    eta = q[3, :]
    topo = eta - h
    return topo


def plot_landfall_gauge(gauge, axes, landfall=0.0, style='b', kwargs={}):
    """Plot gauge data on the axes provided

    This will transform the plot so that it is relative to the landfall value
    provided.
    """
    axes = plt.gca()

    # Add GeoClaw gauge data
    t = sec2days(gauge.t - landfall)
    axes.plot(t, gauge.q[3, :], style, **kwargs)


# ========================================================================
#  Surge related helper functions
# ========================================================================
def days_figure_title(current_data, land_fall=0.0):
    t = (current_data.t - land_fall) / (60**2 * 24)
    days = int(t)
    hours = (t - int(t)) * 24.0

    title = current_data.plotaxes.title
    plt.title('%s at day %3i, hour %2.1f' % (title, days, hours))


def surge_afteraxes(current_data, track, land_fall=0.0, plot_direction=False,
                    style='ro', kwargs={}):
    """Default surge plotting after axes function

    Includes changing the title to something relative to landfall and plotting
    the location of the storm eye according to the track object.
    """
    track_data = track.get_track(current_data.frameno)

    if track_data[0] is not None and track_data[1] is not None:
        axes = plt.gca()
        axes.plot(track_data[0], track_data[1], style, **kwargs)
        if plot_direction:
            axes.quiver(track_data[0], track_data[1],
                        np.cos(track_data[2]), np.sin(track_data[2]))
    days_figure_title(current_data, land_fall)


def friction(cd):
    return cd.aux[friction_field, :, :]


def wind_x(cd):
    return cd.aux[wind_field, :, :]


def wind_y(cd):
    return cd.aux[wind_field+1, :, :]


def wind_speed(cd):
    return np.sqrt(wind_x(cd)**2 + wind_y(cd)**2)


def pressure(cd):
    # The division by 100.0 is to convert from Pa to millibars
    return cd.aux[pressure_field, :, :] / 100.0


# ========================================================================
#  Water helper functions
# ========================================================================
def b(cd):
    return cd.aux[bathy_index, :, :]


def extract_eta(h, eta, DRY_TOL=1e-3):
    index = np.nonzero((np.abs(h) < DRY_TOL) + (h == np.nan))
    eta[index[0], index[1]] = np.nan
    return eta


def extract_velocity(h, hu, DRY_TOL=1e-8):
    u = np.zeros(hu.shape)
    index = np.nonzero((np.abs(h) > DRY_TOL) * (h != np.nan))
    u[index[0], index[1]] = hu[index[0], index[1]] / h[index[0], index[1]]
    return u


def eta(cd):
    return extract_eta(cd.q[0, :, :], cd.q[3, :, :])


def water_u(cd):
    return extract_velocity(cd.q[0, :, :], cd.q[1, :, :])


def water_v(cd):
    return extract_velocity(cd.q[0, :, :], cd.q[2, :, :])


def water_speed(current_data):
    u = water_u(current_data)
    v = water_v(current_data)

    return np.sqrt(u**2+v**2)


# ========================================================================
#  Plot items
# ========================================================================
def add_surface_elevation(plotaxes, plot_type='pcolor', bounds=None,
                          contours=None, shrink=1.0):
    """Add plotitem representing the sea surface."""

    if plot_type == 'pcolor' or plot_type == 'imshow':
        plotitem = plotaxes.new_plotitem(name='surface', plot_type='2d_pcolor')
        plotitem.plot_var = geoplot.surface_or_depth

        if bounds is not None:
            if bounds[0] == 0.0:
                plotitem.pcolor_cmap = plt.get_cmap('OrRd')
            else:
                plotitem.pcolor_cmap = surface_cmap
            plotitem.pcolor_cmin = bounds[0]
            plotitem.pcolor_cmax = bounds[1]
        plotitem.add_colorbar = True
        plotitem.colorbar_shrink = shrink
        plotitem.colorbar_label = "Surface Height (m)"
        plotitem.amr_celledges_show = [0] * 10
        plotitem.amr_patchedges_show = [1, 1, 1, 0, 0, 0, 0]

    elif plot_type == 'contour':
        plotitem = plotaxes.new_plotitem(name='surface',
                                         plot_type='2d_contour')
        if bounds is None:
            plotitem.contour_levels = [-2.5, -1.5, -0.5, 0.5, 1.5, 2.5]

        plotitem.plot_var = geoplot.surface_or_depth
        plotitem.amr_contour_show = [1] * 10
        plotitem.amr_celledges_show = [0] * 10
        plotitem.amr_patchedges_show = [1, 1, 1, 1, 0, 0, 0]
        plotitem.amr_contour_colors = 'k'

    elif plot_type == 'contourf':
        plotitem = plotaxes.new_plotitem(name='surface',
                                         plot_type='2d_contourf')
        plotitem.plot_var = geoplot.surface_or_depth
        if bounds is not None:
            contours = numpy.linspace(bounds[0], bounds[1], 11)
            plotitem.contour_levels = contours
            plotitem.fill_cmin = bounds[0]
            plotitem.fill_cmax = bounds[1]
        elif contours is not None:
            plotitem.contour_levels = contours
            plotitem.fill_cmin = min(contours)
            plotitem.fill_cmax = max(contours)

        plotitem.add_colorbar = True
        plotitem.fill_cmap = geoplot.tsunami_colormap
        plotitem.colorbar_shrink = shrink
        plotitem.colorbar_label = "Surface Height (m)"
        plotitem.fill_cmap = plt.get_cmap('OrRd')
        if any((value < 0 for value in plotitem.contour_levels)):
            plotitem.fill_cmap = surface_cmap

        plotitem.amr_contour_show = [1] * 10
        plotitem.amr_celledges_show = [0] * 10
        plotitem.amr_patchedges_show = [1, 1, 1, 1, 0, 0, 0]
        plotitem.amr_contour_colors = 'k'


def add_speed(plotaxes, plot_type='pcolor', bounds=None,  contours=None,
              shrink=1.0):
    """Add plotitem representing speed of the water."""
    if plot_type == 'pcolor' or plot_type == 'imshow':
        plotitem = plotaxes.new_plotitem(name='speed', plot_type='2d_pcolor')
        plotitem.plot_var = water_speed
        # plotitem.plot_var = 1
        plotitem.pcolor_cmap = speed_cmap
        if bounds is not None:
            plotitem.pcolor_cmin = bounds[0]
            plotitem.pcolor_cmax = bounds[1]
        plotitem.add_colorbar = True
        plotitem.colorbar_shrink = shrink
        plotitem.colorbar_label = "Current (m/s)"
        plotitem.amr_celledges_show = [0] * 10
        plotitem.amr_patchedges_show = [1, 1, 1, 0, 0, 0, 0]

    elif plot_type == 'contour':
        plotitem = plotaxes.new_plotitem(name='speed', plot_type='2d_contour')
        if bounds is None:
            plotitem.contour_levels = [0.5, 1.5, 3, 4.5, 6.0]
        plotitem.kwargs = {'linewidths': 1}

        plotitem.plot_var = water_speed
        plotitem.amr_contour_show = [1] * 10
        plotitem.amr_celledges_show = [0] * 10
        plotitem.amr_patchedges_show = [1, 1, 1, 1, 1, 0, 0]
        plotitem.amr_contour_colors = 'k'

    elif plot_type == 'contourf':

        plotitem = plotaxes.new_plotitem(name='speed', plot_type='2d_contourf')

        plotitem.add_colorbar = True
        plotitem.colorbar_label = "Current (m/s)"
        plotitem.colorbar_shrink = shrink
        plotitem.fill_cmap = plt.get_cmap('PuBu')
        if bounds is not None:
            plotitem.contour_levels = numpy.linspace(bounds[0], bounds[1], 11)
            plotitem.fill_cmin = bounds[0]
            plotitem.fill_cmap = bounds[1]
        elif contours is not None:
            plotitem.contour_levels = contours
            plotitem.fill_cmin = min(contours)
            plotitem.fill_cmax = max(contours)

        # Modify the 'extends' plot attribute as we don't want this to extend
        # below 0
        plotitem.kwargs['extend'] = 'max'

        plotitem.plot_var = water_speed
        plotitem.amr_contour_show = [1] * 10
        plotitem.amr_celledges_show = [0] * 10
        plotitem.amr_patchedges_show = [1, 1, 1, 1, 1, 0, 0]
        plotitem.amr_contour_colors = 'k'


def add_friction(plotaxes, bounds=None, plot_type='pcolor', shrink=1.0):
    """Add plotitem for the friction field"""

    if plot_type == 'pcolor' or plot_type == 'imshow':
        plotitem = plotaxes.new_plotitem(name='friction',
                                         plot_type='2d_pcolor')
        plotitem.plot_var = friction
        plotitem.pcolor_cmap = friction_cmap
        plotitem.colorbar_shrink = shrink
        if bounds is not None:
            plotitem.pcolor_cmin = bounds[0]
            plotitem.pcolor_cmax = bounds[1]
        plotitem.add_colorbar = True
        plotitem.colorbar_shrink = shrink
        plotitem.colorbar_label = "Manning's-$n$ Coefficient"
        plotitem.amr_celledges_show = [0] * 10
        plotitem.amr_patchedges_show = [0] * 10


def add_wind(plotaxes, bounds=None, plot_type='pcolor', shrink=1.0):
    """Add plotitem for the wind speed."""

    if plot_type == 'pcolor' or plot_type == 'imshow':
        plotitem = plotaxes.new_plotitem(name='wind', plot_type='2d_pcolor')
        plotitem.plot_var = wind_speed
        plotitem.pcolor_cmap = wind_cmap
        if bounds is not None:
            plotitem.pcolor_cmin = bounds[0]
            plotitem.pcolor_cmax = bounds[1]
        plotitem.add_colorbar = True
        plotitem.colorbar_shrink = shrink
        plotitem.colorbar_label = "Wind Speed (m/s)"
        plotitem.amr_celledges_show = [0] * 10
        plotitem.amr_patchedges_show = [1, 1, 1, 1, 1, 0, 0]

    elif plot_type == 'contour':
        plotitem = plotaxes.new_plotitem(name='wind', plot_type='2d_contour')
        plotitem.plot_var = wind_speed
        plotitem.contour_nlevels = len(surge_data.wind_refine)
        plotitem.countour_min = surge_data.wind_refine[0]
        plotitem.patchedges_show = 1


def add_pressure(plotaxes, bounds=None, plot_type='pcolor', shrink=1.0):
    """Add plotitem for the pressure field."""

    if plot_type == 'pcolor' or plot_type == 'imshow':
        plotitem = plotaxes.new_plotitem(name="pressure",
                                         plot_type='2d_pcolor')
        plotitem.plot_var = pressure
        plotitem.colorbar_shrink = shrink
        plotitem.pcolor_cmap = pressure_cmap
        if bounds is not None:
            plotitem.pcolor_cmin = bounds[0]
            plotitem.pcolor_cmax = bounds[1]
        plotitem.add_colorbar = True
        plotitem.colorbar_shrink = shrink
        plotitem.colorbar_label = "Pressure (mbar)"
        plotitem.amr_celledges_show = [0] * 10
        plotitem.amr_patchedges_show = [1, 1, 1, 1, 1, 0, 0]
    elif plot_type == 'contour':
        pass


def add_land(plotaxes, plot_type='pcolor', bounds=[-10, 10]):
    """Add plotitem for land"""

    if plot_type == 'pcolor':
        plotitem = plotaxes.new_plotitem(name='land', plot_type='2d_pcolor')
        plotitem.show = True
        plotitem.plot_var = geoplot.land
        plotitem.pcolor_cmap = land_cmap
        plotitem.pcolor_cmin = bounds[0]
        plotitem.pcolor_cmax = bounds[1]
        plotitem.add_colorbar = False
        plotitem.amr_celledges_show = [0] * 10
        plotitem.amr_patchedges_show = [1, 1, 1, 1, 1, 0, 0]

    elif plot_type == 'contour':
        plotitem = plotaxes.new_plotitem(name="land", plot_type='2d_contour')
        plotitem.plot_var = geoplot.land
        plotitem.contour_nlevels = 40
        plotitem.contour_min = bounds[0]
        plotitem.contour_max = bounds[1]
        plotitem.amr_contour_colors = ['g']  # color on each level
        plotitem.amr_patch_bgcolor = ['#ffeeee', '#eeeeff', '#eeffee']
        plotitem.celledges_show = 0
        plotitem.patchedges_show = 0


def add_bathy_contours(plotaxes, contour_levels=None, color='k'):
    """Add plotitem to plot contours of the topography"""

    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.plot_var = geoplot.topo
    if contour_levels is None:
        contour_levels = [0.0]
    plotitem.contour_levels = contour_levels
    plotitem.amr_contour_colors = [color]
    plotitem.kwargs = {'linestyles': 'solid', 'linewidths': 2}
    plotitem.amr_contour_show = [1] * 10
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0


# ===== Storm related plotting =======
def sec2days(seconds):
    """Converst seconds to days."""
    return seconds / (60.0**2 * 24.0)


def plot_track(t, x, y, wind_radius, wind_speed, Pc, name=None):
    r"""Plot hurricane track given a storm.data file"""

    if name is None:
        name = ""
    else:
        name = " - %s" % name

    colors = ['r', 'b']
    divide = (np.max(Pc) + np.min(Pc)) / 2.0

    fig = plt.figure(1)
    axes = fig.add_subplot(111)
    indices = Pc < divide
    axes.scatter(x[indices], y[indices], color='r', marker='o')
    indices = Pc >= divide
    axes.scatter(x[indices], y[indices], color='b', marker='o')
    axes.set_title("Track%s" % name)

    fig = plt.figure(2, figsize=(24, 6))
    axes = fig.add_subplot(131)
    axes.plot(sec2days(t), wind_speed)
    axes.set_title("Maximum Wind Speed%s" % name)

    axes = fig.add_subplot(132)
    axes.plot(sec2days(t), wind_radius)
    axes.set_title("Maximum Wind Radius%s" % names)

    axes = fig.add_subplot(133)
    axes.plot(sec2days(t), Pc)
    axes.plot(sec2days(t), np.ones(t.shape) * divide, 'k--')
    axes.set_title("Central Pressure%s" % name)
