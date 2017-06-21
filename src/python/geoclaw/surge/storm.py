#!/usr/bin/env python
# encoding: utf-8
r"""
Module defines a class and routines for managing parameterized storm input.

:Formats Supported:

"""

from __future__ import print_function
from __future__ import absolute_import

import numpy

# =============================================================================
#  Basic storm class


class Storm(object):
    r"""
    Storm data object

    This object contains a time series of time data that describe a particular
    storm.  This includes the attributes below and the ability to read from
    multiple sources for data such as the U.S. National Hurricane Center (NHC),
    the Japanese Meterological Agency (JMA), and the Indian Meteorlogical
    Department (IMD).  This class can then write out in any of these formats,
    construct the wind and pressure fields using a supported parameterized
    model, or output the GeoClaw supported storm format used for running storm
    surge simulations.

    *TODO:*  Add description of unit handling

    :Attributes:
     - *t* (ndarray(:)) Contains the time at which each entry of the other
       arrays are at.  Default units are seconds.
     - *eye_location* (ndarray(:, :)) location of the eye of the storm.
       Default units are in signed decimcal longitude and latitude.
     - *max_wind_speed (ndarray(:)) Maximum wind speed.  Default units are
       meters/second.
     - *max_wind_radius (ndarray(:)) Radius at which the maximum wind speed
       occurs.  Default units are meters.
     - *central_pressure* (ndarray(:)) Central pressure of storm.  Default
       units are Pascals.
     - *storm_radius* (ndarray(:)) Radius of storm, often defined as the last
       closed iso-bar of pressure.  Default units are meters.

    :Initialization:

    """

    # Define supported formats and models
    _supported_formats = ["HURDAT", "HURDAT2", "JMA", "IMD"]
    _supported_models = ["Holland_80", "Holland_10", "CLE"]

    def __init__(self, path=None):
        r"""Storm Initiatlization Routine

        See :class:`Storm` for more info.
        """

        self.t = None
        self.eye_location = None
        self.max_wind_speed = None
        self.max_wind_radius = None
        self.central_pressure = None
        self.storm_radius = None

        if path is None:
            self.read(path)

    def read(self, path, file_format="hurdat"):
        r""""""

        if file_format.upper() not in _supported_formats:
            raise ValueError("")

    def write(self, path, file_format="geoclaw"):
        r""""""

        pass

    def wind(self, x, t):
        r""""""
        pass

    def pressure(self, t):
        r""""""
        pass

# =============================================================================
# Plotting functions
#
def holland_08_fields():
    r""""""
    pass

# =============================================================================
# Plotting functions
#
def plot_wind(storm, x, t, model='', axes=None):
    r""""""

    if axes is None:
        fig = plt.figure()
        axes = fig.add_subplot(1, 1, 1)

    axes.plot(x, storm.wind(x, t))

    return axes


def plot_pressure(storm, x, t, model='', axes=None):
    r""""""

    if axes is None:
        fig = plt.figure()
        axes = fig.add_subplot(1, 1, 1)

    axes.plot(x, storm.pressure(x, t))

    return axes


if __name__ == '__main__':
    pass
