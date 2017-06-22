#!/usr/bin/env python
# encoding: utf-8
r"""
Module defines a class and routines for managing parameterized storm input.

:Formats Supported:

:Models Supported:

"""

from __future__ import print_function
from __future__ import absolute_import

import numpy

import clawpack.geoclaw.units as units

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
    _supported_formats = ["GEOCLAW", "HURDAT", "HURDAT2", "JMA", "IMD"]

    def __init__(self, path=None, **kwargs):
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
            self.read(path, **kwargs)

    # =========================================================================
    # Read Routines
    def read(self, path, file_format="hurdat2"):
        r""""""

        if file_format.upper() not in _supported_formats:
            raise ValueError("File format %s not available." % file_format)

        getattr(self, 'read_%s' % file_format.lower())(path)

    def read_geoclaw(self, path):
        r""""""
        raise NotImplementedError("GeoClaw format not fully implemented.")

    def read_hurdat(self, path):
        r""""""
        raise NotImplementedError("HURDAT format not fully implemented.")

    def read_hurdat2(self, path):
        r""""""
        raise NotImplementedError("HURDAT2 format not fully implemented.")

    def read_jma(self, path):
        r""""""
        raise NotImplementedError("JMA format not fully implemented.")

    def read_imd(self, path):
        r""""""
        raise NotImplementedError("IMD format not fully implemented.")

    # =========================================================================
    # Write Routines
    def write(self, path, file_format="geoclaw"):
        r""""""

        if file_format.upper() not in _supported_formats:
            raise ValueError("File format %s not available." % file_format)

        getattr(self, 'write_%s' % file_format.lower())(path)

    def write_geoclaw(self, path):
        r""""""
        raise NotImplementedError("GeoClaw format not fully implemented.")

    def write_hurdat(self, path):
        r""""""
        raise NotImplementedError("HURDAT format not fully implemented.")

    def write_hurdat2(self, path):
        r""""""
        raise NotImplementedError("HURDAT2 format not fully implemented.")

    def write_jma(self, path):
        r""""""
        raise NotImplementedError("JMA format not fully implemented.")

    def write_imd(self, path):
        r""""""
        raise NotImplementedError("IMD format not fully implemented.")


# =========================================================================
# Model field construction - Models supported are
#  - Holland 1980 [1]
#  - Holland 2010 [2]
#  - Chavas, Lin, Emmanuel [3]
#
def wind(self, x, t, model="holland_80"):
    r""""""
    pass


def pressure(self, x, t, model="holland_80"):
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
