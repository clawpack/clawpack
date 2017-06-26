#!/usr/bin/env python
# encoding: utf-8
r"""
Module defines a class and routines for managing parameterized storm input.

:Formats Supported:

:Models Supported:

"""

from __future__ import print_function
from __future__ import absolute_import

import sys

import numpy
import datetime 

#import clawpack.geoclaw.units as units
                            
#                           days   s/hour   hour/day                             
days2seconds = lambda days: days * 60.0**2 * 24.0

import clawpack.geoclaw.units as units

# Define supported formats and models
_supported_formats = ["GEOCLAW", "HURDAT", "HURDAT2", "JMA", "IMD"]
_supported_models = ["holland_1980", "holland_2010", "cle_2015"]

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
     1. Read in existing file at *path*.
     2. Construct an empty storm and supply the fields needed.  Note that these
        fields must be converted to the appropriate units.

    :Input:
     - *path* (string) Path to file to be read in if requested.
     - *file_format* (string) Format of file at path.  Default is "hurdata2"
     - *kwargs* (dict) Other key-word arguments are passed to the appropriate
       read routine.
    """

    def __init__(self, path=None, file_format="hurdat2", **kwargs):
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
            self.read(path, file_format=file_format, **kwargs)

    # =========================================================================
    # Read Routines
    def read(self, path, file_format="hurdat2", **kwargs):
        r""""""
        if file_format.upper() not in _supported_formats:
            raise ValueError("File type not one of supproted formats.")
        else:

            data_file = open(path,'r')
 
            date = []
            eye_loc = [] 
            mws = []
            Pc = []
            rrp = []
            mwr = []

            for line in data_file: 
                line = line.split(',')
                date.append(self.date2seconds(line[2].strip()))
                if line[6].strip()[-1] == 'N':
                    lat = (float(line[6].strip()[0:-1])/10) 
                else: 
                    lat = (-1*float(line[6].strip()[0:-1])/10)
                if line[7].strip()[-1] == 'E': 
                    lon = (float(line[7].strip()[0:-1])/10) 
                else: 
                    lon = (-1*float(line[7].strip()[0:-1])/10)
                eye_loc.append([lat,lon]) 
                mws.append(int(line[8].strip()))
                Pc.append(int(line[9].strip()))
                rrp.append(int(line[18].strip()))
                mwr.append(int(line[19].strip()))

            self.t = date 
            self.eye_location = eye_loc
            self.max_wind_speed = mws 
            self.max_wind_radius = mwr 
            self.central_pressure = Pc 
            self.storm_radius = rrp 
            data_file.close() 
 
    def write(self, path, file_format="geoclaw"):
        r""""""
        with open('%s/hurdat.%s' %(path,file_format), 'w') as storm_file: 
            for i in range(len(self.t)): 
                storm_file.write("%s" %(str(self.t[i])))
                storm_file.write("%s" %(" ".rjust(1)))
                storm_file.write("%s %s" %(str(self.eye_location[i][0]).rjust(4), str(self.eye_location[i][1]).rjust(4)))
                storm_file.write("%s" % (" ".rjust(1)))
                storm_file.write("%s" % (str(self.max_wind_speed[i]).rjust(3)))
                storm_file.write("%s" % (" ".rjust(1)))
                storm_file.write("%s" % (str(self.central_pressure[i]).rjust(4)))
                storm_file.write("%s" % (" ".rjust(1)))
                storm_file.write("%s" % (str(self.storm_radius[i]).rjust(3)))
                storm_file.write("%s" % (" ".rjust(1)))
                storm_file.write("%s" % (str(self.max_wind_radius[i]).rjust(3)))
                storm_file.write("\n")
    
    def date2seconds(self,date): 
        r"""
        Parse date to convert to seconds. 
        """
        year = int(date[0:4])
        month = int(date[4:6])
        day = int(date[6:8]) 
        hour = int(date[8:])
        new_year_day = datetime.datetime(year,1,1,0)  
        d = datetime.datetime(year, month,day,hour) - new_year_day 
        d = days2seconds(d.days) + d.seconds
        return d  
            
    def wind(self, x, t):
            raise ValueError("File format %s not available." % file_format)

        getattr(self, 'read_%s' % file_format.lower())(path, **kwargs)

    def read_geoclaw(self, path):
        r"""Read in a GeoClaw formatted storm file

        GeoClaw storm files are read in by the Fortran code and are not meant
        to be human readable.

        :Input:
         - *path* (string) Path to the file to be read.
        """

        with open(path, 'r') as data_file:
            num_casts = int(data_file.readline())
            data = numpy.loadtxt(path)

        num_forecasts = data.shape[0]
        self.t = data[:, 0]
        self.eye_location = data[:, 1]
        self.max_wind_speed = data[:, 2]
        self.max_wind_radius = data[:, 3]
        self.central_pressure = data[:, 4]
        self.storm_radius = data[:, 5]

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
        r"""Write out a GeoClaw formatted storm file

        GeoClaw storm files are read in by the Fortran code and are not meant
        to be human readable.

        :Input:
         - *path* (string) Path to the file to be written.
        """

        with open(path, 'w') as data_file:
            data_file.write("%s\n" % self.t.shape[0])
            for n in range(self.t.shape[0]):
                data_file.write("%s %s %s %s %s %s %s" %
                                                (self.t[n],
                                                 self.eye_location[n, 0],
                                                 self.eye_location[n, 1],
                                                 self.max_wind_speed[n],
                                                 self.max_wind_radius[n],
                                                 self.central_pressure[n],
                                                 self.storm_radius[n]))

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


# =============================================================================
# Model field construction - Models supported are
#  - Holland 1980 ('HOLLAND_1980') [1]
#  - Holland 2010 ('HOLLAND_2010') [2]
#  - Chavas, Lin, Emmanuel ('CLE_2015') [3]
# *TODO* - Add citations
#
# In the case where the field is not rotationally symmetric then the r value
# defines the x and y axis extents.
def construct_fields(storm, r, t, model="holland_1980"):
    r""""""

    if model.lower() not in _supported_models:
        raise ValueError("Model %s not available." % model)

    return getattr(sys.modules[__name__], model.lower())(storm, x, t)


# Specific implementations
def holland_1980(storm, r, t):
    r""""""
    raise NotImplementedError("Holland 1980 model has not been implemeted.")
    return None, None


def holland_2010(storm, r, t):
    r""""""
    raise NotImplementedError("Holland 2010 model has not been implemeted.")
    return None, None


def cle_2015(storm, r, t):
    r""""""
    raise NotImplementedError("CLE 2015 model has not been implemeted.")
    return None, None


# =============================================================================
# Utility functions
def available_formats():
    r"""Construct a string suitable for listing available storm file formats.
    """
    return ""


def available_models():
    r"""Construct a string suitable for listing available storm models.
    """
    return ""


if __name__ == '__main__':
    Test_Hurdat = Storm('test_files/hurdat.test')   
    Test_Hurdat.read('test_files/hurdat.test') 
    Test_Hurdat.write('test_files') 
    # TODO:  Add commandline ability to convert between formats
    construct_fields(None, None, None)
