#!/usr/bin/env python

r"""
Module defines a class and routines for managing parameterized storm input.

:Formats Supported:
    - GeoClaw
    - ATCF
    - HURDAT
    - JMA
    - IMD

:Models Supported:
    - Holland 1980
    - Holland 2010
    - Chavas, Lin, Emmanuel 2015
"""

from __future__ import print_function
from __future__ import absolute_import

from six.moves import range

import sys
import os
import argparse
import datetime

import numpy

import clawpack.geoclaw.units as units
import clawpack.clawutil.data

# =============================================================================
#  Common acronyms across formats

# ATCF basins with their expanded names
# see https://www.nrlmry.navy.mil/atcf_web/docs/database/new/abrdeck.html
ATCF_basins = {"AL": "Atlantic",
               "CP": "Central Pacific",
               "EP": "East Pacific",
               "IO": "North Indian Ocean",
               "SH": "Southern Hemisphere",
               "SL": "Southern Atlantic",
               "LS": "Southern Atlantic",
               "WP": "North West Pacific"}

# Tropical Cyclone Designations
# see https://www.nrlmry.navy.mil/atcf_web/docs/database/new/abrdeck.html
TC_designations = {"DB": "disturbance",
                   "TD": "tropical depression",
                   "TS": "tropical storm",
                   "TY": "typhoon",
                   "ST": "super typhoon",
                   "TC": "tropical cyclone",
                   "HU": "hurricane",
                   "SD": "subtropical depression",
                   "SS": "subtropical storm",
                   "EX": "extratropical systems",
                   "IN": "inland",
                   "DS": "dissipating",
                   "LO": "low",
                   "WV": "tropical wave",
                   "ET": "extrapolated",
                   "XX": "unknown"}

# HURDAT special designations
# see http://www.aoml.noaa.gov/hrd/data_sub/newHURDAT.html
hurdat_special_entries = {"L": "landfall",
                          "W": "max wind",
                          "P": "min pressure",
                          "I": "max intensity",
                          "C": "closest approach",
                          "S": "status change",
                          "G": "genesis",
                          "T": "additional track point"}


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
     - *t* (list(datetime.datetiem)) Contains the time at which each entry of
       the other arrays are at.  These are expected to be *datetime* objects.
       Note that when written some formats require a *time_offset* to be set.
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
     - *time_offset* (datetime.datetime) A date time that as an offset for the
       simulation time.  This will default to the beginning of the first of the
       year that the first time point is found in.

    :Initialization:
     1. Read in existing file at *path*.
     2. Construct an empty storm and supply the fields needed.  Note that these
        fields must be converted to the appropriate units.

    :Input:
     - *path* (string) Path to file to be read in if requested.
     - *file_format* (string) Format of file at path.  Default is "hurdat"
     - *kwargs* (dict) Other key-word arguments are passed to the appropriate
       read routine.
    """

    # Define supported formats and models - keys are function name related and
    # values are the proper name and a citation or URL documenting the format
    _supported_formats = {"geoclaw": ["GeoClaw", "http://www.clawpack.org/storms"],
                          "atcf": ["ATCF", "http://www.nrlmry.navy.mil/atcf_web/docs/database/new/database.html"],
                          "hurdat": ["HURDAT", "http://www.aoml.noaa.gov/hrd/hurdat/Data_Storm.html"],
                          "jma": ["JMA", "http://www.jma.go.jp/jma/jma-eng/jma-center/rsmc-hp-pub-eg/Besttracks/e_format_bst.html"],
                          "imd": ["IMD", "http://www.rsmcnewdelhi.imd.gov.in/index.php"],
                          "tcvitals": ["TC-Vitals", None]}

    def __init__(self, path=None, file_format="ATCF", **kwargs):
        r"""Storm Initiatlization Routine

        See :class:`Storm` for more info.
        """

        self.t = None
        self.time_offset = None
        self.eye_location = None
        self.max_wind_speed = None
        self.max_wind_radius = None
        self.central_pressure = None
        self.storm_radius = None

        # Storm descriptions - not all formats provide these
        self.name = None
        self.basin = None                   # Basin containing storm
        self.ID = None                      # ID code - depends on format
        self.classification = None          # Classification of storm (e.g. HU)
        self.event = None                   # Event (e.g. landfall) - HURDAT

        if path is not None:
            self.read(path, file_format=file_format, **kwargs)

    # ==========================================================================
    #  Basic object support
    def __str__(self):
        r""""""
        output = "Name: %s" % self.name
        output = "\n".join((output, "Dates: %s - %s" % (self.t[0].isoformat(),
                                                        self.t[-1].isoformat())
                            ))
        return output

    def __repr__(self):
        return "None"

    # ==========================================================================
    # Read Routines
    def read(self, path=None, file_format="atcf", **kwargs):
        r"""Read in storm data from *path* with format *file_format*

        :Input:
         - *path* (string) Path to data file.
         - *file_format (string) Format of the data file.  See list of
           supported formats for a list of valid strings.  Defaults to
           "hurdat".
         - *kwargs* (dict) Keyword dictionary for additional arguments that can
           be passed down to the appropriate read functions.  Please refer to
           the specific routine for a list of valid options.

        :Raises:
         - *ValueError* If the *file_format* requested does not match any of
           the available supported formats a *ValueError* is raised.
        """

        # If a path is not provided then we can try and find the relevant
        # database and download it
        if path is None:
            data_str = ("Currently automatic download of storm databases is ", 
                        "not implemented.  Please refer to the URLs below for", 
                        "references as to where you can download storm data",
                        "files:",
                        " - ATCF - http://ftp.nhc.noaa.gov/atcf/archive/",
                        " - HURDAT - http://www.aoml.noaa.gov/hrd/hurdat/Data_Storm.html",
                        " - JMA - http://www.jma.go.jp/jma/jma-eng/jma-center/rsmc-hp-pub-eg/besttrack.html",
                        " - IMD - http://www.rsmcnewdelhi.imd.gov.in/index.php")
            raise NotImplementedError("\n".join(data_str))

        if file_format.lower() not in self._supported_formats.keys():
            raise ValueError("File format %s not available." % file_format)

        getattr(self, 'read_%s' % file_format.lower())(path, **kwargs)

    def read_geoclaw(self, path, verbose=False):
        r"""Read in a GeoClaw formatted storm file

        GeoClaw storm files are read in by the Fortran code and are not meant
        to be human readable.

        :Input:
         - *path* (string) Path to the file to be read.
        """

        with open(path, 'r') as data_file:
            num_casts = int(data_file.readline())
            self.time_offset = datetime.datetime.strptime(
                                                      data_file.readline()[:19],
                                                      '%Y-%m-%dT%H:%M:%S')

        data = numpy.loadtxt(path, skiprows=3)
        num_forecasts = data.shape[0]
        self.eye_location = numpy.empty((2, num_forecasts))
        assert(num_casts == num_forecasts)
        self.t = [self.time_offset + datetime.timedelta(seconds=data[i, 0])
                  for i in range(num_forecasts)]
        self.eye_location[0, :] = data[:, 1]
        self.eye_location[1, :] = data[:, 2]
        self.max_wind_speed = data[:, 3]
        self.max_wind_radius = data[:, 4]
        self.central_pressure = data[:, 5]
        self.storm_radius = data[:, 6]

    def read_atcf(self, path, verbose=False):
        r"""Read in a ATCF formatted storm file

        ATCF format has storm stored individually so there is no support for
        multiple storms in a particular file.

        :Input:
         - *path* (string) Path to the file to be read.
        """

        # No header, can assume storm data
        data_block = []
        with open(path, 'r') as ATCF_file:
            for line in ATCF_file:
                line = line.split(",")
                line = [value.strip() for value in line]
                data_block.append(line)
        num_lines = len(data_block)

        # Parse data block - convert to correct units
        # Conversions:
        #  max_wind_speed - Convert knots to m/s - 0.51444444
        #  max_wind_radius  - convert from nm to m - 1.8520000031807990 * 1000.0
        #  central_pressure - convert from mbar to Pa - 100.0
        #  Radius of last isobar contour - convert from nm to m - 1.852000003180799d0 * 1000.0
        self.t = []
        self.event = numpy.empty(num_lines, dtype=str)
        self.classification = numpy.empty(num_lines, dtype=str)
        self.eye_location = numpy.empty((num_lines, 2))
        self.max_wind_speed = numpy.empty(num_lines)
        self.central_pressure = numpy.empty(num_lines)
        self.max_wind_radius = numpy.empty(num_lines)
        self.storm_radius = numpy.empty(num_lines)

        for (i, line) in enumerate(data_block):
            # Skip empty linies
            if len(line) == 0:
                break
            data = line

            # Mark if this is a shortened line
            short_data = len(data) < 19

            # Create time
            self.t.append(datetime.datetime(int(data[2][:4]),
                                            int(data[2][4:6]),
                                            int(data[2][6:8]),
                                            int(data[2][-2:])))

            # If an event is occuring record it.  If landfall then use as an
            # offset.   Note that if there are multiple landfalls the last one
            # is used as the offset
            if not short_data:
                if len(data[22].strip()) > 0:
                    self.event[i] = data[22].strip()
                    if self.event[i].upper() == "L":
                        self.time_offset = self.t[i]

            # Classification, note that this is not the category of the storm
            self.classification[i] = data[10]

            # Parse eye location - longitude/latitude order
            if data[6][-1] == "N":
                self.eye_location[i, 1] = float(data[6][0:-1]) / 10.0
            else:
                self.eye_location[i, 1] = -float(data[6][0:-1]) / 10.0
            if data[7][-1] == "E":
                self.eye_location[i, 0] = float(data[7][0:-1]) / 10.0
            else:
                self.eye_location[i, 0] = -float(data[7][0:-1]) / 10.0

            # Intensity information
            self.max_wind_speed[i] = float(data[8]) * 0.51444444
            self.central_pressure[i] = float(data[9]) * 100.0
            if not short_data:
                self.storm_radius[i] = (float(data[18]) * 1.852000003180799
                                        * 1000.0)
                self.max_wind_radius[i] = (float(data[19]) * 1.852000003180799
                                           * 1000.0)

    def read_hurdat(self, path, single_storm=True, name=None, year=None,
                                verbose=False):
        r"""Read in HURDAT formatted storm file

        This is the current version of HURDAT data available (HURDAT 2).  Note
        that if the file contains multiple storms only the first will be read
        in unless name and/or year is provided.

        For more details on the HURDAT format and getting data see

        http://www.aoml.noaa.gov/hrd/hurdat/Data_Storm.html

        :Input:
         - *path* (string) Path to the file to be read.
         - *single_storm* (bool) If *True* then this file contains one storm.
           Default is *single_storm = False*.
         - *name* (string) If the file contains multiple storms use *name* to
           search for the correct storm.  If there are multiple storms with the
           same name then the first one encountered is read in.
         - *year* (int) Additional filtering criteria.  If there are multiple
           storms with the same name use the year of the storm to pick out the
           right one.

        :Raises:
         - *ValueError* If the method cannot find the name/year matching the
           storm or they are not provided when *single_storm == False* then a
           value error is risen.
        """

        # If no name and/or year are provided then we read until the end of the
        # file or until we reach another header line
        if not single_storm:
            warning_str = "\n".join("Reading from the database is currently is",
                                    "currently not supported but is planned",
                                    "for a future release.")
            raise NotImplementedError(warning_str)

            if name is None:
                err_msg = "".join(("Input indicated that there was more than ",
                                   "one storm in the file being read.  If ",
                                   "this is the case then a name must be ",
                                   "provided to pick out the storm to be read",
                                   "in."))
                raise ValueError(err_msg)

            with open(path, 'r') as hurdat_file:
                success = False
                for (n, line) in enumerate(hurdat_file):
                    if line[:2] in ATCF_basins.keys():
                        # This is a header line
                        storm_year = int(line.split(",")[0].strip()[4:])
                        storm_name = line.split(',')[1].strip()
                        num_lines = int(line.split(",")[2].strip())

                        if name.lower() == storm_name.lower():
                            if year is not None:
                                if year == storm_year:
                                    # Take this storm
                                    success = True
                                    break
                            else:
                                # Take this storm
                                success = True
                                break

                # Extract data chunk
                if success:
                    self.name = storm_name
                    self.ID = line.split(",")[0].strip()
                    self.basin = self.ID[:2]

                    data_block = ""
                    for n in range(num_lines):
                        line = hurdat_file.readline()
                        data_block = "".join((data_block, line))
                    data_block = data_block[:-1].split('\n')
                    assert len(data_block) == num_lines

                else:
                    # Return error based on failure of criteria
                    err_msg = "".join(("Name %s or year %s " % (name, year),
                                       "did not match available ",
                                       "values.  Please check to make sure",
                                       "the name and year are present in the",
                                       "file you provided."))
                    raise ValueError(err_msg)

        else:
            # No header, just assume storm data
            data_block = []
            with open(path, 'r') as hurdat_file:
                data_block.append(hurdat_file.readlines())

        # Parse data block
        self.t = []
        self.event = numpy.empty(num_lines, dtype=str)
        self.classification = numpy.empty(num_lines, dtype=str)
        self.eye_location = numpy.empty((num_lines, 2))
        self.max_wind_speed = numpy.empty(num_lines)
        self.central_pressure = numpy.empty(num_lines)
        self.max_wind_radius = numpy.empty(num_lines)
        self.storm_radius = numpy.empty(num_lines)
        for (i, line) in enumerate(data_block):
            if len(line) == 0:
                break
            data = [value.strip() for value in line.split(",")]

            # Create time
            self.t.append(datetime.datetime(int(data[0][:4]),
                                            int(data[0][4:6]),
                                            int(data[0][6:8]),
                                            int(data[1][:2]),
                                            int(data[1][2:])))

            # If an event is occuring record it.  If landfall then use as an
            # offset.   Note that if there are multiple landfalls the last one
            # is used as the offset
            if len(data[2].strip()) > 0:
                self.event[i] = data[2].strip()
                if self.event[i].upper() == "L":
                    self.time_offset = self.t[i]

            # Classification, note that this is not the category of the storm
            self.classification[i] = data[3]

            # Parse eye location
            if data[4][-1] == "N":
                self.eye_location[i, 0] = float(data[4][0:-1])
            else:
                self.eye_location[i, 0] = -float(data[4][0:-1])
            if data[5][-1] == "E":
                self.eye_location[i, 1] = float(data[5][0:-1])
            else:
                self.eye_location[i, 1] = -float(data[5][0:-1])

            # Intensity information
            self.max_wind_speed[i] = float(data[6])
            self.central_pressure[i] = float(data[7])
            self.max_wind_radius[i] = float(data[8])
            self.storm_radius[i] = float(data[9])

    def read_jma(self, path, single_storm=True, name=None, year=None,
                             verbose=False):
        r"""Read in JMA formatted storm file

        Note that if the file contains multiple storms only the first will be
        read in unless name and/or year is provided.

        For more details on the JMA format and getting data see

        http://www.jma.go.jp/jma/jma-eng/jma-center/rsmc-hp-pub-eg/Besttracks/e_format_bst.html

        :Input:
         - *path* (string) Path to the file to be read.
         - *single_storm* (bool) If *True* then this file contains one storm.
           Default is *single_storm = False*.
         - *name* (string) If the file contains multiple storms use *name* to
           search for the correct storm.  If there are multiple storms with the
           same name then the first one encountered is read in.
         - *year* (int) Additional filtering criteria.  If there are multiple
           storms with the same name use the year of the storm to pick out the
           right one.
        """

        if not single_storm:
            warning_str = "\n".join("Reading from the database is currently is",
                                    "currently not supported but is planned",
                                    "for a future release.")
            raise NotImplementedError(warning_str)
            if name is None:
                err_msg = "".join(("Input indicated that there was more than ",
                                   "one storm in the file being read.  If ",
                                   "this is the case then a name must be ",
                                   "provided to pick out the storm to be read",
                                   " in."))
                raise ValueError(err_msg)

            with open(path, 'r') as JMA_file:
                success = False
                for (n, line) in enumerate(JMA_file):
                    if line[:5] == "66666":
                        # This is a header line
                        storm_year = int(line[6:8])
                        if storm_year > 50:
                            storm_year += 1900
                        else:
                            storm_year += 2000
                        ID = int(line[8:10])
                        num_lines = int(line[13:15])
                        flag = bool(int(line[26]))
                        storm_name = line[30:51].strip()

                        if name.lower() == storm_name.lower():
                            if year is not None:
                                if year == storm_year:
                                    # Take this storm
                                    success = True
                                    break
                            else:
                                # Take this storm
                                success = True
                                break

                # Extract data chunk
                if success:
                    self.name = storm_name
                    self.ID = ID

                    data_block = ""
                    for n in range(num_lines):
                        line = JMA_file.readline()
                        data_block = "".join((data_block, line))
                    data_block = data_block[:-1].split('\n')
                    assert len(data_block) == num_lines

                else:
                    # Return error based on failure of criteria
                    err_msg = "".join(("Name %s or year %s " % (name, year),
                                       "did not match available ",
                                       "values.  Please check to make sure",
                                       "the name and year are present in the",
                                       "file you provided."))
                    raise ValueError(err_msg)

        else:
            # No header, just assume storm data
            data_block = []
            with open(path, 'r') as JMA_file:
                data_block.append(JMA_file.readlines())

        # Parse data block
        self.t = []
        self.event = numpy.empty(num_lines, dtype=str)
        self.classification = numpy.empty(num_lines, dtype=str)
        self.eye_location = numpy.empty((num_lines, 2))
        self.max_wind_speed = numpy.empty(num_lines)
        self.central_pressure = numpy.empty(num_lines)
        self.max_wind_radius = numpy.empty(num_lines)
        self.storm_radius = numpy.empty(num_lines)
        for (i, line) in enumerate(data_block):
            if len(line) == 0:
                break
            data = [value.strip() for value in line.split(" ")]

            # Create time
            self.t.append(datetime.datetime(int(data[0][:4]),
                                            int(data[0][4:6]),
                                            int(data[0][6:8]),
                                            int(data[1][:2]),
                                            int(data[1][2:])))

            # If an event is occuring record it.  If landfall then use as an
            # offset.   Note that if there are multiple landfalls the last one
            # is used as the offset
            if len(data[2].strip()) > 0:
                self.event[i] = data[2].strip()
                if self.event[i].upper() == "L":
                    self.time_offset = self.t[i]

            # Classification, note that this is not the category of the storm
            self.classification[i] = data[3]

            # Parse eye location
            if data[4][-1] == "N":
                self.eye_location[i, 0] = float(data[4][0:-1])
            else:
                self.eye_location[i, 0] = -float(data[4][0:-1])
            if data[5][-1] == "E":
                self.eye_location[i, 1] = float(data[5][0:-1])
            else:
                self.eye_location[i, 1] = -float(data[5][0:-1])

            # Intensity information
            self.max_wind_speed[i] = float(data[6])
            self.central_pressure[i] = float(data[7])
            self.max_wind_radius[i] = float(data[8])
            self.storm_radius[i] = float(data[9])

    def read_imd(self, path, verbose=False):
        r"""Extract relevant hurricane data from IMD file
            and update storm fields with proper values.

        :Input:
         - *path* (string) Path to the file to be read.

        Return ValueError if format incorrect or if file not IMD.
        """
        raise ValueError("File type not implemented yet.")

    def read_tcvitals(self, path, single_storm=True, name=None, year=None,
                                  verbose=False):
        r"""Extract relevant hurricane data from TCVITALS file
            and update storm fields with proper values.

        :Input:
         - *path* (string) Path to the file to be read.

        """

        try:
            import requests
        except ImportError as e:
            print("The 'requests' module is required to read TCVitals data.")
            raise e

        try:
            from bs4 import BeautifulSoup
        except ImportError as e:
            print("BeautfulSoup is required for reading in TCVitals data.")
            raise e

        if not single_storm:
            return "Multiple storms not implemented for this format yet."
        else:
            if int(year) < 2011:
                err_msg = "Years not contained on this page"
                raise ValueError(err_msg)
            else:
                if int(year) == 2016:
                    file_directory_url = "".join((path))
                else:
                    file_directory_url = "".join((path, year, '/'))
                print('File Directory URL:', file_directory_url)
                storm_directory_page = requests.get(file_directory_url)
                soup = BeautifulSoup(storm_directory_page.content,
                                     'html.parser')
                storm_directory_links = soup.find_all('a')
                storm_files = []
                for link in storm_directory_links:
                    if ".dat" in str(link):
                        if "combined" in str(link):
                            continue
                        else:
                            storm_files.append(str(link)[9:35])

                found = False
                for data_file_name in storm_files:
                    data_file_url = "".join((file_directory_url,
                                             data_file_name))
                    unpack = True
                    data_path = clawpack.clawutil.data.get_remote_file(
                                                  data_file_url, unpack=unpack)
                    with open(data_path, 'r') as data:
                        for line in data:
                            self.name = line.split()[2]
                            if name.upper() == self.name.upper():
                                found = True
                                break
                            else:
                                continue
                    data.close()

                    if found:
                        storm_file_url = data_file_url
                        self.name = name.upper()
                        print('Storm File URL', storm_file_url)
                        os.remove(data_path)
                        break
                    else:
                        os.remove(data_path)
                        continue

                if not found:
                    return("Storm not found for the year you specified.")

                storm_path = clawpack.clawutil.data.get_remote_file(
                                                   storm_file_url, unpack=True)

                data_block = []
                with open(storm_path, 'r') as tcvitals_file:
                    data_block = tcvitals_file.readlines()
                print('data_block', data_block)

                num_lines = len(data_block)

                # Parse data block
                self.t = []
                self.event = numpy.empty(num_lines, dtype=str)
                self.classification = numpy.empty(num_lines, dtype=str)
                self.eye_location = numpy.empty((num_lines, 2))
                self.max_wind_speed = numpy.empty(num_lines)
                self.central_pressure = numpy.empty(num_lines)
                self.max_wind_radius = numpy.empty(num_lines)
                self.storm_radius = numpy.empty(num_lines)

                len_data = []

                for (i, line) in enumerate(data_block):
                    if len(line) == 0:
                        break
                    data = [value.strip() for value in line.split()]

                    self.t.append(datetime.datetime(int(data[3][0:4]),
                                                    int(data[3][4:6]),
                                                    int(data[3][6:8]),
                                                    int(data[4][0:2]),
                                                    int(data[4][2:])))

                    self.event[i] = data[1][-1]
                    if self.event[i] == 'L':
                        self.time_offset = self.t[i]

                    if data[5][-1] == 'N':
                        self.eye_location[i, 0] = float(data[5][0:-1])/10
                    else:
                        self.eye_location[i, 0] = -float(data[5][0:-1])/10
                    if data[6][-1] == "E":
                        self.eye_location[i, 1] = float(data[6][0:-1])/10
                    else:
                        self.eye_location[i, 1] = -float(data[6][0:-1])/10

                    # Intensity Information
                    self.max_wind_speed[i] = float(data[8])
                    self.central_pressure[i] = float(data[9])
                    self.max_wind_radius[i] = float(data[11])
                    self.storm_radius[i] = float(data[13])

    # =========================================================================
    # Write Routines
    def write(self, path, file_format="geoclaw", **kwargs):
        r"""Write out the storm data to *path* in format *file_format*

        :Input:
         - *path* (string) Path to data file.
         - *file_format (string) Format of the data file.  See list of
           supported formats for a list of valid strings.  Defaults to
           "geoclaw".
         - *kwargs* (dict) Keyword dictionary for additional arguments that can
           be passed down to the appropriate write functions.  Please refer to
           the specific routine for a list of valid options.

        :Raises:
         - *ValueError* If the *file_format* requested does not match any of
           the available supported formats a *ValueError* is raised.
        """

        if file_format.lower() not in self._supported_formats.keys():
            raise ValueError("File format %s not available." % file_format)

        getattr(self, 'write_%s' % file_format.lower())(path)

    def write_geoclaw(self, path, verbose=False):
        r"""Write out a GeoClaw formatted storm file

        GeoClaw storm files are read in by the GeoClaw Fortran code.

        :Input:
         - *path* (string) Path to the file to be written.
        """

        # Remove duplicate times
        num_casts = len(self.t)
        for n in range(1, len(self.t)):
            if self.t[n] == self.t[n - 1]:
                num_casts -= 1

        try:
            with open(path, 'w') as data_file:
                data_file.write("%s\n" % num_casts)
                data_file.write("%s\n\n" % self.time_offset.isoformat())
                for n in range(len(self.t)):
                    if n > 0:
                        if self.t[n] == self.t[n - 1]:
                            continue

                    format_string = ("{:19,.8e} " * 7)[:-1] + "\n"
                    data = ((self.t[n] - self.time_offset).total_seconds(),
                            self.eye_location[n, 0],
                            self.eye_location[n, 1],
                            self.max_wind_speed[n],
                            self.max_wind_radius[n],
                            self.central_pressure[n],
                            self.storm_radius[n])
                    data_file.write(format_string.format(*data))

        except Exception as e:
            # Remove possibly partially generated file if not successful
            if os.path.exists(path):
                os.remove(path)
            raise e

    def write_atcf(self, path, verbose=False):
        r"""Write out a ATCF formatted storm file

        :Input:
         - *path* (string) Path to the file to be written
        """
        try:
            with open(path, 'w') as data_file:
                for n in range(self.t.shape[0]):
                    data_file.write("".join((", " * 2,
                                         "%s" % seconds2date(self.t[n]),
                                         ", " * 4,
                                         "%s" % (int(self.eye_location[n, 0] *
                                                     10.0)),
                                         ", ",
                                         "%s" % (int(self.eye_location[n, 1] *
                                                     10.0)),
                                         ", ",
                                         "%s" % self.max_wind_speed[n],
                                         ", ",
                                         "%s" % self.central_pressure[n],
                                         ", ",
                                         ", " * 8,
                                         "%s" % self.storm_radius[n],
                                         ", ",
                                         "%s" % self.max_wind_radius[n],
                                         ", " * 10,
                                         "\n")))
        except Exception as e:
            # Remove possiblly partially generated file if not successful
            if os.path.exists(path):
                os.remove(path)
            raise e

    def write_hurdat(self, path, verbose=False):
        r"""Write out a HURDAT formatted storm file

        :Input:
         - *path* (string) Path to the file to be written
        """
        try:
            with open(path, 'w') as data_file:
                data_file.write('%s %s %s' % ("Date", "Hurricane Name",
                                              "Indicator"))
                for n in range(self.t.shape[0]):

                    latitude = float(self.eye_location[n, 0])
                    longitude = float(self.eye_location[n, 1])

                    # Convert latitude to proper Hurdat format e.g 12.0N
                    if latitude > 0:
                        latitude = str(numpy.abs(latitude)) + 'N'
                    else:
                        latitude = str(numpy.abs(latitude)) + 'S'

                    # Convert longitude to proper Hurdat format e.g 12.0W
                    if longitude > 0:
                        longitude = str(numpy.abs(longitude)) + 'E'
                    else:
                        longitude = str(numpy.abs(longitude)) + 'W'

                    data_file.write("".join(("%s" % self.seconds2date(
                                                              self.t[n])[0:-2],
                                         "%s00" % self.seconds2date(
                                                              self.t[n])[-2:],
                                         ", " * 3,
                                         "%s" % (latitude),
                                         ", ",
                                         "%s" % (longitude),
                                         ", ",
                                         "%s" % self.max_wind_speed[n],
                                         ", ",
                                         "%s" % self.central_pressure[n],
                                         ", ",
                                         "%s" % self.storm_radius[n],
                                         ", ",
                                         "%s" % self.max_wind_radius[n],
                                         ", " * 10,
                                         "\n")))
        except Exception as e:
            # Remove possiblly partially generated file if not successful
            if os.path.exists(path):
                os.remove(path)
            raise e

    def write_jma(self, path, verbose=False):
        r"""Write out a JMA formatted storm file

        :Input:
         - *path* (string) Path to the file to be written
        """
        try:
            with open(path, 'w') as data_file:
                for n in range(self.t.shape[0]):
                    data_file.write("".join(("%s" % self.seconds2date(self.t[n]),
                                         " " * 4,
                                         "%s" % (int(self.eye_location[n, 0] *
                                                     10.0)),
                                         ", ",
                                         "%s" % (int(self.eye_location[n, 1] *
                                                     10.0)),
                                         ", ",
                                         "%s" % self.max_wind_speed[n],
                                         ", ",
                                         "%s" % self.central_pressure[n],
                                         ", ",
                                         ", " * 8,
                                         "%s" % self.storm_radius[n],
                                         ", ",
                                         "%s" % self.max_wind_radius[n],
                                         ", " * 10,
                                         "\n")))
        except Exception as e:
            # Remove possiblly partially generated file if not successful
            if os.path.exists(path):
                os.remove(path)
            raise e

    def write_imd(self, path, verbose=False):
        r"""Write out a IMD formatted storm file

        :Input:
         - *path* (string) Path to the file to be written
        """
        raise NotImplementedError("IMD format not fully implemented.")

    # =========================================================================
    # Other Useful Routines
    def plot(self, axes=None, intensity=False, limits=None, track_color='red',
                   category_colors=None, categorization="NHC"):
        r"""Plot the track and optionally the strength of the storm

        """

        # TODO:  Switch to cartopy plotting
        import matplotlib.pyplot as plt
        from mpl_toolkits.basemap import Basemap

        if axes is None:
            fig = plt.figure()
            axes = fig.add_subplot(1, 1, 1)

        # limits = ((long), (lat))
        if limits is None:
            raise NotImplementedError("Need to do this...")

        if category_color is None:
            category_color = {5: 'red',
                              4: 'yellow',
                              3: 'orange',
                              2: 'green',
                              1: 'blue',
                              0: 'gray'}

        mapping = Basemap()
        longitude, latitude = mapping(self.eye_location[:, 0],
                                      self.eye_location[:, 1])
        category = self.category(categorization=categorization)
        for i in range(len(longitude)):
            if intensity:
                color = category_color[category[i]]
            else:
                color = track_color
            mapping.plot(longitude[i:i + 2], latitude[i:i + 2], color=color)

        mapping.drawcoastlines()
        mapping.drawcountries()
        mapping.fillcontinents()

        return axes

    def category(self, categorization="NHC", cat_names=False):
        r"""Categorizes storm based on relevant storm data

        :Input:
         - *categorization* (string) Type of categorization to use.  Defaults
           to the National Hurricane Center "NHC".
         - *cat_names* (bool) If True returns the category name rather than a
           number.  Default to *False*.

        :Output:
         - (ndarray) Integer array of categories at each time point of the
           storm.
         - (list) Similar to the above but the name of the category as a
           *string*.  This is only returned if *car_names = True*.

        """

        # TODO:  Need to standardize on 1-minute (almost never available) or
        # 10-minute (widely available) - see
        # https://en.wikipedia.org/wiki/Tropical_cyclone#Major_basins_and_related_warning_centers

        if categorization.upper() == "BEAUFORT":
            # Beaufort scale below uses knots
            speeds = units.convert(self.max_wind_speed, "m/s", "knots")
            category = (numpy.zeros(speeds.shape) +
                        (speeds >= 1) * (speeds < 4) * 1 +
                        (speeds >= 4) * (speeds < 7) * 2 +
                        (speeds >= 7) * (speeds < 11) * 3 +
                        (speeds >= 11) * (speeds < 17) * 4 +
                        (speeds >= 17) * (speeds < 22) * 5 +
                        (speeds >= 22) * (speeds < 28) * 6 +
                        (speeds >= 28) * (speeds < 34) * 7 +
                        (speeds >= 34) * (speeds < 41) * 8 +
                        (speeds >= 41) * (speeds < 48) * 9 +
                        (speeds >= 48) * (speeds < 56) * 10 +
                        (speeds >= 56) * (speeds < 64) * 11 +
                        (speeds >= 64) * 12)
            cat_map = { 0: "Calm",
                        1: "Light air",
                        2: "Light breeze",
                        3: "Gentle breeze",
                        4: "Moderate breeze",
                        5: "Fresh breeze",
                        6: "Strong breeze",
                        7: "High wind",
                        8: "Gale",
                        9: "Strong gale",
                       10: "Whole gale",
                       11: "Violent storm",
                       12: "Hurricane"}

        elif categorization.upper() == "NHC":
            # TODO:  Change these to m/s (knots are how these are defined).
            # Definitely not in the correct format now
            # TODO:  Add TD and TS designations
            speeds = units.convert(self.max_wind_speed, "m/s", "knots")
            category = (numpy.zeros(speeds.shape) +
                        (speeds < 30) * -1 +
                        (speeds >= 64) * (speeds < 83) * 1 +
                        (speeds >= 83) * (speeds < 96) * 2 +
                        (speeds >= 96) * (speeds < 113) * 3 +
                        (speeds >= 113) * (speeds < 135) * 4 +
                        (speeds >= 135) * 5)
            cat_map = {-1: "Tropical Depression",
                        0: "Tropical Storm",
                        1: "Category 1 Hurricane",
                        2: "Category 2 Hurricane",
                        3: "Category 3 Hurricane",
                        4: "Category 4 Hurricane",
                        5: "Category 5 Hurricane"}

        elif categorization.upper() == "JTWC":
            raise NotImplementedError("JTWC categorization not implemented.")
        elif categorization.upper() == "JMA":
            raise NotImplementedError("JMA categorization not implemented.")
        elif categorization.upper() == "IMD":
            raise NotImplementedError("IMD categorization not implemented.")
        elif categorization.upper() == "MF":
            raise NotImplementedError("MF categorization not implemented.")
        elif categorization.upper() == "BOM":
            raise NotImplementedError("BOM categorization not implemented.")
        else:
            raise ValueError("Categorization %s not available."
                             % categorization)

        if cat_names:
            category_name = []
            for (i, cat) in enumerate(category):
                category_name.append(cat_map[cat])

            return category, category_name
        else:
            return category


# =============================================================================
# Model field construction - Models supported are
#  - Holland 1980 ('HOLLAND_1980') [1]
#  - Holland 2010 ('HOLLAND_2010') [2]
#  - Chavas, Lin, Emmanuel ('CLE_2015') [3]
# *TODO* - Add citations

# Dictionary of models.  Keys are function names, values are the proper name
# and a citation to the model
_supported_models = {"holland_1980": ["Holland 1980", "Holland, G. J. An Analytic Model of the Wind and Pressure Profiles in Hurricanes. Monthly Weather Review 108, 1212–1218 (1980)."],
                     "holland_2010": ["Holland 2010", "Holland, G. J., Belanger, J. I. & Fritz, A. A Revised Model for Radial Profiles of Hurricane Winds. Monthly Weather Review 138, 4393–4393 (2010)."],
                     "cle_2015": ["Chavas, Lin, Emmanuel 2015", "Chavas, D. R., Lin, N. & Emanuel, K. A Model for the Complete Radial Structure of the Tropical Cyclone Wind Field. Part I: Comparison with Observed Structure*. http://dx.doi.org.ezproxy.cul.columbia.edu/10.1175/JAS-D-15-0014.1 72, 3647–3662 (2015)."]}


# In the case where the field is not rotationally symmetric then the r value
# defines the x and y axis extents.
def construct_fields(storm, r, t, model="holland_1980"):
    r""""""

    if model.lower() not in _supported_models.keys():
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
    output = "Available Formats: (Function, Name, Citation)\n"
    for (model, values) in Storm._supported_formats.items():
        output = "".join((output, "%s: %s %s\n" % (values[0], model,
                                                   values[1])))
    return output


def available_models():
    r"""Construct a string suitable for listing available storm models.
    """
    output = "Function, Name, Citation\n"
    for (model, values) in _supported_models.items():
        output = "".join((output, "%s: %s %s\n" % (values[0], model,
                                                   values[1])))
    return output


# =============================================================================
# Ensmeble Storm Formats
def load_emmanuel_storms(path, mask_distance=None, mask_coordinate=(0.0, 0.0),
                               mask_category=None, categorization="NHC"):
    r"""Load storms from a Matlab file containing storms

    This format is based on the format Prof. Emmanuel uses to generate storms.

    :Input:
     - *path* (string) Path to the file to be read in
     - *mask_distance* (float) Distance from *mask_coordinate* at which a storm
       needs to in order to be returned in the list of storms.  If
       *mask_distance* is *None* then no masking is used.  Default is to
       use no *mask_distance*.
     - *mask_coordinate* (tuple) Longitude and latitude coordinates to measure
       the distance from.  Default is *(0.0, 0.0)*.
     - *mask_category* (int) Category or highter a storm needs to be to be
       included in the returned list of storms.  If *mask_category* is *None*
       then no masking occurs.  The categorization used is controlled by
       *categorization*.  Default is to use no *mask_category*.
     - *categorization* (string) Categorization to be used for the
       *mask_category* filter.  Default is "NHC".

    :Output:
     - (list) List of Storm objects that have been read in and were not
       filtered out.
    """

    # Load the mat file and extract pertinent data
    import scipy.io
    mat = scipy.io.loadmat(path)

    lon = mat['longstore']
    lat = mat['latstore']
    hour = mat['hourstore']
    day = mat['daystore']
    month = mat['monthstore']
    year = mat['yearstore']
    max_wind_radius = mat['rmstore']
    max_wind_speed = mat['vstore']
    central_pressure = mat['pstore']

    # Convert into storms and truncate zeros
    storms = []
    for n in range(lon.shape[0]):
        m = len(lon[n].nonzero()[0])

        storm = Storm()
        storm.t = [datetime.datetime(year[0, n],
                                     month[n, i],
                                     day[n, i],
                                     hour[n, i]) for i in range(m)]

        storm.eye_location = numpy.empty((m, 2))
        storm.max_wind_speed = numpy.empty(m)
        storm.max_wind_radius = numpy.empty(m)
        storm.central_pressure = numpy.empty(m)

        storm.eye_location[:, 0] = lon[n, :m]
        storm.eye_location[:, 1] = lat[n, :m]
        storm.max_wind_speed = max_wind_speed[n, :m]
        storm.max_wind_radius = max_wind_radius[n, :m]
        storm.central_pressure = central_pressure[n, :m]
        storm.storm_radius = numpy.ones(m) * 300e3

        include_storm = True
        if mask_distance is not None:
            distance = numpy.sqrt((storm.eye_location[:, 0] -
                                   mask_coord[0])**2 +
                                  (storm.eye_location[:, 1] -
                                   mask_coord[1])**2)
            inlcude_storm = numpy.any(distance < mask_distance)
        if mask_category is not None:
            raise NotImplementedError("Category masking not implemented.")

        if include_storm:
            storms.append(storm)

    return storms


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    # Positional argument
    parser.add_argument("path", help="Path to storm file to be read in")

    # Optional arguments
    parser.add_argument("-f", "--from", default="atcf", dest="input_format",
                        help="Format to convert from, defaults to 'atcf'")
    parser.add_argument("-o", "--output", default="geoclaw.storm",
                        dest="output_path",
                        help="Output path, default to 'geoclaw.storm'")
    parser.add_argument("-t", "--to", default="geoclaw",
                        dest="output_format",
                        help="Format to convert to, defaults to 'geoclaw'")
    parser.add_argument("-v", "--verbose",
                        help="Increase verbosity of output",
                        action="store_true")

    args = parser.parse_args()
    input_storm = Storm(args.path, file_format=args.input_format,
                        verbose=args.verbose)
    input_storm.write(args.output_path, file_format=args.output_format,
                      verbose=args.verbose)
