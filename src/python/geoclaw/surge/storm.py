#!/usr/bin/env python

r"""
Module defines a class and routines for managing storm best-track type input and
testing reconstructed wind and pressure fields.  Additionally some support for
ensembles of storms from various providers is also included.

The primary class of interest in the module is the `Storm` class that
facilitates dealing with various best-track formats found around the world and
the expected GeoClaw storm format that is read into the FORTRAN code.  The basic
workflow in a `setrun.py` file would do the following:

1. Create a `Storm` object by reading in from a file::

    storm = clawpack.geoclaw.surge.storm.Storm("my_storm.txt", file_format='ATCF')

2. Write out the storm object created into the GeoClaw format::

    storm.write("my_geoclaw_storm.txt", file_format="geoclaw")

3. Specify the path to the GeoClaw formatted storm file, in this case
   "my_geoclaw_storm.txt".

:Formats Supported:
    - GeoClaw (fully)
    - ATCF (reading only)
    - HURDAT (reading only)
    - IBTrACS (reading only)
    - JMA (reading only)
    - IMD (planned)
    - tcvitals (reading only)
"""

from __future__ import print_function
from __future__ import absolute_import
from six.moves import range

import warnings
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

# TCVitals basins with their expanded names
# see http://www.emc.ncep.noaa.gov/HWRF/tcvitals-draft.html
TCVitals_Basins = {"L": "North Atlantic",
                   "E": "North East Pacific",
                   "C": "North Central Pacific",
 		   "W": "North West Pacific",
		   "B": "Bay of Bengal (North Indian Ocean)",
		   "A": "Arabian Sea (North Indian Ocean)",
		   "Q": "South Atlantic",
		   "P": "South Pacific",
		   "S": "South Indian Ocean"}

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


# Warning for formats that have yet to have a default way to determine crticial
# radii from the input data
missing_data_warning_str = """*** Cannot yet automatically determine the
    maximum wind radius.  Will write out GeoClaw
    formats but note that these will not work
    when running GeoClaw currently."""


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
     - *max_wind_speed* (ndarray(:)) Maximum wind speed.  Default units are
       meters/second.
     - *max_wind_radius* (ndarray(:)) Radius at which the maximum wind speed
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
                          "ibtracs": ["IBTrACS", "ftp://filsrv.cicsnc.org/kknapp/ibtracs/testing/hotel1/provisional"],
                          "jma": ["JMA", "http://www.jma.go.jp/jma/jma-eng/jma-center/rsmc-hp-pub-eg/Besttracks/e_format_bst.html"],
                          "imd": ["IMD", "http://www.rsmcnewdelhi.imd.gov.in/index.php"],
                          "tcvitals": ["TC-Vitals", "http://www.emc.ncep.noaa.gov/mmb/data_processing/tcvitals_description.htm"]}

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
         - *file_format* (string) Format of the data file.  See list of
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
                        " - IBTrACS - ftp://filsrv.cicsnc.org/kknapp/ibtracs/testing/hotel1/provisional"
                        " - JMA - http://www.jma.go.jp/jma/jma-eng/jma-center/rsmc-hp-pub-eg/besttrack.html",
                        " - IMD - http://www.rsmcnewdelhi.imd.gov.in/index.php",
                        " - TCVITALS - http://www.emc.ncep.noaa.gov/mmb/data_processing/tcvitals_description.htm")
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
         - *verbose* (bool) Output more info regarding reading.
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
         - *verbose* (bool) Output more info regarding reading.
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
        self.classification = numpy.empty(num_lines, dtype=str)
        self.eye_location = numpy.empty((num_lines, 2))
        self.max_wind_speed = numpy.empty(num_lines)
        self.central_pressure = numpy.empty(num_lines)
        self.max_wind_radius = numpy.empty(num_lines)
        self.storm_radius = numpy.empty(num_lines)

        for (i, data) in enumerate(data_block):
            # End at an empty lines - skips lines at the bottom of a file
            if len(data) == 0:
                break

            # Grab data regarding basin and cyclone number if we are starting
            if i == 0:
                self.basin = ATCF_basins[data[0]]
                self.ID = int(data[1])

            # Create time
            self.t.append(datetime.datetime(int(data[2][:4]),
                                            int(data[2][4:6]),
                                            int(data[2][6:8]),
                                            int(data[2][-2:])))

            # # If an event is occuring record it.  If landfall then use as an
            # # offset.   Note that if there are multiple landfalls the last one
            # # is used as the offset
            # if not short_data:
            #     if len(data[22].strip()) > 0:
            #         self.event[i] = data[22].strip()
            #         if self.event[i].upper() == "L":
            #             self.time_offset = self.t[i]

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
            self.max_wind_speed[i] = units.convert(float(data[8]), 'knots', 'm/s')
            self.central_pressure[i] = units.convert(float(data[9]), 'mbar', 'Pa')

            # Mark if this is a shortened line - does not contain max wind
            # radius and outer storm radius - set those to -1 to mark them as
            # missing
            if len(data) < 19:
                self.storm_radius[i] = -1
                self.max_wind_radius[i] = -1
            else:
                self.storm_radius[i] = units.convert(float(data[18]), 'nmi', 'm')
                self.max_wind_radius[i] = units.convert(float(data[19]), 'nmi', 'm')

    def read_hurdat(self, path, verbose=False):
        r"""Read in HURDAT formatted storm file

        This is the current version of HURDAT data available (HURDAT 2).  Note
        that this assumes there is only one storm in the file (includes the
        header information though).  Future features will be added that will allow for
        a file to be read with multiple storms defined.

        For more details on the HURDAT format and getting data see

        http://www.aoml.noaa.gov/hrd/hurdat/Data_Storm.html

        :Input:
         - *path* (string) Path to the file to be read.
         - *verbose* (bool) Output more info regarding reading.

        :Raises:
         - *ValueError* If the method cannot find the name/year matching the
           storm or they are not provided when *single_storm == False* then a
           value error is risen.
        """

        with open(path, 'r') as hurdat_file:
            # Extract header
            data = [value.strip() for value in
                    hurdat_file.readline().split(',')]
            self.basin = data[0][:2]
            self.name = data[1]
            self.ID = data[2]

            # Store rest of data
            data_block = hurdat_file.readlines()

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
                self.eye_location[i, 1] = float(data[4][0:-1])
            else:
                self.eye_location[i, 1] = -float(data[4][0:-1])
            if data[5][-1] == "E":
                self.eye_location[i, 0] = float(data[5][0:-1])
            else:
                self.eye_location[i, 0] = -float(data[5][0:-1])

            # Intensity information - radii are not included directly in this
            # format and instead radii of winds above a threshold are included
            self.max_wind_speed[i] = units.convert(float(data[6]), 'knots', 'm/s')
            self.central_pressure[i] = units.convert(float(data[7]), 'mbar', 'Pa')
            warnings.warn(missing_data_warning_str)
            self.max_wind_radius[i] = -1
            self.storm_radius[i] = -1

    def read_ibtracs(self, path, storm_name, year):
        r"""Read in IBTrACS formatted storm file

        This reads in the netcdf-formatted IBTrACS v4 BETA data. The .nc
        file passed as *path* must contain a storm matching *storm_name*
        and *year*. This function will be updated, if needed, once the BETA
        version becomes an operational release.

        NOTE: Thus far, only the reading of hurdat/atcf-based best tracks (i.e. USA
        tracks) is supported.

        TODO: account for data formats from other reporting agencies

        For more details on the IBTrACS v4 BETA format and getting data see

        ftp://filsrv.cicsnc.org/kknapp/ibtracs/testing/hotel1/provisional/

        :Input:
         - *path* (string) Path to the file to be read.
         - *storm_name* (string) name of storm of interest (NAME field in IBTrACS).
         - *year* (int) year of storm of interest

        :Raises:
         - *ValueError* If the method cannot find the name/year matching the
           storm then a value error is risen.
        """

        # imports that you don't need for other read functions
        try:
            import xarray as xr
            from pandas import to_datetime
        except ImportError as e:
            print("IBTrACS currently requires xarray and pandas to work.")
            raise e

        storm_name = storm_name.upper()
        with xr.open_dataset(path,drop_variables=['time']) as ds:

            ## SLICE IBTRACS DATASET

            # match on storm-name and year
            storm_match = (ds.name == storm_name.encode())
            dts = xr.DataArray.from_series(to_datetime(ds.iso_time.astype(str).to_series())).values
            ds.iso_time.values = dts
            year_match = (ds.iso_time.dt.year == year).any(dim='time')
            ds = ds.sel(storm=(year_match & storm_match)).squeeze()
            # make sure
            if ('storm' in ds.dims.keys()) and (ds.storm.shape[0] == 0):
                raise ValueError('Storm/year not found in provided file')

            # include only valid time points for this storm
            # i.e. when we have max wind values
            ds = ds.sel(time=(ds.wmo_wind>=0))


            ## CONVERT TO GEOCLAW FORMAT

            # assign basin to be the basin where track originates
            # in case track moves across basins
            self.basin = ds.basin.values[0].astype(str)
            self.name = storm_name
            self.ID = ds.sid.astype(str).item()

            # convert datetime64 to datetime.datetime
            self.t = []
            for d in ds.iso_time:
                t = d.dt
                self.t.append(datetime.datetime(t.year,t.month,t.day,t.hour,t.minute,t.second))

            ## events
            self.event = ds.usa_record.values.astype(str)

            # time offset
            self.time_offset = numpy.array(self.t)[self.event=='L'][-1]

            # Classification, note that this is not the category of the storm
            self.classification = ds.usa_status.values
            self.eye_location = numpy.array([ds.lon,ds.lat]).T

            # Intensity information - for now, including only common, basic intensity
            # info.
            # TODO: add more detailed info for storms that have it
            self.max_wind_speed = units.convert(ds.wmo_wind,'knots','m/s').values
            self.central_pressure = units.convert(ds.wmo_pres,'mbar','Pa').values
            self.max_wind_radius = numpy.where(ds.usa_rmw >= 0,
                units.convert(ds.usa_rmw,'nmi','m'),-1)
            self.storm_radius = numpy.where(ds.usa_roci >=0,
                units.convert(ds.usa_roci.values,'nmi','m'),-1)


    def read_jma(self, path, verbose=False):
        r"""Read in JMA formatted storm file

        Note that only files that contain one storm are currently supported.

        For more details on the JMA format and getting data see

        http://www.jma.go.jp/jma/jma-eng/jma-center/rsmc-hp-pub-eg/Besttracks/e_format_bst.html

        :Input:
         - *path* (string) Path to the file to be read.
         - *verbose* (bool) Output more info regarding reading.

        :Raises:
         - *ValueError* If the method cannot find the name/year matching the
           storm or they are not provided when *single_storm == False* then a
           value error is risen.
        """

        data_block = []
        with open(path, 'r') as JMA_file:
            # Extract header
            data = JMA_file.readline()
            self.ID = data[6:10]
            num_lines = int(data[12:14])
            self.name = data[30:51].strip()

            data_block = JMA_file.readlines()
        assert(num_lines == len(data_block))

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
            data = [value.strip() for value in line.split()]

            # Create time
            self.t.append(datetime.datetime(int(data[0][:2]),
                                            int(data[0][2:4]),
                                            int(data[0][4:6]),
                                            int(data[0][6:])))

            # Classification, note that this is not the category of the storm
            self.classification[i] = int(data[1])

            # Parse eye location - Always N latitude and E longitude
            self.eye_location[i, 0] = float(data[4]) / 10.0
            self.eye_location[i, 1] = float(data[3]) / 10.0

            # Intensity information - current the radii are not directly given
            # Available data includes max/min of radius of winds of 50 and
            # 30 kts instead
            self.central_pressure[i] = units.convert(float(data[5]), 'hPa', 'Pa')
            self.max_wind_speed[i] = units.convert(float(data[6]), 'knots', 'm/s')
            warnings.warn(missing_data_warning_str)
            self.max_wind_radius[i] = -1
            self.storm_radius[i] = -1


    def read_imd(self, path, verbose=False):
        r"""Extract relevant hurricane data from IMD file
            and update storm fields with proper values.

        :Input:
         - *path* (string) Path to the file to be read.

        Return ValueError if format incorrect or if file not IMD.
        """
        raise NotImplementedError(("Reading in IMD files is not ",
                                   "implemented yet but is planned for a ",
                                   "future release."))


    def read_tcvitals(self, path, verbose=False):
        r"""Extract relevant hurricane data from TCVITALS file
            and update storm fields with proper values.

        :Input:
         - *path* (string) Path to the file to be read.
         - *verbose* (bool) Output more info regarding reading.

        """

        # read in TCVitals_file
        data_block = []
        with open(path, 'r') as TCVitals_file:
            data = TCVitals_file.readlines()
            for line in data:
                line = line.split()
                line = [value.strip() for value in line]
                data_block.append(line)
        num_lines = len(data_block)

        # Parse data block - convert to correct units
        # Conversions:
        #  max_wind_radius  - convert from km to m - 1000.0
        #  Central_pressure - convert from mbar to Pa - 100.0
        #  Radius of last isobar contour - convert from km to m - 1000.0
        self.t = []
        self.classification = numpy.empty(num_lines, dtype=str)
        self.eye_location = numpy.empty((num_lines, 2))
        self.max_wind_speed = numpy.empty(num_lines)
        self.central_pressure = numpy.empty(num_lines)
        self.max_wind_radius = numpy.empty(num_lines)
        self.storm_radius = numpy.empty(num_lines)

        for (i, data) in enumerate(data_block):
            # End at an empty lines - skips lines at the bottom of a file
            if len(data) == 0:
                break

            # Grab data regarding basin and cyclone number if we are starting
            if i == 0:
                self.basin = TCVitals_Basins[data[1][2:]]
                self.ID = int(data[1][:2])

            # Create time
            self.t.append(datetime.datetime(int(data[3][0:4]),
                                            int(data[3][4:6]),
                                            int(data[3][6:]),
                                            int(data[4][:2])))

            # Parse eye location - longitude/latitude order
            if data[5][-1] == 'N':
                self.eye_location[i, 1] = float(data[5][0:-1])/10.0
            else:
                self.eye_location[i, 1] = -float(data[5][0:-1])/10.0
            if data[6][-1] == "E":
                self.eye_location[i, 0] = float(data[6][0:-1])/10.0
            else:
                self.eye_location[i, 0] = -float(data[6][0:-1])/10.0

            # Intensity Information
            self.max_wind_speed[i] = float(data[12])
            self.central_pressure[i] = units.convert(float(data[9]), 'mbar', 'Pa')
            self.max_wind_radius[i] = units.convert(float(data[13]), 'km', 'm')
            self.storm_radius[i] = units.convert(float(data[11]), 'km', 'm')


    # =========================================================================
    # Write Routines
    def write(self, path, file_format="geoclaw", **kwargs):
        r"""Write out the storm data to *path* in format *file_format*

        :Input:
         - *path* (string) Path to data file.
         - *file_format* (string) Format of the data file.  See list of
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

    def write_geoclaw(self, path, verbose=False, max_wind_radius_fill=None,
                            storm_radius_fill=None):
        r"""Write out a GeoClaw formatted storm file

        GeoClaw storm files are read in by the GeoClaw Fortran code.

        :Input:
         - *path* (string) Path to the file to be written.
         - *verbose* (bool) Print out additional information when writing.
         - *max_wind_radius_fill* (func) Function that can be used to fill in
           missing data for `max_wind_radius` values.  This defaults to simply
           setting the value to -1.  The function signature should be
           `max_wind_radius(t, storm)` where t is the time of the forecast and
           `storm` is the storm object.  Note that if this or `storm_radius`
           field remains -1 that this data line will be assumed to be redundant
           and not be written out.
         - *storm_radius_fill* (func) Function that can be used to fill in
           missing data for `storm_radius` values.  This defaults to simply
           setting the value to -1.  The function signature should be
           `storm_radius(t, storm)` where t is the time of the forecast and
           `storm` is the storm object.  Note that if this or `max_wind_radius`
           field remains -1 that this data line will be assumed to be redundant
           and not be written
        """

        if max_wind_radius_fill is None:
            max_wind_radius_fill = lambda t, storm: -1
        if storm_radius_fill is None:
            storm_radius_fill = lambda t, storm: -1

        # Create list for output
        # Leave this first line blank as we need to count the actual valid lines
        # that will be left in the file below
        num_casts = 0
        data_string = [""]
        if self.time_offset is None:
            # Use the first time in sequence if not provided
            self.time_offset = self.t[0]
        data_string.append("%s\n\n" % self.time_offset.isoformat())
        for n in range(len(self.t)):
            # Remove duplicate times
            if n > 0:
                if self.t[n] == self.t[n - 1]:
                    continue

            format_string = ("{:19,.8e} " * 7)[:-1] + "\n"
            data = []
            data.append((self.t[n] - self.time_offset).total_seconds())
            data.append(self.eye_location[n, 0])
            data.append(self.eye_location[n, 1])
            data.append(self.max_wind_speed[n])
            # Allow custom function to set max wind radius if not
            # available
            if self.max_wind_radius[n] == -1:
                new_wind_radius = max_wind_radius_fill(self.t[n], self)
                if new_wind_radius == -1:
                    continue
                else:
                    data.append(new_wind_radius)
            else:
                data.append(self.max_wind_radius[n])

            data.append(self.central_pressure[n])

            # Allow custom function to set storm radius if not available
            if self.storm_radius[n] == -1:
                new_storm_radius = storm_radius_fill(self.t[n], self)
                if new_storm_radius == -1:
                    continue
                else:
                    data.append(new_storm_radius)
            else:
                data.append(self.storm_radius[n])

            data_string.append(format_string.format(*data))
            num_casts += 1


        # Write to actual file now that we know exactly how many lines it will
        # contain
        try:
            # Update number of forecasts here
            data_string[0] = "%s\n" % num_casts
            with open(path, "w") as data_file:
                for data_line in data_string:
                    data_file.write(data_line)

        except Exception as e:
            # Remove possibly partially generated file if not successful
            if os.path.exists(path):
                os.remove(path)
            raise e

    def write_atcf(self, path, verbose=False):
        r"""Write out a ATCF formatted storm file

        :Input:
         - *path* (string) Path to the file to be written.
         - *verbose* (bool) Print out additional information when writing.
        """
        raise NotImplementedError(("Writing out ATCF files is not implemented ",
                                   "yet but is planned for a future release."))
        try:
            with open(path, 'w') as data_file:
                for n in range(len(self.t)):
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
         - *path* (string) Path to the file to be written.
         - *verbose* (bool) Print out additional information when writing.
        """
        raise NotImplementedError(("Writing out hurdat files is not ",
                                   "implemented yet but is planned for a ",
                                   "future release."))
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
         - *path* (string) Path to the file to be written.
         - *verbose* (bool) Print out additional information when writing.
        """
        raise NotImplementedError(("Writing out JMA files is not implemented ",
                                   "yet but is planned for a future release."))
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
        r"""Write out an IMD formatted storm file

        :Input:
         - *path* (string) Path to the file to be written.
         - *verbose* (bool) Print out additional information when writing.
        """
        raise NotImplementedError(("Writing out IMD files is not implemented ",
                                   "yet but is planned for a future release."))

    def write_tcvitals(self, path, verbose=False):
        r"""Write out an TCVITALS formatted storm file

        :Input:
         - *path* (string) Path to the file to be written.
         - *verbose* (bool) Print out additional information when writing.
         """

        raise NotImplementedError(("Writing in TCVITALS files is not",
                                   "implemented yet but is planned for a ",
                                   "future release."))

    # =========================================================================
    # Other Useful Routines
    def plot(self, axes=None, intensity=False, limits=None, track_color='red',
                   category_colors=None, categorization="NHC"):
        r"""Plot the track and optionally the strength of the storm

        """

        # TODO:  Switch to cartopy plotting
        import matplotlib.pyplot as plt
        try:
            from mpl_toolkits.basemap import Basemap
        except ImportError as e:
            print("Plotting of storms is dependent on the Basemap package.")
            print("Given that this package has been end-of-lifed this will be")
            print("replaced by another package in the future.")
            raise e

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
_supported_models = {"holland_1980": ["Holland 1980", "Holland, G. J. An Analytic Model of the Wind and Pressure Profiles in Hurricanes. Monthly Weather Review 108, 1212-1218 (1980)."],
                     "holland_2010": ["Holland 2010", "Holland, G. J., Belanger, J. I. & Fritz, A. A Revised Model for Radial Profiles of Hurricane Winds. Monthly Weather Review 138, 4393-4393 (2010)."],
                     "cle_2015": ["Chavas, Lin, Emmanuel 2015", "Chavas, D. R., Lin, N. & Emanuel, K. A Model for the Complete Radial Structure of the Tropical Cyclone Wind Field. Part I: Comparison with Observed Structure*. http://dx.doi.org.ezproxy.cul.columbia.edu/10.1175/JAS-D-15-0014.1 72, 3647-3662 (2015)."]}


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
