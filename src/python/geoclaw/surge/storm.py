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

#                           days   s/hour   hour/day                             
days2seconds = lambda days: days * 60.0**2 * 24.0

import clawpack.geoclaw.units as units

# Define supported formats and models
_supported_formats = ["GEOCLAW", "HURDAT", "HURDAT2", "JMA", "IMD", "TCVITALS"]
#_supported_formats = ["GEOCLAW"]
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
        self.new_year = None  

        if path is None:
            self.read(path, file_format=file_format, **kwargs)
        else: 
            self.read(path, file_format=file_format, **kwargs)  

    # =========================================================================
    # Read Routines
 
    def read(self, path, file_format, **kwargs):
#        raise ValueError("File format %s not available." % file_format)

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
        self.eye_location[0, :] = data[:, 1]
        self.eye_location[1, :] = data[:, 2]
        self.max_wind_speed = data[:, 3]
        self.max_wind_radius = data[:, 4]
        self.central_pressure = data[:, 5]
        self.storm_radius = data[:, 6]

    def read_hurdat(self, path, file_format):
        r"""Read in a Hurdat formatted storm file and 
        extract relevant storm fields. 

        :Input:
         - *path* (string) Path to the file to be read.
         - file_format (string) File format to be used. 
        
        """
        print('file_format', file_format) 
        if file_format.upper() not in _supported_formats:
            raise ValueError("File type not one of supported formats.")
        else:
            with open(path,'r') as data_file:
                # Collect data from columns of the same type 
                data = numpy.genfromtxt(path,delimiter=',',dtype=None,usecols=(8,9,18,19))
                self.max_wind_speed = data[:,0]                # Col 8 Hurdat 
                self.central_pressure = data[:,1]              # Col 9 Hurdat 
                self.storm_radius = data[:,2]                  # Col 18 Hurdat
                self.max_wind_radius = data[:,3]               # Col 19 Hurdat 

                # Convert Col 2 Hurdat Date into seconds 
                date = numpy.genfromtxt(path,delimiter=',',dtype=None,usecols=(2)) 
                self.new_year = int(str(date[0])[0:4])  
                for i in range(date.shape[0]):                 # Convert date into seconds 
                    date[i] = self.date2seconds(str(date[i]))
                self.t = date
                
                self.eye_location = numpy.genfromtxt(path,dtype=None,usecols=(6,7),delimiter=',')  
                for n in range(self.eye_location.shape[0]): 
                    lat = self.eye_location[n,0]
                    lon = self.eye_location[n,1]
                    if lat[-1] == 'N':
                        lat = float(lat[0:-1])/10 
                    else: 
                        lat = -1*float(lat[0:-1])/10
                    if lon == 'E': 
                        lon = float(lon[0:-1])/10 
                    else: 
                        lon = -1*float(lon[0:-1])/10
                    self.eye_location[n,0] = lat 
                    self.eye_location[n,1] = lon

    def read_hurdat2(self, path, file_format):
        r"""Extract relevant hurricane data from Hurdat2 file
            and update storm fields with proper values. 
        
        :Input: 
         - *path* (string) Path to the file to be read.

        Return ValueError if format incorrect or if file not Hurdat2.
        """
        if file_format.upper() not in _supported_formats:
            raise ValueError("File type not one of supported formats.")
        else:
            with open(path,'r') as data_file:
                # Collect data from columns of the same type 
                data = numpy.genfromtxt(path,delimiter=',',skip_header=1,dtype=None,usecols=(6,7,8,9))
                self.max_wind_speed = data[:,0]                # Col 6 Hurdat 
                self.central_pressure = data[:,1]              # Col 7 Hurdat 
                self.max_wind_radius = data[:,2]               # Col 8 Hurdat 
                self.storm_radius = data[:,3]                  # Col 9 Hurdat

                date = numpy.genfromtxt(path,delimiter=',',skip_header=1,dtype=str,usecols=(0,1))
                self.new_year = int(str(date[0][0])[0:4])  
                for i in range(date.shape[0]): 
                    temp_date = "%s%s" %(date[i][0],date[i][1][0:2])
                    temp_date = self.date2seconds(temp_date)
                    date[i][0] = temp_date 
                self.t = date[:,0] 

                self.eye_location = numpy.genfromtxt(path,delimiter=',',skip_header=1,dtype=None,usecols=(4,5))  
                for n in range(self.eye_location.shape[0]): 
                    lat = self.eye_location[n,0]
                    lon = self.eye_location[n,1]
                    if lat[-1] == 'N':
                        lat = float(lat[0:-1]) 
                    else: 
                        lat = -1*float(lat[0:-1])
                    if lon == 'E': 
                        lon = float(lon[0:-1]) 
                    else: 
                        lon = -1*float(lon[0:-1])
                    self.eye_location[n,0] = lat 
                    self.eye_location[n,1] = lon
                
    def read_jma(self, path, file_format):
        r"""Extract relevant hurricane data from JMA file
            and update storm fields with proper values. 
        
        :Input: 
         - *path* (string) Path to the file to be read.

        Return ValueError if format incorrect or if file not JMA.
        """
        if file_format.upper() not in _supported_formats:
            raise ValueError("File type not one of supported formats.")
        else:
            with open(path,'r') as data_file:
                # Collect data from columns of the same type 
                data = numpy.genfromtxt(path,skip_header=1,dtype=None,usecols=(5,6,7,9))
                self.central_pressure = data[:,0]             # Col 5 JMA  
                self.max_wind_speed  = data[:,1]              # Col 6 JMA 
                self.max_wind_radius = data[:,2]              # Col 7 JMA 
                self.storm_radius = data[:,3]                 # Col 9 JMA 

                date = numpy.genfromtxt(path,skip_header=1,dtype=None,usecols=(0))
                self.new_year = int('20'+str(date[0])[0:2])  
                for i in range(date.shape[0]): 
                    temp_date = "%s%s" %('20',date[i][0])
                    temp_date = self.date2seconds(temp_date[0:-2])
                    date[i][0] = temp_date 
                self.t = date[:,0] 
                self.eye_location = numpy.genfromtxt(path,skip_header=1,dtype=None,usecols=(4,5))  
                for n in range(self.eye_location.shape[0]):
					lat = self.eye_location[n,0]
					lon = self.eye_location[n,1] 
					lat = float(lat)*0.1
					lon = float(lon)*0.1
					self.eye_location[n,0] = lat 
					self.eye_location[n,1] = lon  
 
    def read_imd(self, path):
        r"""Extract relevant hurricane data from IMD file
            and update storm fields with proper values. 
        
        :Input: 
         - *path* (string) Path to the file to be read.

        Return ValueError if format incorrect or if file not IMD.
        """
        raise ValueError("File type not implemented yet.") 
 
    def read_tcvitals(self, path):
        r"""Extract relevant hurricane data from TCVITALS file
            and update storm fields with proper values. 
        
        :Input: 
         - *path* (string) Path to the file to be read.

        Return ValueError if format incorrect or if file not TCVITALS.
        """
        with open(path,'r') as data_file:
            # Collect data from columns of the same type 
            data = numpy.genfromtxt(path,dtype=None,usecols=(8,9,11,13))
            self.max_wind_speed = data[:,0]                # Col 8 TCVITALS
            self.central_pressure = data[:,1]              # Col 9 TCVITALS
            self.storm_radius = data[:,2]                  # Col 11 TCVITALS
            self.max_wind_radius = data[:,3]               # Col 13 TCVITALS

            yymmdd_date = numpy.genfromtxt(path,dtype='string',usecols=(3))
            hhmm_date = numpy.genfromtxt(path,dtype='string',usecols=(4))
            self.new_year = int(str(yymmdd_date[0])[0:4]) 
            for i in range(yymmdd_date.shape[0]): 
                temp_date = "%s%s" %(yymmdd_date[i],hhmm_date[i])
                temp_date = self.date2seconds(temp_date[0:-2])
                yymmdd_date[i] = temp_date 
            self.t = yymmdd_date 

            self.eye_location = numpy.genfromtxt(path,dtype=None,usecols=(5,6))  
            for n in range(self.eye_location.shape[0]): 
                lat = self.eye_location[n,0]
                lon = self.eye_location[n,1]
                if lat[-1] == 'N':
                    lat = float(lat[0:-1])/10 
                else: 
                    lat = -1*float(lat[0:-1])/10
                if lon == 'E': 
                    lon = float(lon[0:-1])/10
                else: 
                    lon = -1*float(lon[0:-1])/10
                self.eye_location[n,0] = lat 
                self.eye_location[n,1] = lon



    # =========================================================================
    # Write Routines
    def write(self, path, file_format="geoclaw"):
        r"""Take in the file_format and employ the appropriate 
		write function. 
		"""
        if file_format.upper() not in _supported_formats:
            raise ValueError("File format %s not available." % file_format)
        else:
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
                data_file.write("%s %s %s %s %s %s %s %s" %
                                                (self.t[n],
                                                 self.eye_location[n, 0],
                                                 self.eye_location[n, 1],
                                                 self.max_wind_speed[n],
                                                 self.max_wind_radius[n],
                                                 self.central_pressure[n],
                                                 self.storm_radius[n], 
                                                 "\n"))

    def write_hurdat(self, path):
#        r"""Rewrite the storm in the format 
#        used by Hurdat 
#        """
#        with open(path, 'w') as data_file:
#            for n in range(self.t.shape[0]):
#                data_file.write(", , %s, , , , %s, %s, %s, %s, , , , , , , , , %s, %s, , , , , , , , , , \n" %
#                                                 (str(self.seconds2date(self.t[n])),
#                                                 str(self.return_lat(float(self.eye_location[n, 0])*10)),
#                                                 str(self.return_lon(float(self.eye_location[n, 1])*10)),
#                                                 str(self.max_wind_speed[n]),
#                                                 str(self.central_pressure[n]),
#                                                 str(self.storm_radius[n]), 
#                                                 str(self.max_wind_radius[n])
#                               ))
        r"""Rewrite the storm in the format 
        used by Hurdat 
        """
        with open(path, 'w') as data_file:
            for n in range(self.t.shape[0]):
                #data_file.write("%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s\n" %
                data_file.write("%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s\n" %
                                                 (" ".rjust(8),
                                                str(self.seconds2date(self.t[n])).rjust(10),
                                                " ".rjust(6), 
                                                "BEST".rjust(4), 
                                                " ".rjust(2), 
                                                "000".rjust(3), 
                                                " ".rjust(1), 
                                                str(self.return_lat(float(self.eye_location[n, 0])*10)).rjust(5),
                                                " ".rjust(2), 
                                                str(self.return_lon(float(self.eye_location[n, 1])*10)).rjust(5),
                                                " ".rjust(2), 
                                                str(self.max_wind_speed[n]).rjust(3),
                                                " ".rjust(2), 
                                                str(self.central_pressure[n]).rjust(4),
                                                " ".rjust(47), 
                                                str(self.storm_radius[n]).rjust(3), 
                                                " ".rjust(2), 
                                                str(self.max_wind_radius[n]).rjust(3)
                               ))

    def write_hurdat2(self, path):
        r"""Rewrite the storm in the format 
        used by Hurdat2  
        """
        with open(path, 'w') as data_file:
            data_file.write("%s,    %s,    %s" %
                             (self.seconds2date(self.t[n]),
                             'Hurricane Name', 
                             'Indicator'))
            for n in range(self.t.shape[0]):
                if self.return_lat < 0: 
                    latitude_direction = 'S'
                else: 
                    latitude_direction = 'N'
                if self.return_lon < 0: 
                    longitude_direction = 'W'
                else: 
                    longitude_direction = 'E' 
                data_file.write("%s, %s00, , ,%s%s,%s%s , %s, %s, %s, %s, , , , , , , , , , \n" %
                                                 (self.seconds2date(self.t[n])[0:-2],
											     self.seconds2date(self.t[n])[-2], 
                                                 self.return_lat(numpy.abs(float(self.eye_location[n, 0])),
											     latitude_direction, 
                                                 self.return_lon(numpy.abs(float(skelf.eye_location[n, 1]))),
												 longitude_direction, 
                                                 self.max_wind_speed[n],
                                                 self.central_pressure[n],
                                                 self.storm_radius[n], 
                                                 self.max_wind_radius[n]
                               ))) 
    

    def write_jma(self, path):
        r"""Rewrite the storm in the format 
		used by JMA. 
		"""
        with open(path, 'w') as data_file:
            for n in range(self.t.shape[0]):
                data_file.write("%s    %s %s   %s %s  %s %s %s" %
                               					(self.seconds2date(self.t[n]),
                                                 self.return_lat(float(self.eye_location[n, 0])*10),
                                                 self.return_lon(float(self.eye_location[n, 1])*10),
                                                 self.max_wind_speed[n],
                                                 self.central_pressure[n],
                                                 self.storm_radius[n], 
                                                 self.max_wind_radius[n], 
												 '\n'))

    def write_imd(self, path):
        r""""""
        raise NotImplementedError("IMD format not fully implemented.")

    # =========================================================================
    # =========================================================================
    # TO BE DELETED SOME TIME SOON  
    
    def date2seconds(self,date): 
        r"""Helper function to transform dates into seconds. 

        Input: Date in format YYYYMMDDHH
        
        Output: Difference between date and beginning of the year 
                returned in units of seconds. 
        """
        # Parse data to collect year, month, day, and hour as integers 
        year = int(date[0:4])
        month = int(date[4:6])
        day = int(date[6:8]) 
        hour = int(date[8:])

        new_year_day = datetime.datetime(year,1,1,0) # Determine first day of year  
        delta_date = datetime.datetime(year,month,day,hour) - new_year_day
        date_in_seconds = days2seconds(delta_date.days) + delta_date.seconds
        return date_in_seconds 
    
    def seconds2date(self,seconds): 
        r"""Helper function to transform second into appropriate date.
        
        Input: Date in seconds 
        
        Output: Date as YYYYMMDDHH 
        """
        new_year_day = datetime.datetime(self.new_year,1,1,0)
        seconds = int(seconds) 
        days = seconds/86400
        secs = seconds%86400
        delta_date = datetime.timedelta(days,secs)
        date=new_year_day+delta_date
        year = str(date.year).rjust(4,'0') 
        month = str(date.month).rjust(2,'0') 
        day = str(date.day).rjust(2,'0') 
        hour = str(date.hour).rjust(2,'0') 
        date_to_string = "%s%s%s%s" %(year.rjust(4),month.rjust(2),day.rjust(2),hour.rjust(2))  
        return date_to_string
    
    def return_lat(self, lat): 
        r"""Helper function to transfrom latitude.
        """
        if lat > 0: 
            lat=str(int(lat))+'N'
        else:
            lat = lat*-1 
            lat=str(int(lat))+'S'
        return lat
  
    def return_lon(self,lon): 
        r"""Helper function to transform longitude.
        """
        if lon > 0: 
            lon=str(int(lon))+'E'
        else:
            lon = lon*-1  
            lon=str(int(lon))+'W'
        return lon
     
    # =========================================================================
    # =========================================================================
    
    # =========================================================================
    # Other Useful Routines
    def plot(self, axes=None, intensity=False, limits=None, track_color='red',
                   category_colors=None, categorization="NHC"):
        r"""Plot the track and optionally the strength of the storm

        """

        import matplotlib.pyplot as plt
        from mpl_toolkits.basemap import Basemap

        if axes is None:
            fig = plt.figure()
            axes = fig.add_subplot(1, 1, 1)

        #limits = ((long), (lat))
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
        # Not sure how to do this automatically yet
        #mapping.drawparallels(limits[])
        # mapping.drawparallels((0.0, 20.0), labels=[1, 1])
        # mapping.drawmeridians(numpy.arange(coord[0][0], coord[1][0], 20),
        #                       labels=[0, 1, 1, 1])

        return axes

    def category(self, categorization="NHC"):
        r"""Categorizes storm based on relevant storm data

        :Category Mappings:
         - "NHC":  -1 (Tropical Depression),
                    0 (Tropical Storm),
                    1 (Category 1 Hurricane),
                    2 (Category 2 Hurricane),
                    3 (Category 3 Hurricane),
                    4 (Category 4 Hurricane),
                    5 (Category 5 Hurricane)
         - "JMA":
         - "IMD":

        :Input:
         - *categorization* (string) Type of categorization to use.  Defaults to
           the National Hurricane Center "NHC".
         - *names* (bool) If True returns the category name rather than a
           number.  Default to *False*.

        :Output:
         - (ndarray) Integer array of categories at each time point of the storm
         - (list) Similar to the above but the name of the category as a 
           *string*

        """

        if categorization.upper() == "NHC":
            # TODO:  Check to see if these are in knots or mph.  Definitely not
            #        in the correct format now
            # TODO:  Add TD and TS designations
            category = numpy.zeros(self.max_wind_speed) + \
                       (self.max_wind_speed >= 64) * (self.max_wind_speed < 83) * 1 + \
                       (self.max_wind_speed >= 83) * (self.max_wind_speed < 96) * 2 + \
                       (self.max_wind_speed >= 96) * (self.max_wind_speed < 113) * 3 + \
                       (self.max_wind_speed >= 113) * (self.max_wind_speed < 135) * 4 + \
                       (self.max_wind_speed >= 135) * 5
            cat_map = {-1: "Tropical Depression",
                        0: "Tropical Storm",
                        1: "Category 1 Hurricane",
                        2: "Category 2 Hurricane",
                        3: "Category 3 Hurricane",
                        4: "Category 4 Hurricane",
                        5: "Category 5 Hurricane"}

        elif categorization.upper() == "JMA":
            raise NotImplementedError("JMA categorization not implemented.")
        elif categorization.upper() == "IMD":
            raise NotImplementedError("IMD categorization not implemented.")
        else:
            raise ValueError("Categorization %s not available."
                             % categorization)

        if names:
            category_name = []
            for (i, cat) in enumerate(category):
                category_name.append(cat_map[cat])

        return category, category_name


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


# =============================================================================
# Ensmeble Storm Formats
#def load_emmanuel_storms(path, mask_distance=None, mask_coordinate=(0.0, 0.0),
#                               mask_category=None, categorization="NHC"):
#    r"""Load storms from a Matlab file containing storms
#
#    This format is based on the format Prof. Emmanuel uses to generate storms.
#
#    :Input:
#     - *path* (string) Path to the file to be read in
#     - *mask_distance* (float) Distance from *mask_coordinate* at which a storm
#       needs to in order to be returned in the list of storms.  If
#       *mask_distance* is *None* then no masking is used.  Default is to
#       use no *mask_distance*.
#     - *mask_coordinate* (tuple) Longitude and latitude coordinates to measure
#       the distance from.  Default is *(0.0, 0.0)*.
#     - *mask_category* (int) Category or highter a storm needs to be to be
#       included in the returned list of storms.  If *mask_category* is *None*
#       then no masking occurs.  The categorization used is controlled by
#       *categorization*.  Default is to use no *mask_category*.
#     - *categorization* (string) Categorization to be used for the
#       *mask_category* filter.  Default is "NHC".
#
#    :Output:
#     - (list) List of Storm objects that have been read in and were not filtered
#       out.
#    """
#
#    # Load the mat file and extract pertinent data
#    import scipy.io
#    mat = scipy.io.loadmat(path)
#
#    lon = mat['longstore']
#    lat = mat['latstore']
#    hour = mat['hourstore']
#    day = mat['daystore']
#    month = mat['monthstore']
#    year = mat['yearstore']
#    radius_max_winds = mat['rmstore']
#    max_winds = mat['vstore']
#    central_pressure = mat['pstore']
#
#    # Convert into storms and truncate zeros
#    storms = []
#    for n in xrange(lon.shape[0]):
#        m = len(lon[n].nonzero()[0])
#
#        storm = Storm()
#        storm.t = [datetime.datetime(year[0, n],
#                                     month[n, i],
#                                     day[n, i],
#                                     hour[n, i]) for i in xrange(m)]
#        storm.eye_location[:, 0] = lon[n, :m]
#        storm.eye_location[:, 1] = lat[n, :m]
#        storm.max_wind_speed = max_winds[n, :m]
#        storm.radius_max_winds = radius_max_winds[n, :m]
#        storm.central_pressure = central_pressure[n, :m]
#
#        include_storm = True
#        if mask_distance is not None:
#            distance = numpy.sqrt((storm.eye_location[:, 0] -
#                                   mask_coord[0])**2 +
#                                  (storm.eye_location[:, 1] -
#                                   mask_coord[1])**2)
#            inlcude_storm = numpy.any(distance < mask_distance)
#        if mask_category is not None:
#            include storm = include_storm and numpy.any(
#                  storm.category(categorization=categorization) > mask_category)
#
#        if include_storm:
#            storms.append(storm)
#
#    return storms


if __name__ == '__main__':
    # TODO:  Add commandline ability to convert between formats
    #construct_fields(None, None, None)
    base = Storm('al092012-tcvitals-arch.dat', 'tcvitals') 
    base.write_hurdat('issac.storm')  
