#!/usr/bin/env python

r"""
GeoClaw util Module  `$CLAW/geoclaw/src/python/geoclaw/util.py`

Module provides provides utility functions.


:Functions:

 - dms2decimal - Convert (degrees, minutes, seconds) to decimal degrees
 - dist_meters2latlong - Convert dx, dy distance in meters to degrees
 - dist_latlong2meters - Convert dx, dy distance in degrees to meters
 - haversine - Calculate the haversine based great circle distance
 - inv_haversine - Inverts the haversine distance
 - fetch_noaa_tide_data - Fetches water levels and tide predictions
"""

from __future__ import absolute_import
from __future__ import print_function
import io
import os
import os.path

import numpy
from six.moves.urllib.parse import urlencode
from six.moves.urllib.request import urlopen

# ==============================================================================
#  Constants
# ==============================================================================
from clawpack.geoclaw.data import Rearth, DEG2RAD, RAD2DEG, LAT2METER

NOAA_API_URL = 'https://tidesandcurrents.noaa.gov/api/datagetter'

# ==============================================================================
#  Functions for calculating Distances
# ==============================================================================
def dms2decimal(d,m,s,coord='N'):
    r"""Convert coordinates in (degrees, minutes, seconds) to decimal form.  
    
    If coord == 'S' or coord == 'W' then value is negated too.

    :Example: 

        >>> topotools.dms2decimal(7,30,36,'W')
        -7.51

    (Note that you might want to add 360 to resulting W coordinate
    if using E coordinates everywhere in a computation spanning date line.)

    :Returns: float

    """

    deg = d + m / 60.0 + s / 3600.0
    if coord in ['S','W']:
        deg = -deg

    return deg


def dist_latlong2meters(dx, dy, latitude=0.0):
    """Convert distance from degrees longitude-latitude to meters.

    Takes the distance described by *dx* and *dy* in degrees and converts it into
    distances in meters.

    returns (float, float) 

    """

    dym = Rearth * DEG2RAD * dy
    dxm = Rearth * numpy.cos(latitude * DEG2RAD) * dx * DEG2RAD

    return dxm,dym


def dist_meters2latlong(dx, dy, latitude=0.0):
    """Convert distance from meters to degrees of longitude-latitude.

    Takes the distance described by *dx* and *dy* in meters and converts it into
    distances in the longitudinal and latitudinal directions in degrees.  

    returns (float, float)

    """

    dxd = dx / (Rearth * numpy.cos(latitude * DEG2RAD)) * RAD2DEG
    dyd = dy * RAD2DEG / Rearth

    return dxd, dyd


def haversine(x0, y0, x1=None, y1=None, units='degrees'):

    """
    x0,y0 is assumed to be a point (or an array with the same shapes as x1,y1)
    x1,y1 is a point or two arrays of points (of the same dimension)
    returns array with same shape as x1 and y1 containing distance of each point
    from (x0,y0).

    For backward compatibility, also allows x0,y0 to be 2-tuples specifying
    two points, but this is not suggested since the notation is not consistent.
    """
    
    if x1 is None:
        # for backward compatibility, assume in this case that x0 and y0 
        # are tuples for the two desired points:
        assert len(x0)==len(y0)==2, "*** Unexpected input"
        x1,y1 = y0
        x0,y0 = x0

    if units == 'degrees':
        # convert to radians:
        x0 = x0*DEG2RAD
        y0 = y0*DEG2RAD
        x1 = x1*DEG2RAD
        y1 = y1*DEG2RAD

    dx = x1 - x0
    dy = y1 - y0

    # angle subtended by two points, using Haversine formula:
    dsigma = 2.0 * numpy.arcsin( numpy.sqrt( numpy.sin(0.5 * dy)**2   \
            + numpy.cos(y0) * numpy.cos(y1) * numpy.sin(0.5 * dx)**2))

    return Rearth * dsigma


def inv_haversine(d,x1,y1,y2,Rsphere=Rearth,units='degrees'):
    r"""Invert the Haversine function to find dx given a distance and point.


    Invert the haversine function to find dx given distance d and (x1,y1) and y2.
    The corresponding x2 can be x1+dx or x1-dx.
    May return NaN if no solution.
    """

    if units=='degrees':
        # convert to radians:
        x1 = x1 * RAD2DEG
        y1 = y1 * RAD2DEG
        y2 = y2 * RAD2DEG
    elif units != 'radians':
        raise Exception("unrecognized units")
    dsigma = d / Rsphere
    cos_dsigma = (numpy.cos(dsigma) - numpy.sin(y1)*numpy.sin(y2)) / (numpy.cos(y1)*numpy.cos(y2))
    dx = numpy.arccos(cos_dsigma)
    if units=='degrees':
        dx = dx * RAD2DEG
    return dx


def fetch_noaa_tide_data(station, begin_date, end_date, time_zone='GMT',
                         datum='STND', units='metric', cache_dir=None,
                         verbose=True):
    """Fetch water levels and tide predictions at given NOAA tide station.

    The data is returned in 6 minute intervals between the specified begin and
    end dates/times.  A complete specification of the NOAA CO-OPS API for Data
    Retrieval used to fetch the data can be found at:

        https://tidesandcurrents.noaa.gov/api/

    By default, retrieved data is cached in the geoclaw scratch directory
    located at:

        $CLAW/geoclaw/scratch

    :Required Arguments:
      - station (string): 7 character station ID
      - begin_date (datetime): start of date/time range of retrieval
      - end_date (datetime): end of date/time range of retrieval

    :Optional Arguments:
      - time_zone (string): see NOAA API documentation for possible values
      - datum (string): see NOAA API documentation for possible values
      - units (string): see NOAA API documentation for possible values
      - cache_dir (string): alternative directory to use for caching data
      - verbose (bool): whether to output informational messages

    :Returns:
      - date_time (numpy.ndarray): times corresponding to retrieved data
      - water_level (numpy.ndarray): preliminary or verified water levels
      - prediction (numpy.ndarray): tide predictions
    """
    # use geoclaw scratch directory for caching by default
    if cache_dir is None:
        if 'CLAW' not in os.environ:
            raise ValueError('CLAW environment variable not set')
        claw_dir = os.environ['CLAW']
        cache_dir = os.path.join(claw_dir, 'geoclaw', 'scratch')

    def fetch(product, expected_header, col_idx, col_types):
        noaa_params = get_noaa_params(product)
        cache_path = get_cache_path(product)

        # use cached data if available
        if os.path.exists(cache_path):
            if verbose:
                print('Using cached {} data for station {}'.format(
                    product, station))
            return parse(cache_path, col_idx, col_types, header=True)

        # otherwise, retrieve data from NOAA and cache it
        if verbose:
            print('Fetching {} data from NOAA for station {}'.format(
                product, station))
        full_url = '{}?{}'.format(NOAA_API_URL, urlencode(noaa_params))
        with urlopen(full_url) as response:
            text = response.read().decode('utf-8')
            with io.StringIO(text) as data:
                # ensure that received header is correct
                header = data.readline().strip()
                if header != expected_header or 'Error' in text:
                    # if not, response contains error message
                    raise ValueError(text)

                # if there were no errors, then cache response
                save_to_cache(cache_path, text)

                return parse(data, col_idx, col_types, header=False)

    def get_noaa_params(product):
        noaa_date_fmt = '%Y%m%d %H:%M'
        noaa_params = {
            'product': product,
            'application': 'NOS.COOPS.TAC.WL',
            'format': 'csv',
            'station': station,
            'begin_date': begin_date.strftime(noaa_date_fmt),
            'end_date': end_date.strftime(noaa_date_fmt),
            'time_zone': time_zone,
            'datum': datum,
            'units': units
        }
        return noaa_params

    def get_cache_path(product):
        cache_date_fmt = '%Y%m%d%H%M'
        dates = '{}_{}'.format(begin_date.strftime(cache_date_fmt),
                               end_date.strftime(cache_date_fmt))
        filename = '{}_{}_{}'.format(time_zone, datum, units)
        abs_cache_dir = os.path.abspath(cache_dir)
        return os.path.join(abs_cache_dir, product, station, dates, filename)

    def save_to_cache(cache_path, data):
        # make parent directories if they do not exist
        parent_dir = os.path.dirname(cache_path)
        if not os.path.exists(parent_dir):
            os.makedirs(parent_dir)

        # write data to cache file
        with open(cache_path, 'w') as cache_file:
            cache_file.write(data)

    def parse(data, col_idx, col_types, header):
        # read data into structured array, skipping header row if present
        a = numpy.genfromtxt(data, usecols=col_idx, dtype=col_types,
                             skip_header=int(header), delimiter=',',
                             missing_values='')

        # return tuple of columns
        return tuple(a[col] for col in a.dtype.names)

    # only need first two columns of data; first column contains date/time,
    # and second column contains corresponding value
    col_idx = (0, 1)
    col_types = 'datetime64[m], float'

    # fetch water levels and tide predictions
    date_time, water_level = fetch(
        'water_level', 'Date Time, Water Level, Sigma, O, F, R, L, Quality',
        col_idx, col_types)
    date_time2, prediction = fetch('predictions', 'Date Time, Prediction',
                                   col_idx, col_types)

    # ensure that date/time ranges are the same
    if not numpy.array_equal(date_time, date_time2):
        raise ValueError('Received data for different times')

    return date_time, water_level, prediction
