#!/usr/bin/env python

r"""
GeoClaw util Module  `$CLAW/geoclaw/src/python/geoclaw/util.py`

Module provides provides utility functions.


:Functions:

 - strip_archive_extensions - strip off things like .tar or .gz
 - get_remote_file - Fetch a file from a url, e.g. a topo file
 - dms2decimal - Convert (degrees, minutes, seconds) to decimal degrees
 - dist_meters2latlong - Convert dx, dy distance in meters to degrees
 - dist_latlong2meters - Convert dx, dy distance in degrees to meters
 - haversine - Calculate the haversine based great circle distance
 - inv_haversine - Inverts the haversine distance
"""

import os
import urllib2
import tarfile
import numpy
import clawpack.geoclaw.data

# ==============================================================================
#  Constants
# ==============================================================================
from data import Rearth, DEG2RAD, RAD2DEG, LAT2METER

# ======================
#  Remote file handling
# ======================
def strip_archive_extensions(path, extensions=["tar", "tgz", "bz2", "gz"]):
    r"""
    Strip off archive extensions defined in *extensions* list.

    Return stripped path calling this function recursively until all splitext
    does not provide an extension in the *extensions* list.

    """

    if os.path.splitext(path)[-1][1:] in extensions:
        return strip_archive_extensions(os.path.splitext(path)[0])
    else:
        return path


def get_remote_file(url, output_dir=None, file_name=None, force=False,  
                         verbose=False, ask_user=False, unpack=True):
    r"""Fetch file located at *url* and store at *output_dir*.

    :Input:
    
     - *url* (path) - URL to file to be downloaded.
     - *output_dir* (path) - Directory that the remote file will be downloaded
       to.  Defaults to the GeoClaw sratch directory defined by 
       *os.path.join(os.environ['CLAW'], 'geoclaw', 'scratch')*.
     - *file_name* (string) - Name of local file.  This defaults to the name of
       the remote file.
     - *force* (bool) - Force downloading of remote file regardless of whether
       it exists locally or not.  Default is *False*
     - *verbose* (bool) - Print out status information.  Default is *False*
     - *ask_user* (bool) - Whether to ask the user if it is ok to download the
       file before proceeding.  Default is *False*

    :Raises:
     
    Exceptions are raised from the *urllib2* module having to do with errors
    fetching the remote file.  Please see its documentation for more details of
    the exceptions that can be raised.

    returns *unarchived_output_path*
    """

    if output_dir is None:
        output_dir = os.path.join(os.environ['CLAW'], 'geoclaw', 'scratch')

    if file_name is None:
        file_name = os.path.basename(url)
        
    output_path = os.path.join(output_dir, file_name)
    unarchived_output_path = strip_archive_extensions(output_path)

    if not os.path.exists(unarchived_output_path) or force:

        if ask_user:
            ans = raw_input("  Ok to download topo file and save as %s?  \n"
                            % unarchived_output_path,
                            "     Type y[es], n[o].")
            if ans.lower() in ['y', 'yes']:
                if verbose:
                    print "*** Aborting download."
                return None
            
        if not os.path.exists(output_path):
            # Fetch remote file, will raise a variety of exceptions depending on
            # the retrieval problem if it happens
            if verbose:
                print "Downloading %s to %s..." % (url, output_path)
            with open(output_path, "w") as output_file:
                remote_file = urllib2.urlopen(url)
                output_file.write(remote_file.read())
            if verbose:
                print "Done downloading."
        elif verbose:
            print "File already exists, not downloading"

        if tarfile.is_tarfile(output_path) and unpack:
            if verbose:
                print "Un-archiving %s to %s..." % (output_path, 
                                                    unarchived_output_path)
            with tarfile.open(output_path, mode="r:*") as tar_file:
                tar_file.extractall(path=output_dir)
            if verbose:
                print "Done un-archiving."
    else:
        if verbose:
            print "Skipping %s " % url
            print "  because file already exists: %s" % output_path
        return None

    if unpack:
        return unarchived_output_path
    else:
        return output_path


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


def haversine(x, y, units='degrees'):
    r"""Compute the great circle distance on the earth between points x and y.


    """
    if units == 'degrees':
        # convert to radians:
        x *= DEG2RAD
        y *= DEG2RAD

    delta = [x[0] - y[0], x[1] - y[1]]

    # angle subtended by two points, using Haversine formula:
    dsigma = 2.0 * numpy.arcsin( numpy.sqrt( numpy.sin(0.5 * delta[1])**2   \
            + numpy.cos(x[1]) * numpy.cos(y[1]) * numpy.sin(0.5 * delta[0])**2))

    # alternative formula that may have more rounding error:
    #dsigma2 = arccos(sin(y1)*sin(y2)+ cos(y1)*cos(y2)*cos(dx))
    #print "max diff in dsigma: ", abs(dsigma-dsigma2).max()

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
