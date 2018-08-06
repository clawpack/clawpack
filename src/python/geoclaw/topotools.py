#!/usr/bin/env python
# encoding: utf-8

r"""
GeoClaw topotools Module  `$CLAW/geoclaw/src/python/geoclaw/topotools.py`

Module provides several functions for reading, writing and manipulating
topography (bathymetry) files.

:Classes:
 - Topography 

:Functions:

 - determine_topo_type
 - create_topo_func
 - topo1writer
 - topo2writer 
 - topo3writer 
 - swapheader


:TODO:
 - Add sub and super sampling capababilities
 - Add functions for creating topography based off a topo function, incorporate
   the create_topo_func into Topography class, maybe allow more broad 
   initialization ability to the class to handle this?
 - Fix `in_poly` function
 - Add remove/fill no data value
 - Add more robust plotting capabilities
"""

from __future__ import absolute_import
from __future__ import print_function
import os

import numpy

import clawpack.geoclaw.util as util
import clawpack.clawutil.data
import clawpack.geoclaw.data
import six
from six.moves import range

# ==============================================================================
#  Topography Related Functions
# ==============================================================================
def determine_topo_type(path, default=None):
    r"""Using the file suffix of path, attempt to deterimine the topo type.

    :Input:

     - *path* (string) - Path to the file.  Can include archive extensions (they
       will be stripped off). 
     - *default* (object) - Value returned if no suitable topo type was 
       determined.  Default is *None*.

    returns integer between 1-3 or *default* if nothing matches.
    
    """

    extension = os.path.splitext(
                  clawpack.clawutil.data.strip_archive_extensions(path))[-1][1:]
    
    topo_type = default
    if extension[:2] == "tt" or extension[:8] == 'topotype':
        topo_type = int(extension[-1])
    elif extension == 'xyz':
        topo_type = 1
    elif extension == 'asc':
        topo_type = 3
    elif extension == 'txyz':
        topo_type = 1
    elif extension == 'nc':
        topo_type = 4

    return topo_type


def create_topo_func(loc,verbose=False):
    """
    Given a 1-dimensional topography profile specfied by a set of (x,z) 
    values, create a lambda function that when evaluated will give the 
    topgraphy at the point (x,y).  (The resulting function is constant in y.)
    
    :Example: 

        >>> f = create_topo_func(loc)
        >>> b = f(x, y)
    
    :Input:
     - *loc* (list) - Create a topography file with the profile denoted by the
       tuples inside of loc.  A sample set of points are shown below.  Note 
       that the first value of the list is the x location and the second is 
       the height of the topography. ::


        z (m)
        ^                                                  o loc[5]  o
        |                                                    
        |                                          loc[4]   
        |--------------------------------------------o-----> x (m) (sea level)
        |                                            
        |                                o loc[2] o loc[3]
        |                         
        |                         
        |                           o loc[1]
        |           
        |                               
        |__________________o loc[0]
        0.0               
        
        
    """
    
    cmd_str = "lambda x,y: (x <= %s) * %s" % (loc[0][0],loc[0][1])
    for i in range(0,len(loc)-1):
        loc_str = " + (%s < x) * (x <= %s)" % (loc[i][0],loc[i+1][0])
        loc_str = "".join((loc_str," * ((%s - %s) " % (loc[i][1],loc[i+1][1])))
        loc_str = "".join((loc_str," / (%s - %s)" % (loc[i][0],loc[i+1][0])))
        loc_str = "".join((loc_str," * (x - %s) + %s)" % (loc[i][0],loc[i][1])))
        cmd_str = "".join((cmd_str,loc_str))
    cmd_str = "".join((cmd_str," + (%s < x) * %s" % (loc[-1][0],loc[-1][1])))
    
    if verbose:
        print(cmd_str)
    return eval(cmd_str)


def topo1writer (outfile,topo,xlower,xupper,ylower,yupper,nxpoints,nypoints):
    """
    Function topo1writer will write out the topofiles by evaluating the
    function topo on the grid specified by the other parameters.

    Assumes topo can be called on arrays X,Y produced by numpy.meshgrid.

    Output file is of "topotype1," which we use to refer to a file with
    (x,y,z) values on each line, progressing from upper left corner across
    rows, then down.
    """
    topography = Topography(topo_func=topo)

    topography.x = numpy.linspace(xlower,xupper,nxpoints)
    topography.y = numpy.linspace(ylower,yupper,nypoints)
    
    topography.write(outfile, topo_type=1)


def topo2writer (outfile,topo,xlower,xupper,ylower,yupper,nxpoints,nypoints, \
                 nodata_value=-99999):
    r"""Write out a topo type 2 file by evaluating the function *topo*.

    This routine is here for backwards compatibility and simply creates a new
    topography object and writes it out.

    """

    topography = Topography(topo_func=topo)

    topography.x = numpy.linspace(xlower,xupper,nxpoints)
    topography.y = numpy.linspace(ylower,yupper,nypoints)

    topography.write(outfile, topo_type=2)


def topo3writer (outfile,topo,xlower,xupper,ylower,yupper,nxpoints,nypoints, \
                 nodata_value=-99999):
    r"""Write out a topo type 3 file by evaluating the function *topo*.

    This routine is here for backwards compatibility and simply creates a new
    topography object and writes it out.

    """

    topography = Topography(topo_func=topo)

    topography.x = numpy.linspace(xlower,xupper,nxpoints)
    topography.y = numpy.linspace(ylower,yupper,nypoints)

    topography.write(outfile, topo_type=3)


def fetch_topo_url(url, local_fname=None, force=None, verbose=False, 
                        ask_user=False):
    """
    DEPRECATED:  Use *clawpack.clawutil.data.get_remote_file* instead (see note below).

    Replaces get_topo function.

    Download a topo file from the web, provided the file does not
    already exist locally.

    :Input:
        - *url* (str) URL including file name
        - *local_fname* (str) name of local file to create.  
          If *local_fname == None*, take file name from URL
        - *force* (bool) If False, prompt user before downloading.

    For GeoClaw examples, some topo files can be found in
    `http://www.geoclaw.org/topo`_
    See that website for a list of archived topo datasets.

    If force==False then prompt the user to make sure it's ok to download,

    If force==None then check for environment variable CLAW_TOPO_DOWNLOAD
    and if this exists use its value.  This is useful for the script
    python/run_examples.py that runs all examples so it won't stop to prompt.
    
    This routine has been deprecated in favor of 
    *clawpack.clawutil.data.get_remote_file*.  All the functionality should be 
    the same but calls the other routine internally.
    """

    if force is None:
        CTD = os.environ.get('CLAW_TOPO_DOWNLOAD', None)
        force = (CTD in [True, 'True'])

    if local_fname is not None:
        output_dir = os.path.dirname(local_fname)
        file_name = os.path.basename(local_fname)

    clawpack.clawutil.data.get_remote_file(url, output_dir=output_dir, 
                                                file_name=file_name,
                                                force=force, 
                                                verbose=verbose, 
                                                ask_user=ask_user)


def get_topo(topo_fname, remote_directory, force=None):
    """
    DEPRECATED:  Use *clawpack.geoclaw.util.get_remote_file* instead

    Download a topo file from the web, provided the file does not
    already exist locally.

    remote_directory should be a URL.  For GeoClaw data it may be a
    subdirectory of  http://www.clawpack.org/geoclaw/topo
    See that website for a list of archived topo datasets.

    If force==False then prompt the user to make sure it's ok to download,
    with option to first get small file of metadata.

    If force==None then check for environment variable CLAW_TOPO_DOWNLOAD
    and if this exists use its value.  This is useful for the script
    python/run_examples.py that runs all examples so it won't stop to prompt.
    """

    url = remote_directory + '/' + topo_fname
    clawpack.clawutil.data.get_remote_file(url, force=force)


def swapheader(inputfile, outputfile):
    r"""Swap the order of key and value in header to value first.

    Note that this is a wrapper around functionality in the Topography class.

    """
    topo = Topography(inputfile)
    topo.write(outputfile)



# ==============================================================================
#  Topography class
# ==============================================================================
class Topography(object):

    r"""Base topography class.

    A class representing a single topography file.

    :Properties:
    
    Note: Modified to check the `grid_registration` when reading or writing
    topo files and properly deal with `llcorner` registration in which case
    the x,y data should be offset by dx/2, dy/2 from the lower left corner 
    specified in the header of a DEM file.

    :Initialization:
         - 

    :Examples:

        >>> import clawpack.geoclaw.topotools as topo
        >>> topo_file = topo.Topography()
        >>> topo_file.read('./topo.tt3', topo_type=3)
        >>> topo_file.plot()

    """

    @property
    def z(self):
        r"""A representation of the data as an 1d array."""
        if (self._z is None) and self.unstructured:
            self.read(mask=False)
        return self._z
    @z.setter
    def z(self, value):
        self._z = value
    @z.deleter
    def z(self):
        del self._z

    @property
    def Z(self):
        r"""A representation of the data as a 2d array."""
        if self._Z is None:
            self.generate_2d_topo(mask=False)
        return self._Z
    @Z.setter
    def Z(self, value):
        self._Z = value
    @Z.deleter
    def Z(self):
        del self._Z

    @property
    def x(self):
        r"""One dimensional coorindate array in x direction."""
        if self._x is None:
            self.read(mask=False)
        return self._x
    @x.setter
    def x(self, value):
        self._extent = None
        self._x = value
    @x.deleter
    def x(self):
        del self._x

    @property
    def X(self):
        r"""Two dimensional coordinate array in x direction."""
        if self._X is None:
            self.generate_2d_coordinates(mask=False)
        return self._X
    @X.setter
    def X(self, value):
        self._extent = None
        self._X = value
        self._x = numpy.nan
    @X.deleter
    def X(self):
        del self._X

    @property
    def y(self):
        r"""One dimensional coordinate array in y direction."""
        if self._y is None:
            self.read(mask=False)
        return self._y
    @y.setter
    def y(self, value):
        self._extent = None
        self._y = value
    @y.deleter
    def y(self):
        del self._y

    @property
    def Y(self):
        r"""Two dimensional coordinate array in y direction."""
        if self._Y is None:
            self.generate_2d_coordinates(mask=False)
        return self._Y
    @Y.setter
    def Y(self, value):
        self._extent = None
        self._Y = value
        self._y = numpy.nan
    @Y.deleter
    def Y(self):
        del self._Y

    @property
    def extent(self):
        r"""Extent of the topography."""
        if self._extent is None:
            self._extent = ( numpy.min(self.x), numpy.max(self.x), 
                             numpy.min(self.y), numpy.max(self.y) )
        return self._extent
    @extent.setter
    def extent(self, value):
        self._extent = value

    @property
    def delta(self):
        r"""Spacing of data points."""
        if self._delta is None:
            if self.unstructured:

                # Calculate the smallest spacing between grid points            
                dx = numpy.infty
                dy = numpy.infty
                num_comparisons = self.x.shape[0] - 1
                for i in range(self.x.shape[0]):
                    for j in range(num_comparisons):
                        dx = min(dx, numpy.abs(self.x[i + j + 1] - self.x[i]))
                        dy = min(dy, numpy.abs(self.y[i + j + 1] - self.y[i]))

                    num_comparisons -= 1
                self._delta = [dx, dy]
            else:
                # All other topography types should have equally spaced grid
                # points in each direction
                begin_delta = numpy.array([abs(self.x[1] - self.x[0]), 
                                           abs(self.y[1] - self.y[0])])
                end_delta =   numpy.array([abs(self.x[-2] - self.x[-1]), 
                                           abs(self.y[-2] - self.y[-1])])
                if not numpy.allclose(begin_delta, end_delta, 1e-8):
                    raise ValueError("Grid spacing delta not constant, ",
                                     "%s != %s." % (begin_delta, end_delta))
                       
                dx = numpy.round(begin_delta[0], 15) 
                dy = numpy.round(begin_delta[1], 15) 
                self._delta = (dx, dy)
        return self._delta


    def __init__(self, path=None, topo_func=None, topo_type=None, 
                       unstructured=False):
        r"""Topography initialization routine.
        
        See :class:`Topography` for more info.

        """

        super(Topography, self).__init__()

        self.path = path
        self.topo_func = topo_func
        self.topo_type = topo_type

        self.unstructured = unstructured
        self.no_data_value = -9999

        # Data storage for only calculating array shapes when needed
        self._z = None
        self._Z = None
        self._x = None
        self._X = None
        self._y = None
        self._Y = None
        self._extent = None
        self._delta = None

        self.coordinate_transform = lambda x,y: (x,y)

        # RJL: should we read in by default if path is specified?
        #      If not, why include all these parameters in __init__?
        #if path:
        #    self.read(path=path, topo_type=topo_type, unstructured=unstructured,
        #     mask=mask, filter_region=filter_region)


    def generate_2d_topo(self, mask=False):
        r"""Generate a 2d array of the topo."""

        # Check to see if we need to generate these
        if self._Z is None:

            if self.unstructured:
                # Really no way to do this here with performing interpolation via
                # extract.  Note that if the interpolation is performed these
                # arrays are already stored in self._X and self._Y
                raise ValueError("Unstructured data does not allow for use of" \
                                 + " 2d arrays, first interpolate the data and" \
                                 + " try to perform this operation again.") 

            if self.path is not None:
                # RJL: why do we expect 1d z?
                if self._z is None:
                # Try to read the data, may not have done this yet
                    if self.topo_type is None:
                        self.topo_type = determine_topo_type(self.path)
                    if self.topo_type is None:
                        raise ValueError("topo_type must be specified")
                    self.read(path=self.path, topo_type=self.topo_type, mask=mask)
                    if self._Z is not None:
                        # We are done, the read function did our work
                        return

                # See if self._X and self._Y are already computed and use them if
                # available, otherwise just use self._x and self._y
                if self._X is not None and self._Y is not None:
                    new_shape = self._X.shape
                else:
                    new_shape = (self._x.shape[0], self._y.shape[0])
                # Reshape, note that the mask follows along with the new array
                self._Z = numpy.reshape(self._z, new_shape)

            elif self.topo_func is not None:
                # Generate topo via topo_func
                ## self._Z = numpy.flipud(self.topo_func(self.X, self.Y))
                ## RJL:  Don't flip -- leave so Z[i,j] has same dimensions as X,Y
                ## Othewise does not plot properly.
                self._Z = self.topo_func(self.X, self.Y)


    def generate_2d_coordinates(self, mask=False):
        r"""Generate 2d coordinate arrays."""

        # Check to see if we need to generate these
        if self._X is None and self._Y is None:

            # RJL: Added this to generate from _x and _y if available.
            # Correct?
            if (self._x is not None) and (self._y is not None):
                self._X,self._Y = numpy.meshgrid(self._x, self._y)

        if self._X is None and self._Y is None:
            if self.unstructured:
                # Really no way to do this here with performing interpolation via
                # extract.  Note that if the interpolation is performed these
                # arrays are already stored in self._X and self._Y
                raise ValueError("Unstructured data does not allow for use of" \
                                 + " 2d coordinates, first interpolate the data" \
                                 + " and try to perform this operation again.")

            if self.path is not None:
                if abs(self.topo_type) == 1:
                    # Reading this topo_type should produce the X and Y arrays
                    self.read(mask=mask)
                elif abs(self.topo_type) in [2,3]:
                    if self._x is None or self._y is None:
                        # Try to read the data to get these, may not have been done yet
                        self.read(mask=mask)
                    # Generate arrays
                    self._X, self._Y = numpy.meshgrid(self._x, self._y)
                else:
                    raise ValueError("Unrecognized topo_type: %s" % self.topo_type)

            elif self.topo_func is not None:
                if self._x is None or self._y is None:
                    raise ValueError("The x and y arrays must be set to ",
                                     "create 2d coordinate arrays.")
                self._X, self._Y = numpy.meshgrid(self._x, self._y)

            
            # If masking has been requested try to get the mask first from 
            # self._Z and then self._z
            if mask:
                if self._Z is None:
                    # Check to see if we really need to do anything here
                    if isinstance(self._z, numpy.ma.MaskedArray):
                        # Try to create self._Z
                        self.generate_2d_topo(mask=mask)

                if isinstance(self._Z, numpy.ma.MaskedArray):
                    # Use Z's mask for the X and Y coordinates
                    self._X = numpy.ma.MaskedArray(self._X, mask=self._Z.mask, 
                                                                     copy=False)
                    self._Y = numpy.ma.MaskedArray(self._Y, mask=self._Z.mask, 
                                                                     copy=False)


    def read(self, path=None, topo_type=None, unstructured=False, 
             mask=False, filter_region=None, force=False, stride=[1, 1],
             nc_params={}):
        r"""Read in the data from the object's *path* attribute.

        Stores the resulting data in one of the sets of *x*, *y*, and *z* or 
        *X*, *Y*, and *Z*.  

        :Input:
         - *path* (str)  file to read
         - *topo_type* (int) - GeoClaw format topo_type 
         - *unstructured* (bool) - default is False for lat-long grids.
         - *mask* (bool) - whether to store as masked array for missing
           values (default if False)
         - *filter_region* (tuple)
         - *stride* (list) - List of strides for the x and y dimensions
           respectively.  Default is *[1, 1]*.  Note that this is only
           implemented for NetCDF reading currently.
         - *nc_params* (dict) - 

        The first three might have already been set when instatiating object.

        """

        if (path is None) and (self.path is None):
            raise ValueError("*** Need to set path for file to read")

        if path:
            self.path = path   # set or perhaps reset
            self.topo_type = None  # force resetting below

        if unstructured:
            self.unstructured = unstructured

        # Check if the path is a URL and fetch data if needed or forced
        #if "http" in self.path:
        #    fetch_topo_url(self.path)
        # RJL: should switch to util.get_remote_file, but after fetching
        # still need to read it in, which that routine does not do.
        # Do we really want to support this?  Seems better for user
        # to fetch and store as desired filename and then read file.
            

        if self.topo_type is None:
            if topo_type is not None:
                self.topo_type = topo_type
            else:
                # Try to look at suffix for type
                self.topo_type = determine_topo_type(self.path)
                if self.topo_type is None:
                    #self.topo_type = 3
                    raise ValueError("topo_type must be specified")

        if self.unstructured:
            # Read in the data as series of tuples
            data = numpy.loadtxt(self.path)
            points = []
            values = []

            # Filter region if requested
            if filter_region is not None:
                for coordinate in data:
                    if filter_region[0] <= coordinate[0] <= filter_region[1]:
                        if filter_region[2] <= coordinate[1] <= filter_region[3]:
                            points.append(coordinate[0:2])
                            values.append(coordinate[2])

                if len(points) == 0:
                    raise Exception("No points were found inside requested " \
                                  + "filter region.")

                # Cast lists as ndarrays
                self._x = numpy.array(points[:,0])
                self._y = numpy.array(points[:,1])
                self._z = numpy.array(values)

            else:
                self._x = data[:,0]
                self._y = data[:,1]
                self._z = data[:,2]

        else:
            # Data is in one of the GeoClaw supported formats
            if abs(self.topo_type) == 1:
                data = numpy.loadtxt(self.path)
                N = [0,0]
                y0 = data[0,1]
                for (n, y) in enumerate(data[1:,1]):
                    if y != y0:
                        N[1] = n + 1
                        break
                N[0] = data.shape[0] // N[1]

                self._x = data[:N[1],0]
                self._y = data[::N[1],1]
                self._Z = numpy.flipud(data[:,2].reshape(N))
                dx = self.X[0,1] - self.X[0,0]
                dy = self.Y[1,0] - self.Y[0,0]
                self._delta = (dx,dy)

            elif abs(self.topo_type) in [2,3]:
                # Get header information
                N = self.read_header()  # note this also sets self._extent
                                        # self._x, self._y, self._delta, 
                                        # and  self.grid_registration

                if abs(self.topo_type) == 2:
                    # Data is read in as a single column, reshape it
                    self._Z = numpy.loadtxt(self.path, skiprows=6).reshape(N[1],N[0])
                    self._Z = numpy.flipud(self._Z)
                elif abs(self.topo_type) == 3:
                    # Data is read in starting at the top right corner
                    self._Z = numpy.flipud(numpy.loadtxt(self.path, skiprows=6))
        
                if mask:
                    self._Z = numpy.ma.masked_values(self._Z, self.no_data_value, copy=False)

            elif abs(self.topo_type) == 4:
                import netCDF4

                # NetCDF4 GEBCO topography
                with netCDF4.Dataset(self.path, 'r', format="NETCDF4") as nc_file:
                    x_var = nc_params.get('x_var', None)
                    y_var = nc_params.get('y_var', None)
                    z_var = nc_params.get('z_var', None)
                    for (key, var) in six.iteritems(nc_file.variables):
                        if 'axis' in var.ncattrs():
                            if var.axis.lower() == "x" and x_var is None:
                                x_var = key
                            elif var.axis.lower() == "y" and y_var is None:
                                y_var = key
                        else:
                            if z_var is None:
                                z_var = key

                    if x_var is None or y_var is None or z_var is None:
                        raise IOError("Could not automatically determine ",
                                      "variable ids.  Please check if the file",
                                      " has the 'axis' attribute attached to",
                                      " the appropriate x and y variables or ",
                                      "specify the variables directly via the",
                                      " *nc_params* dictionary.")


                    self._x = nc_file.variables[x_var][::stride[0]]
                    self._y = nc_file.variables[y_var][::stride[1]]
                    self._Z = nc_file.variables[z_var][::stride[0], 
                                                       ::stride[1]]

                if mask:
                    self._Z = numpy.ma.masked_values(self._Z, self.no_data_value, copy=False)
                    
            else:
                raise IOError("Unrecognized topo_type: %s" % self.topo_type)
                
            if self.topo_type < 0:
                # positive Z means distance below sea level for these
                # topo_type's, contrary to our convention, so negate:
                self._Z = -self._Z
                
            # Make sure these are set to None to force re-generating:
            self._X = None
            self._Y = None
                
            # Perform region filtering
            if filter_region is not None:
                # Find indices of region
                region_index = [None, None, None, None]
                region_index[0] = (self.x >= filter_region[0]).nonzero()[0][0]
                region_index[1] = (self.x <= filter_region[1]).nonzero()[0][-1]
                region_index[2] = (self.y >= filter_region[2]).nonzero()[0][0]
                region_index[3] = (self.y <= filter_region[3]).nonzero()[0][-1]

                self._x = self._x[region_index[0]:region_index[1]]
                self._y = self._y[region_index[2]:region_index[3]]

                # Force regeneration of 2d coordinate arrays and extent
                if self._X is not None or self._Y is not None:
                    del self._X, self._Y
                    self._X = None
                    self._Y = None
                self._extent = None

                # Modify Z array as well
                self._Z = self._Z[region_index[2]:region_index[3],
                                  region_index[0]:region_index[1]]


    def read_header(self):
        r"""Read in header of topography file at path.

        If a value returns numpy.nan then the value was not retrievable.  Note
        that this routine can read in headers whose values and labels are 
        swapped.

        """

        if abs(self.topo_type) in [2,3]:

            # Default values to track errors
            num_cells = [numpy.nan,numpy.nan]
            self._extent = [numpy.nan,numpy.nan,numpy.nan,numpy.nan]
            self._delta = numpy.nan

            with open(self.path, 'r') as topo_file:
                # Check to see if we need to flip the header values
                first_line = topo_file.readline()
                try:
                    num_cells[0] = int(first_line.split()[0])
                except ValueError:
                    # Assume the header is flipped from what we expect
                    num_cells[0] = int(first_line.split()[-1])
                    value_index = -1
                    label_index = 0
                else:
                    value_index = 0
                    label_index = -1

                num_cells[1] = int(topo_file.readline().split()[value_index])

                xline = topo_file.readline().split()
                xll = float(xline[value_index])
                # drop 'x' character and convert remaining string to lower case:
                x_registration = xline[label_index][1:].lower()

                yline = topo_file.readline().split()
                yll = float(yline[value_index])
                # drop 'y' character and convert remaining string to lower case:
                y_registration = yline[label_index][1:].lower()
                
                if x_registration == y_registration:
                    self.grid_registration = x_registration
                    # expect registration in ['llcorner', 'llcenter', 'lower']
                else:
                    raise IOError("x_registration and y_registration don't " \
                        + "match: %s,%s" % (x_registration, y_registration))

                # parse line allowing possibility of dx and dy (or just dx=dy)
                line = topo_file.readline()
                tokens = line.split() 
                values = []
                for token in tokens:
                    try:
                        v = float(token)
                        values.append(v)
                    except:
                        pass
                dx = values[0]
                if len(values) == 1:
                    dy = dx   # only dx given
                elif len(values) == 2:
                    dy = values[1]
                    self._delta = (values[0], values[1])  # if dx,dy on line
                else:
                    raise IOError("Cannot parse dx,dy line: %s" % line)
                self._delta = (dx, dy)
                    
                        
                self.no_data_value = float(topo_file.readline().split()[value_index])

                x = numpy.linspace(xll, xll+(num_cells[0]-1)*dx, num_cells[0])
                y = numpy.linspace(yll, yll+(num_cells[1]-1)*dy, num_cells[1])
                if self.grid_registration in ['lower', 'llcenter']:
                    # x,y are cell center / data locations:
                    self._x = x
                    self._y = y
                elif self.grid_registration == 'llcorner':
                    # x,y are lower left corners:
                    # data points are offset by dx/2, dy/2
                    self._x = x + dx/2.
                    self._y = y + dy/2.
                    print('*** Note: since grid registration is llcorner,')
                    print('    will shift x,y values by (dx/2, dy/2) to cell centers')
                else:
                    # assume that x,y are cell center / data locations:
                    self._x = x
                    self._y = y
                    print('*** Warning: Unrecognized grid_registration: %s' \
                                    % self.grid_registration)
                    print('    Assuming x,y at grid points')
                
                # set extent based on data locations (not lower corner for 'llcorner')
                self._extent = [self._x[0],self._x[-1],self._y[0],self._y[-1]]

        else:
            raise IOError("Cannot read header for topo_type %s" % self.topo_type)
            
        return num_cells

    def write(self, path, topo_type=None, no_data_value=None, masked=True, 
                header_style='geoclaw', Z_format="%15.7e", grid_registration=None):
        r"""Write out a topography file to path of type *topo_type*.

        Writes out a topography file of topo type specified with *topo_type* or
        inferred from the output file's extension, defaulting to 3, to path
        from data in Z.  The rest of the arguments are used to write the header
        data.

        :Input:
         - *path* (str)  - file to write
         - *topo_type* (int) - GeoClaw format topo_type 
           **Note:** this is second positional argument, agreeing with
           the read function in this class.  It was the third argument in
           GeoClaw version 5.3.1 and earlier.  
         - *no_data_value* - values used to indicate missing data
         - *masked* (bool) - unused??
         - *header_style* (str) - indicates format of header lines
             'geoclaw' or 'default'  ==> write value then label 
                                     with grid_registration == 'lower' as default
             'arcgis' or 'asc' ==> write label then value  
                                   with grid_registration == 'llcorner' as default
                                   (needed for .asc files in ArcGIS)

         - *Z_format* (str) - string format to use for Z values
           The default format "%15.7e" gives at least millimeter precision
           for topography with abs(Z) < 10000 and results in
           smaller files than the previous default of "%22.15e" used in
           GeoClaw version 5.3.1 and earlier.  A shorter format can be used
           if the user knows there are fewer significant digits, e.g.
           etopo1 data is integers and so has a resolution of 1 meter.
           In this case a cropped or coarsened version might be written
           with `Z_format = "%7i"`, for example.
         - *grid_registration* (str) - 'lower', 'llcorner', 'llcenter' 
                or None for defaults described above.

        """

        # Determine topo type if not specified
        if topo_type is None:
            # Look at the the suffix of the path and the object's topo_type
            # attribute to try to deterimine which to use, default to the path
            # version unless it did not work
            path_topo_type = determine_topo_type(path, default=-1)
            
            if self.topo_type is not None and path_topo_type == -1:
                topo_type = self.topo_type
            elif path_topo_type != -1:
                topo_type = path_topo_type
            else:
                # Default to 3 if all else fails
                topo_type = 3

        # Default to this object's no_data_value if the passed is None, 
        # otherwise the argument will override the object's value or it will 
        # default to -9999 (default for the class)
        if no_data_value is None:
            no_data_value = self.no_data_value

        # Check to see if masks have been applied to topography, if so use them
        # if masked is True
        if isinstance(self.Z, numpy.ma.MaskedArray) and masked:
            pass
        else:
            pass

        if self.unstructured:
            with open(path, 'w') as outfile:
                for (i, topo) in enumerate(self.z):
                    outfile.write("%s %s %s\n" % (self.x[i], self.y[i], topo))

        elif topo_type == 1:
            # longitudes = numpy.linspace(lower[0], lower[0] + delta * Z.shape[0], Z.shape[0])
            # latitudes = numpy.linspace(lower[1], lower[1] + delta * Z.shape[1], Z.shape[1])

            with open(path, 'w') as outfile:
                for j in range(len(self.y)-1, -1, -1):
                    latitude = self.y[j]
                    for (i, longitude) in enumerate(self.x):
                        outfile.write("%s %s %s\n" % (longitude, latitude, self.Z[j,i]))

        elif topo_type == 2 or topo_type == 3:

            if grid_registration is None:
                if header_style in ['geoclaw','default']:
                    grid_registration = 'lower'
                elif header_style in ['arcgis','asc']:
                    grid_registration = 'llcorner'
                else:
                    raise ValueError("*** Unrecognized header_style")

            if grid_registration in ['lower','llcenter']:
                xlower = self.x[0]
                ylower = self.y[0]
            elif grid_registration == 'llcorner':
                xlower = self.x[0] - self.delta[0]/2.
                ylower = self.y[0] - self.delta[1]/2.
            else:
                raise ValueError('Unrecognized grid_registration: %s' \
                                % grid_registration)
            xlabel = 'x' + grid_registration
            ylabel = 'y' + grid_registration

            with open(path, 'w') as outfile:
                # Write out header
                if header_style in ['geoclaw','default']:
                    outfile.write('%6i                              ncols\n' % self.Z.shape[1])
                    outfile.write('%6i                              nrows\n' % self.Z.shape[0])
                    outfile.write('%22.15e              %s\n' % (xlower,xlabel))
                    outfile.write('%22.15e              %s\n' % (ylower,ylabel))
                    if abs(self.delta[0] - self.delta[1])/self.delta[0] < 1e-8:
                        # write only dx in usual case:
                        outfile.write('%22.15e              cellsize\n' \
                                % self.delta[0])
                    else:
                        # write both dx and dy if they differ:
                        outfile.write('%22.15e    %22.15e          cellsize\n' \
                                % (self.delta[0], self.delta[1]))
                    outfile.write('%10i                          nodata_value\n' % no_data_value)
                elif header_style in ['arcgis','asc']:
                    outfile.write('ncols  %6i\n' % self.Z.shape[1])
                    outfile.write('nrows  %6i\n' % self.Z.shape[0]) 
                    outfile.write('%s  %22.15e\n' % (xlabel,xlower))
                    outfile.write('%s  %22.15e\n' % (ylabel,ylower))
                    outfile.write('cellsize %22.15e\n'  % self.delta[0])
                    outfile.write('nodata_value  %10i\n' % no_data_value)
                else:
                    raise ValueError("*** Unrecognized header_style")

                masked_Z = isinstance(self.Z, numpy.ma.MaskedArray)

                # Write out topography data
                if topo_type == 2:
                    Z_format = Z_format + "\n"
                    if masked_Z:
                        Z_filled = numpy.flipud(self.Z.filled())
                    else:
                        Z_filled = numpy.flipud(self.Z)
                    for i in range(self.Z.shape[0]):
                        for j in range(self.Z.shape[1]):
                            outfile.write(Z_format % Z_filled[i,j])
                    if masked_Z:
                        del Z_filled
                elif topo_type == 3:
                    Z_format = Z_format + " "
                    if masked_Z:
                        Z_flipped = numpy.flipud(self.Z.filled())
                    else:
                        Z_flipped = numpy.flipud(self.Z)
                    for i in range(self.Z.shape[0]):
                        for j in range(self.Z.shape[1]):
                            outfile.write(Z_format % (Z_flipped[i,j]))
                        outfile.write("\n")
                    if masked_Z:
                        del Z_flipped
        elif topo_type == 4:
            # Write out netCDF4 topography
            import netCDF4
            
            with netCDF4.Dataset(path, 'w') as outfile:
                # Add root attributes
                outfile.Conventions = "CF-1.6"
                outfile.title = "Topography Data"
                outfile.institution = "Unknown"
                outfile.source = "Unknown"
                outfile.history = "" # TODO: Add current date here
                outfile.references = ""
                outfile.comment = "Created by GeoClaw"

                outfile.createDimension("lon", self.x.shape[0])
                outfile.createDimension("lat", self.y.shape[0])
                
                lon = outfile.createVariable('lon', 'f8', ('lon',))
                lon.standard_name = "longitude"
                lon.long_name = "longitude"
                lon.units = "degrees_east"
                lon.axis = "X"
                lon[:] = self.x

                lat = outfile.createVariable('lat', 'f8', ('lat',))
                lat.standard_name = "latitude"
                lat.long_name = "latitude"
                lat.units = "degrees_north"
                lat.axis = "Y"
                lat[:] = self.y

                elevation = outfile.createVariable('elevation', 'f8',
                    ('lat','lon',))
                elevation.standard_name  = "height_above_reference_ellipsoid"
                elevation.long_name  = "Elevation relative to sea level"
                elevation.units  = "m"
                elevation.scale_factor  = 1.0
                elevation.add_offset  = 0.0
                elevation.sdn_parameter_urn  = "SDN:P01::BATHHGHT"
                elevation.sdn_parameter_name  = "Sea floor height (above mean sea level) {bathymetric height}"
                elevation.sdn_uom_urn  = "SDN:P06:ULAA"
                elevation.sdn_uom_name  = "Metres"
                elevation.positive = "up"
                elevation[:, :] = self.Z

                elevation.no_data_value = self.no_data_value


        else:
            raise NotImplemented("Output type %s not implemented." % topo_type)


    def plot(self, axes=None, contour_levels=None, contour_kwargs={}, 
             limits=None, cmap=None, add_colorbar=True, 
             plot_box=False, fig_kwargs={}, data_break=0.):
        r"""Plot the topography.

        :Input:
         - *axes* (matplotlib.pyplot.axes) - If passed in, plot will be
           added to this axes.  Otherwise a new plot figure will be created
           (using *fig_kwargs*) and a new *axes* object created and returned.
         - *contour_levels* (list) - levels for contour lines if these are 
           to be added (default *None*).  Set to [0.] to plot shoreline.
         - *contour_kwargs* (dict) - keyword arguments to be passed to
           contour command, e.g. {'colors':'r', 'linestyles': '-'}.
           Default is empty dict.
         - *limits* (list) - (min, max) of topo values for color map.  
           Defaults to None, in which case (self.Z.min(), self.Z.max()) used.
         - *cmap* (matplotlib.colors.Colormap) - colormap, defaults to
           specialized map with blues for bathymetry and green/browns for topo.
         - *fig_kwargs* (dict) - keyword arguments to be passed to figure.
         - *plot_box* (bool or color specifier) - If evaluates to True, plot
           a box around limits of this topo. 
         - *data_break* (float) - when default cmap is used, the value to use
           to break between water and land colormaps.
           Defaults to 0., but for some topo files may need to use e.g. 0.01
           Or may want to show plots at different tide stage.

        :Output:
         - *axes* (matplotlib.pyplot.axes) - the axes on which plot created.

        Note that:
          - if *type(self.Z)* is *numpy.ma.MaskedArray* then *pcolor* is used,
          - if *type(self.Z)* is *numpy.ndarray* then *imshow* is used.
            (This is faster for large files)
        """

        import matplotlib.pyplot as plt
        import matplotlib.colors as colors

        import clawpack.visclaw.colormaps as colormaps

        # Create axes if needed
        if axes is None:
            fig = plt.figure(**fig_kwargs)
            axes = fig.add_subplot(111)
        
        # Turn off annoying offset
        axes.ticklabel_format(format="plain", useOffset=False)
        for label in axes.get_xticklabels():
            label.set_rotation(20)

        region_extent = self.extent


        if limits is None:
            if self.unstructured:
                topo_extent = (numpy.min(self.z), numpy.max(self.z))
            else:
                topo_extent = (numpy.min(self.Z), numpy.max(self.Z))
        else:
            topo_extent = limits

        # Create color map - assume shore is at z = 0.0
        if cmap is None:
            land_cmap = colormaps.make_colormap({ 0.0:[0.1,0.4,0.0],
                                                 0.25:[0.0,1.0,0.0],
                                                  0.5:[0.8,1.0,0.5],
                                                  1.0:[0.8,0.5,0.2]})
            sea_cmap = plt.get_cmap('Blues_r')
            if topo_extent[0] >= 0.0:
                cmap = land_cmap
                norm = colors.Normalize(vmin=0.0, vmax=topo_extent[1])
            elif topo_extent[1] <= 0.0:
                cmap = sea_cmap
                norm = colors.Normalize(vmin=topo_extent[0], vmax=0.0)
            else:
                cmap, norm = colormaps.add_colormaps((land_cmap, sea_cmap), 
                                                     data_limits=topo_extent,
                                                     data_break=data_break)
        else:
            norm = colors.Normalize(vmin=topo_extent[0], vmax=topo_extent[1])

        # Plot data
        if self.unstructured:
            plot = axes.scatter(self.x, self.y, c=self.z, cmap=cmap, norm=norm,
                                                marker=',', linewidths=(0.0,))
        elif isinstance(self.Z, numpy.ma.MaskedArray):
            # Adjust coordinates so color pixels centered at X,Y locations
            plot = axes.pcolor(self.X - self.delta[0] / 2.0, 
                               self.Y - self.delta[1] / 2.0, 
                               self.Z,
                               norm=norm,
                               cmap=cmap)
        else:
            plot = axes.imshow(self.Z, norm=norm, extent=region_extent, 
                                       cmap=cmap, origin='lower',
                                       interpolation='nearest')
        if add_colorbar:
            cbar = plt.colorbar(plot, ax=axes)
            cbar.set_label("Topography (m)")
        # levels = range(0,int(-numpy.min(Z)),500)

        if (contour_levels is not None) and (not self.unstructured):
            axes.contour(self.X, self.Y, self.Z, levels=contour_levels,
                 **contour_kwargs)

        axes.set_xlim(region_extent[0:2])
        axes.set_ylim(region_extent[2:])

        x1,x2,y1,y2 = self.extent
        if plot_box:
            # plot a box around this topography region
            if type(plot_box) is bool:
                color = 'm'
            else:
                # assume plot_box is a valid color:
                color = plot_box
            plt.plot([x1,x2,x2,x1,x1], [y1,y1,y2,y2,y1], color=color)

        mean_lat = 0.5 * (region_extent[3] + region_extent[2])
        axes.set_aspect(1.0 / numpy.cos(numpy.pi / 180.0 * mean_lat))

        return axes


    def interp_unstructured(self, fill_topo, extent=None, method='nearest',
                                   delta=None, delta_limit=20.0, 
                                   no_data_value=-9999, buffer_length=100.0,
                                   proximity_radius=100.0, 
                                   resolution_limit=2000):
        r"""Interpolate unstructured data on to regular grid.

        Function to interpolate the unstructured data in the topo object onto a
        structured grid.  Utilizes a bounding box plus a buffer of size 
        *buffer_length* (meters) containing all data unless *extent is not None*
        is *True*.  Then uses the fill topography *fill_topo* to fill in the
        gaps in the unstructured data.  By default this is done by masking the 
        fill data with the extents, the value *no_data_value* and if 
        *proximity_radius* (meters) is not 0, by a radius of *proximity_radius* 
        from all grid points in the object.  Stores the 
        result in the *self.X*, *self.Y* and *self.Z* object attributes.  The
        resolution of the final grid is determined by calculating the minimum
        distance between all *x* and *y* data with a hard lower limit of 
        *delta_limit* (meters).

        Note that the function *scipy.interpolate.griddata* does not respect
        masks so a call to *numpy.ma.MaskedArray.compressed()* must be made to 
        remove the masked data.

        :Input:
         - *fill_topo* (list) - List of Topography objects to use as fill data
           in the projection.
         - *extent* (tuple) - A tuple defining the rectangle of the sub-section.  
           Must be in the form (x lower,x upper,y lower, y upper).
         - *method* (string) - Method used for interpolation, valid methods are
           found in *scipy.interpolate.griddata*.  Default is *nearest*.
         - *delta* (tuple) - Directly set the grid spacing of the interpolation
           rather than determining it from the data itself.  Defaults to *None*
           which causes the method to determine this value itself.
           Should be a 2-tuple of floats (delta_x, delta_y).
         - *delta_limit* (float) - Limit of finest horizontal resolution, 
           default is 20 meters.
         - *no_data_value* (float) - Value to use if no data was found to fill in a 
           missing value, ignored if `method = 'nearest'`. Default is `-9999`.
         - *buffer_length* (float) - Buffer around bounding box, only applicable
           when *extent* is None.  Default is `100.0` meters.
         - *proximity_radius* (float) - Radius every unstructured data point
           used to mask the fill data with.  Default is `100.0` meters.
         - *resolution_limit* (int) - Limit the number of grid points in a
           single dimension.  Raises a *ValueError* if the limit is violated.
           Default value is `2000`.

        Sets this object's *unstructured* attribute to *False* if successful.

        """

        import scipy.interpolate as interpolate

        # Convert meter inputs to degrees
        mean_latitude = numpy.mean(self.y)
        buffer_degrees = util.dist_meters2latlong(buffer_length, 0.0, 
                                                  mean_latitude)[0]
        delta_degrees = util.dist_meters2latlong(delta_limit, 0.0, 
                                                 mean_latitude)[0]
        if proximity_radius > 0.0:
            proximity_radius_deg = util.dist_meters2latlong(proximity_radius, 
                                                            0.0,
                                                            mean_latitude)[0]
            
        # Calculate new grid coordinates
        if extent is None:
            extent = [ numpy.min(self.x) - buffer_degrees, 
                       numpy.max(self.x) + buffer_degrees, 
                       numpy.min(self.y) - buffer_degrees, 
                       numpy.max(self.y) + buffer_degrees ]
        if delta is None:
            delta_x = max( numpy.abs(self.x[1:] - self.x[:-1]).min(), delta_degrees)
            delta_y = max( numpy.abs(self.y[1:] - self.y[:-1]).min(), delta_degrees)
        else:   
            try:
                delta_x, delta_y = delta   # tuple provided
            except:
                delta_x = delta_y = delta  # assume float provided
                
        N = ( numpy.ceil((extent[1] - extent[0]) / delta_x),
              numpy.ceil((extent[3] - extent[2]) / delta_y) )
        if not numpy.all(N[:] < numpy.ones((2)) * resolution_limit):
            ValueError("Calculated resolution too high, N=%s!" % str(N))
        self._X, self._Y = numpy.meshgrid( 
                                     numpy.linspace(extent[0], extent[1], N[0]),
                                     numpy.linspace(extent[2], extent[3], N[1]))

        # Add the unstructured points to the data
        points = numpy.array([self.x, self.y]).transpose()
        values = self.z

        # Mask fill topography and flatten the arrays if needed
        if not isinstance(fill_topo, list):
            fill_topo = [fill_topo]
        for topo in fill_topo:
            if topo.unstructured:
                x_fill = topo.x
                y_fill = topo.y
                z_fill = topo.z

                extent_mask = extent[0] > x_fill
                extent_mask = numpy.logical_or(extent_mask,extent[1] < x_fill)
                extent_mask = numpy.logical_or(extent_mask,extent[2] > y_fill)
                extent_mask = numpy.logical_or(extent_mask,extent[3] < y_fill)
                
                # Create fill no-data value mask
                no_data_mask = numpy.logical_or(extent_mask, z_fill == no_data_value)

                all_mask = numpy.logical_or(extent_mask, no_data_mask)

                # Create proximity mask
                if proximity_radius > 0.0:
                    indices = (~all_mask).nonzero()
                    for n in range(indices[0].shape[0]):
                        i = indices[0][n]
                        all_mask[i] = numpy.any(numpy.sqrt((self.x - x_fill[i])**2 
                                                         + (self.y - y_fill[i])**2)
                                                     < proximity_radius_deg)

                x_fill_masked = numpy.ma.masked_where(all_mask, x_fill)
                y_fill_masked = numpy.ma.masked_where(all_mask, y_fill)
                z_fill_masked = numpy.ma.masked_where(all_mask, z_fill)    

                # Add the fill bathymetry to points and values
                fill_points = numpy.column_stack((x_fill_masked.compressed(), 
                                                  y_fill_masked.compressed()))

                points = numpy.concatenate((fill_points, points))
                values = numpy.concatenate((z_fill_masked.compressed(), values))

            else:
                # Structured fill data
                X_fill = topo.X
                Y_fill = topo.Y
                Z_fill = topo.Z

                # Create extent mask
                extent_mask = extent[0] > X_fill
                extent_mask = numpy.logical_or(extent_mask,extent[1] < X_fill)
                extent_mask = numpy.logical_or(extent_mask,extent[2] > Y_fill)
                extent_mask = numpy.logical_or(extent_mask,extent[3] < Y_fill)
                
                # Create fill no-data value mask
                no_data_mask = numpy.logical_or(extent_mask, Z_fill == no_data_value)

                all_mask = numpy.logical_or(extent_mask, no_data_mask)

                # Create proximity mask
                if proximity_radius > 0.0:
                
                    indices = (~all_mask).nonzero()
                    for n in range(indices[0].shape[0]):
                        i = indices[0][n]
                        j = indices[1][n]
                        all_mask[i,j] = numpy.any(numpy.sqrt((self.x - X_fill[i,j])**2 
                                                           + (self.y - Y_fill[i,j])**2)
                                                     < proximity_radius_deg)

                X_fill_masked = numpy.ma.masked_where(all_mask, X_fill)
                Y_fill_masked = numpy.ma.masked_where(all_mask, Y_fill)
                Z_fill_masked = numpy.ma.masked_where(all_mask, Z_fill)    

                # Add the fill bathymetry to points and values
                fill_points = numpy.column_stack((X_fill_masked.compressed(),
                                                  Y_fill_masked.compressed()))

                points = numpy.concatenate((fill_points, points))
                values = numpy.concatenate((Z_fill_masked.compressed(), values))

        # Use specified interpolation
        self._Z = interpolate.griddata(points, values, (self.X, self.Y), 
                                                                  method=method)

        self._extent = extent
        self._delta = (delta_x, delta_y)
        self.unstructured = False


    def in_poly(self, polygon):
        r"""Mask points (x,y) that are not in the specified polygon.

        Uses simple ray casting algorithm for speed so beware of corner cases!

        :Input:
        
         - *polygon* (list) List of points that comprise the polygon.  Note that
           order of the points will effect if this works (positive versus negative
           winding order).  Points should be in counter-clockwise arrangement.

        :Returns:
        
         - *X_mask* (numpy.ma.MaskedArray) Masked array of X coordinates where those
           points outside of the polygon have been masked.
         - *Y* (numpy.ndarray) Coordinates in y direction in a meshgrid type of
           configuration.

        """
        raise NotImplemented("This function is not quite working yet, please "+\
                             "try again later")

        TOLERANCE = 1e-6

        # Flatten the input arrays to make this a bit easier
        x = self.X.flatten()
        y = self.Y.flatten()

        # Construct edges
        edges = []
        for edge in range(len(polygon) - 1):
            edges.append([polygon[edge], polygon[edge+1]])
        edges.append([polygon[-1], polygon[0]])

        # Check for intersections
        num_intersections = numpy.zeros(x.shape[0])

        for edge in edges:
            # Check for a vertical line
            if numpy.abs(edge[0][0] - edge[1][0]) < TOLERANCE:
                x_intersect = edge[0][0]        
            else:
                edge_slope = (edge[0][1] - edge[1][1]) / (edge[0][0] - edge[1][0])
                x_intersect = (y - edge[0][1]) / edge_slope + edge[0][0]

            num_intersections += (min(edge[0][1], edge[1][1]) <= y) * \
                                 (max(edge[0][1], edge[1][1]) >= y) * \
                                 (x_intersect <= x)
                                 

        # General intersection of two lines
        intersect = (numpy.mod(num_intersections, numpy.ones(x.shape) * 2) != 1)

        # Return masked arrays that are reshaped back to the input shapes
        return numpy.ma.masked_where(intersect, x, copy=False).reshape(self.X.shape), \
               numpy.ma.masked_where(intersect, y, copy=False).reshape(self.Y.shape)


    def replace_values(self, indices, value=numpy.nan, method='fill'):
        r"""Replace the values at *indices* by the specified method

        :Methods:
         - "fill"
         - "nearest"
        """

        # Average surrounding good data
        if method == 'fill':
            for index in indices:
                r = 0
                point_replaced = False
                while not point_replaced and r < max(self.Z.shape):
                    r = r + 1
                    i_range = list(range(max(0, index[0] - r), min(index[0] + r + 1, self.Z.shape[0])))
                    j_range = list(range(max(0, index[1] - r), min(index[1] + r + 1, self.Z.shape[1])))
                    num_points = 0
                    summation = 0.0
                    for i in i_range:
                        for j in j_range:
                            if (i,j) not in indices:
                                summation += self.Z[i,j]
                                num_points += 1
                    if num_points > 0:
                        self.Z[index[0], index[1]] = summation / num_points
                        point_replaced = True

        elif method == "nearest":
            pass


    def replace_no_data_values(self, method='fill'):
        r"""Replace *no_data_value* with other values as specified by *method*.

        self.no_data_value

        :Input:
         - *method* can be one of:

             - *fill* - Fill in all *no_data_value* locations with *value*
             - *nearest* - Fill in *no_data_value* locations with 
               average of nearest neighbors.

        """
        raise NotImplemented("This functionality has not been added yet.")
        no_data_value_indices = (self.Z == self.no_data_value).nonzero()
        self.replace_values(no_data_value_indices, method=method)

        # nrows= shape(Z)[0]
        # ncols= shape(Z)[1]
        # npts = nrows*ncols

        # xi=X[0,:]
        # yi=Y[:,0]

        # X.np.reshape(npts)
        # Y.np.reshape(npts)
        # Z.np.reshape(npts)

        # ind=np.where(Z!=nodata_value)
        # X=X[ind]
        # Y=Y[ind]
        # Z=Z[ind]

        # ptsremove=npts-len(Z)
        # if ptsremove>0:
        #     print("Removing %s nodata_value points" % ptsremove)

        # Z = pylab.griddata(X,Y,Z,xi,yi)
        # (X,Y)=np.meshgrid(xi,yi)

        # griddata2topofile(X,Y,Z,outputfile,topotypeout,nodata_value,nodata_value)


    def smooth_data(self, indices, r=1):
        r"""Filter topo data at *indices* by averaging surrounding data.

        Surrounding data is considered within the ball of radius *r* in the 
        inf-norm.  Acts as a low-band pass filter and removes oscillatory data.

        :Input:
         - *indices* (list)
         - *r* (int) 

        :Output:
         None
        """

        index_range = [None, None]
        for index in indices:
            for n in range(2):
                index_range[n] = list(range(max(0, index[n] - r), 
                                       min(index[n] + r + 1, self.Z.shape[n])))
            num_points = 0
            summation = 0.0
            for i in index_range[0]:
                for j in index_range[1]:
                    summation += self.Z[i,j]
                    num_points += 1
            if num_points > 0:
                self.Z[index[0], index[1]] = summation / num_points


    def crop(self, filter_region=None, coarsen=1):
        r"""Crop region to *filter_region*

        Create a new Topography object that is identical to this one but cropped
        to the region specified by filter_region

        :TODO:
         - Currently this does not work for unstructured data, could in principle
         - This could be a special case of in_poly although that routine could
           leave the resulting topography as unstructured effectively.
        """
        
        if self.unstructured:
            raise NotImplemented("*** Cannot currently crop unstructured topo")

        if filter_region is None:
            # only want to coarsen, so this is entire region:
            filter_region = [self.x[0],self.x[-1],self.y[0],self.y[-1]]

        # Find indices of region
        region_index = [None, None, None, None]
        region_index[0] = (self.x >= filter_region[0]).nonzero()[0][0]
        region_index[1] = (self.x <= filter_region[1]).nonzero()[0][-1] + 1
        region_index[2] = (self.y >= filter_region[2]).nonzero()[0][0]
        region_index[3] = (self.y <= filter_region[3]).nonzero()[0][-1] + 1
        newtopo = Topography()

        newtopo._x = self._x[region_index[0]:region_index[1]:coarsen]
        newtopo._y = self._y[region_index[2]:region_index[3]:coarsen]

        # Force regeneration of 2d coordinate arrays and extent if needed
        newtopo._X = None
        newtopo._Y = None
        newtopo._extent = None

        # Modify Z array as well
        newtopo._Z = self._Z[region_index[2]:region_index[3]:coarsen,
                          region_index[0]:region_index[1]:coarsen]

        newtopo.unstructured = self.unstructured
        newtopo.topo_type = self.topo_type

        # print "Cropped to %s by %s array"  % (len(newtopo.x),len(newtopo.y))
        return newtopo

    def make_shoreline_xy(self, sea_level=0):
        r"""
        Returns an array *shoreline_xy* with 2 columns containing x and y values
        for all segements of the shoreline (defined to be the contour 
        where self.z = sea_level) separated by [nan,nan] pairs.  
        This allows all shorelines to be quickly plotted via:

            >>> plot(shoreline_xy[:,0], shoreline_xy[:,1])

        The shoreline can be saved as a binary *.npy* file via:

            >>> numpy.save(filename, shoreline_xy)

        which is much smaller than the original topography file. 
        Reload via:

            >>> shoreline_xy = numpy.load(filename)
        """

        import matplotlib.pyplot as plt

        x = self.x
        y = self.y
        Z = self.Z
        c = plt.contour(x,y,Z,[sea_level])  
        # c is the level 0 contour as list of arrays, one for each segement
        # catenate these together separated by array([nan,nan]):
        shoreline_xy = c.allsegs[0][0]  # first segment
        for k in range(1,len(c.allsegs[0])):
            shoreline_xy = numpy.vstack((shoreline_xy, \
                           numpy.array([numpy.nan,numpy.nan]), c.allsegs[0][k]))
        return shoreline_xy



# Define convenience dictionary of URLs for some online DEMs in netCDF form:
remote_topo_urls = {}

# global 1 arcminute topography:
remote_topo_urls['etopo1'] = \
    'https://www.ngdc.noaa.gov/thredds/dodsC/global/ETOPO1_Ice_g_gmt4.nc'

# some 1/3 arcsecond coastal modeling DEMs:
server = 'https://www.ngdc.noaa.gov/thredds/dodsC/regional/'
remote_topo_urls['astoria'] = server + 'astoria_13_mhw_2012.nc'
remote_topo_urls['puget_sound'] = server + 'puget_sound_13_mhw_2014.nc'
remote_topo_urls['port_townsend'] = server + 'port_townsend_13_mhw_2011.nc'
remote_topo_urls['strait_of_juan_de_fuca'] = \
    server + 'strait_of_juan_de_fuca_13_navd88_2015.nc'



def read_netcdf(path, zvar=None, extent='all', coarsen=1, return_topo=True, 
                return_xarray=False, verbose=False):

    """
    :Input:

     - *path* (str) - Path to the file to read, or url to remote file,
       or a key into the topotools.remote_topo_urls dictionary.
     - *zvar* (str) - variable to read as Z=elevation. 
       if None, will try 'Band1', 'z', 'elevation'.
     - *extent* - [x1,x2,y1,y2] for desired subset, or 'all' for entire file
     - *coarsen* (int) - factor to coarsen by, 1 by default.
     - *return_topo* (bool) - if True, return a topotools.Topography object.
       default is True
     - *return_xarray* (bool) - if True, return an xarray.Dataset object.
       default is False

    :Output:
     - topo and/or xarray_ds depending on what was requested.
       (either a single object or a tuple of two objects.)

    If `return_xarray == True` then `xarray` is used to read the data,
    otherwise `netCDF4` is used directly.

    Sample usage:

        from clawpack.geoclaw import topotools
        extent = [-126,-122,46,49]
        path = 'etopo1'
        topo = topotools.read_netcdf(path, extent=extent, coarsen=2, \
                                     verbose=True)

        # to plot:
        topo.plot()

        # to save topofile for input to GeoClaw:
        topo.write('etopo_sample_2min.tt3', topo_type=3, Z_format='%.0f')  

    This should give a 2-minute resolution DEM of the Western Washington coast.
    Note that etopo1 Z values are integers (vertical resolution is 1 meter)
    and using `Z_format='%.0f'` will save as integers to minimize file size.
    """
    
    from clawpack.geoclaw import topotools
    from numpy import array
    import netCDF4
    if return_xarray:
        import xarray

    # check if path is a key in the remote_topo_urls dictionary:
    if path in remote_topo_urls.keys():
        path = remote_topo_urls[path]

    if verbose:
        print("Will read netCDF data from \n    %s" % path)

    if return_xarray:
        f = xarray.open_dataset(path)
    else: 
        f = netCDF4.Dataset(path, 'r')

    x = f.variables['lon']
    y = f.variables['lat']

    # for selecting subset based on extent, convert to arrays if netCDF4 used:
    if not return_xarray:
        x = array(x)
        y = array(y)
    
    if zvar is None:
        if 'Band1' in f.variables:
            zvar = 'Band1'
        elif 'z' in f.variables:
            zvar = 'z'
        elif 'elevation' in f.variables:
            zvar = 'elevation'
        else:
            raise ValueError("*** Unrecognized zvar in netCDF file")

    if extent == 'all':
        xs = x[::coarsen]
        ys = y[::coarsen]
        Zs = f.variables[zvar][::coarsen,::coarsen]
    else:
        x1,x2,y1,y2 = extent
        # find indices of x,y arrays for points lying within extent:
        iindex = [i for i,xi in enumerate(x) if (x1 <= xi <= x2)]
        jindex = [j for j,yj in enumerate(y) if (y1 <= yj <= y2)]
        i1 = iindex[0]
        i2 = iindex[-1] + 1
        j1 = jindex[0]
        j2 = jindex[-1] + 1

        # create new xarray object with this (possibly coarsened) subset:

        xs = x[i1:i2:coarsen]
        ys = y[j1:j2:coarsen]
        Zs = f.variables[zvar][j1:j2:coarsen, i1:i2:coarsen]

    if verbose:
        print('Returning a DEM with shape = %s' \
                % str(Zs.shape))
        print('x ranges from %.5f to %.5f with dx = %.8f' \
                % (xs[0], xs[-1], (xs[1]-xs[0])))
        print('y ranges from %.5f to %.5f with dy = %.8f' \
                % (ys[0], ys[-1], (ys[1]-ys[0])))
    output = None

    if return_topo:
        topo = topotools.Topography()
        topo._x = xs
        topo._y = ys
        topo._Z = Zs
        topo.generate_2d_coordinates()

        output = topo
            
    if return_xarray:
        # Create a new xarray.Dataset with this subsampled, coarsened data:
        dims = (len(xs),len(ys))
        xarray_ds = xarray.Dataset({'z':(dims,Zs)}, coords={'lon':xs, 'lat':ys})
        if output is None:
            output = xarray_ds
        else:
            output = (topo, xarray_ds)

    return output


