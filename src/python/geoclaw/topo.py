#!/usr/bin/env python
# encoding: utf-8

r"""
This module contains classes for representing topography files (topo) and time 
dependent topography files (dtopo) along with various operations one can perform
on or with them including reading, writing, transformations and plotting.

:TODO:
 - Implement tests (no testing really at all right now)
 - Implement topotype == 2
 - Implement subclass TimeDependentTopography

"""

import numpy

import matplotlib.pyplot as plt
import matplotlib.colors as colors 

import clawpack.visclaw.colormaps as colormaps
import clawpack.geoclaw.topotools as topotools

# ==============================================================================
#  New topography class
# ==============================================================================
class Topography(object):
    r"""
    Base topography class

    A class representing a single topography file.

    :Properties:

    :Initialization:
         - 

    :Examples:

        >>> import clawpack.geoclaw.topo as topo
        >>> topo_file = topo.Topography('./topo.tt3')
        >>> topo_file.plot()

    """

    @property
    def z(self):
        r"""A representation of the data as an 1d array."""
        if self._z is None:
            self.read()
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
            self.generate_2d_depths(mask=True)
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
            self.read(mask=True)
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
        if self._x is None:
            self.generate_2d_coordinates(mask=True)
        return self._X
    @X.deleter
    def X(self):
        del self._X

    @property
    def y(self):
        r"""One dimensional coordinate array in y direction."""
        if self._y is None:
            self.read(mask=True)
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
            self.generate_2d_coordinates(mask=True)
        return self._Y
    @Y.deleter
    def Y(self):
        del self._Y

    @property
    def extent(self):
        r"""Extents of the topography."""
        if self._extent is None:
            self._extent = ( numpy.min(self.x), numpy.max(self.x), 
                             numpy.min(self.y), numpy.max(self.y) )
        return self._extent
    @extent.setter
    def extent(self, value):
        self._extent = value

    @property
    def delta(self):
        r"""Spacing of data points"""
        if self._delta is None:
            if self.unstructured:

                # Calculate the smallest spacing between grid points            
                dx = numpy.infty
                dy = numpy.infty
                num_comparisons = self.x.shape[0] - 1
                for i in xrange(self.x.shape[0]):
                    for j in xrange(num_comparisons):
                        dx = min(dx, numpy.abs(self.x[i + j + 1] - self.x[i]))
                        dy = min(dy, numpy.abs(self.y[i + j + 1] - self.y[i]))

                    num_comparisons -= 1
                self._delta = [dx, dy]
            else:
                # All other topography types should have equally spaced grid
                # points in each direction
                self._delta = [self.x[1] - self.x[0], self.y[1] - self.y[0]]
                check_delta = [self.x[-2] - self.x[-1], self.y[-2] - self.y[-1]]
                assert self._delta[0] == check_delta,                  \
                       "Grid spacing delta not constant, %s != %s." %  \
                       (self._delta, check_delta)
        return self._delta


    def __init__(self, path, topo_type=3, unstructured=False):
        r"""
        Topography initialization routine
        
        See :class:`Topography` for more info.
        """

        super(Topography, self).__init__()

        self.path = path
        self.topo_type = topo_type
        self.unstructured = False

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


    def read(self, mask=True, filter_region=None):
        r"""
        Read in the data from the object's *path* attribute

        Stores the resulting data in one of the sets of *x*, *y*, and *z* or 
        *X*, *Y*, and *Z*.  

        :Input:
         - *mask* (bool)
         - *filter_region* (tuple)


        """

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
            if self.topo_type == 1:
                pass
            if self.topo_type == 2 or self.topo_type == 3:

                # Get header information
                N, extent, delta, no_data_value = self.read_header()
                self._x = numpy.linspace(extent[0],extent[1],N[0])
                self._y = numpy.linspace(extent[3],extent[2],N[1])

                if self.topo_type == 2:
                    raise NotImplemented("Topography type 2 unsupported" + \
                                         " right now.")
                elif self.topo_type == 3:
        
                    # Data is read in starting at the top right corner
                    self._Z = numpy.loadtxt(self.path, skiprows=6)
    
                if mask:
                    self._Z = numpy.ma.masked_values(self._Z, no_data_value, copy=False)


    def read_header(self):
        r"""Read in header of topography file at path.

        If a value returns numpy.nan then the value was not retrievable.
        """

        if self.topo_type != 2 or self.topo_type != 3:
            raise ValueError("The topography type must either be 2 or 3 to" + \
                             " read in a header.")

        # Default values to track errors
        num_cells = [numpy.nan,numpy.nan]
        extent = [numpy.nan,numpy.nan,numpy.nan,numpy.nan]
        delta = numpy.nan
        no_data_value = numpy.nan

        with open(self.path, 'r') as bathy_file:
            num_cells[0] = int(bathy_file.readline().split()[0])
            num_cells[1] = int(bathy_file.readline().split()[0])
            extent[0] = float(bathy_file.readline().split()[0])
            extent[2] = float(bathy_file.readline().split()[0])
            delta = float(bathy_file.readline().split()[0])
            no_data_value = float(bathy_file.readline().split()[0])
            
            extent[1] = extent[0] + num_cells[0] * delta
            extent[3] = extent[2] + num_cells[1] * delta

        return num_cells, extent, delta, no_data_value

    
    def generate_2d_depths(self, mask=True):
        r"""Generate a 2d array of the depths"""

        # Check to see if we need to generate these
        if self._Z is None:

            if self.unstructured:
                # Really no way to do this here with performing a projection via
                # extract.  Note that if the projection is performed these
                # arrays are already stored in self._X and self._Y
                raise ValueError("Unstructured data does not allow for use of" \
                                 + " 2d arrays, first project the data and" \
                                 + " try to perform this operation again.") 

            if self._z is None:
            # Try to read the data, may not have done this yet
                self.read(mask=mask)

            # See if self._X and self._Y are already computed and use them if
            # available, otherwise just use self._x and self._y
            if self._X is not None and self._Y is not None:
                new_shape = self._X.shape
            else:
                new_shape = (self._x.shape[0], self._y.shape[0])
            # Reshape, note that the mask follows along with the new array
            self._Z = numpy.reshape(self._z, new_shape)


    def generate_2d_coordinates(self, mask=True):
        r"""Generate 2d coordinate arrays."""

        # Check to see if we need to generate these
        if self._X is None and self._Y is None:

            if self.unstructured:
                # Really no way to do this here with performing a projection via
                # extract.  Note that if the projection is performed these
                # arrays are already stored in self._X and self._Y
                raise ValueError("Unstructured data does not allow for use of" \
                                 + " 2d coordinates, first project the data" \
                                 + " and try to perform this operation again.")

            if self._x is None or self._y is None:
                # Try to read the data to get these, may not have been done yet
                self.read(mask=mask)

            # Generate arrays
            self._X, self._Y = numpy.meshgrid(self._x, self._y)
            
            # If masking has been requested try to get the mask first from 
            # self._Z and then self._z
            if mask:
                if self._Z is None:
                    # Check to see if we really need to do anything here
                    if isinstance(self._z, numpy.ma.MaskedArray):
                        # Try to create self._Z
                        self.generate_2d_depths(mask=mask)

                if isinstance(self._Z, numpy.ma.MaskedArray):
                    # Use Z's mask for the X and Y coordinates
                    self._X = numpy.ma.MaskedArray(self._X, mask=self._Z.mask, 
                                                                     copy=False)
                    self._Y = numpy.ma.MaskedArray(self._Y, mask=self._Z.mask, 
                                                                     copy=False)


    def plot(self, axes=None, region_extent=None, contours=None, 
                   coastlines=True, limits=None, cmap=plt.get_cmap('terrain')):
        r"""
        Plot the topography
        """

        # Create axes if needed
        if axes is None:
            fig = plt.figure()
            axes = fig.add_subplot(111)
        
        # Turn off annoying offset
        axes.ticklabel_format(format="plain", useOffset=False)

        # Generate limits if need be
        if region_extent is None:
            region_extent = ( numpy.min(self.X), numpy.max(self.X),
                              numpy.min(self.Y), numpy.max(self.Y) )
        mean_lat = 0.5 * (region_extent[3] - region_extent[2])
        axes.set_aspect(1.0 / numpy.cos(numpy.pi / 180.0 * mean_lat))
        if limits is None:
            depth_extent = (numpy.min(self.Z),numpy.max(self.Z))
        else:
            depth_extent = limits

        # Create color map
        cmap = colormaps.make_colormap({-1:[0.3,0.2,0.1],
                                           -0.00001:[0.95,0.9,0.7],
                                           0.00001:[.5,.7,0],
                                           1:[.2,.5,.2]})
        color_norm = colors.Normalize(depth_extent[0],depth_extent[1],clip=True)

        # Plot data
        if contours is not None:
            plot = axes.contourf(self.X, self.Y, self.Z, contours,cmap=cmap)
        elif isinstance(self.Z, numpy.ma.MaskedArray):
            plot = axes.pcolor(self.X, self.Y, self.Z, vmin=depth_extent[0], 
                                                       vmax=depth_extent[1],
                                                       cmap=cmap, 
                                                       norm=color_norm)
        else:
            plot = axes.imshow(self.Z, vmin=depth_extent[0], 
                                       vmax=depth_extent[1],
                                       extent=region_extent, 
                                       cmap=cmap, 
                                       norm=color_norm)
        cbar = plt.colorbar(plot, ax=axes)
        cbar.set_label("Depth (m)")
        # levels = range(0,int(-numpy.min(Z)),500)

        # Plot coastlines
        if coastlines:
            axes.contour(self.X, self.Y, self.Z, levels=[0.0],colors='r')

        axes.set_xlim(region_extent[0:2])
        axes.set_ylim(region_extent[2:])

        return axes


    def write(self, path, no_data_value=999999, topotype=3):
        r"""
        Write out a topography file to path of type topotype

        Writes out a bathymetry file of type 3 to path from data in Z.  The rest of
        the arguments are used to write the header data.
        """

        with open(path,'w') as outfile:
            if self.unstructured:
                for (i, depth) in enumerate(self.z):
                    outfile.write("%s %s %s\n" % (self.x[i], self.y[i], depth))

            elif topotype == 1:
                # longitudes = numpy.linspace(lower[0], lower[0] + delta * Z.shape[0], Z.shape[0])
                # latitudes = numpy.linspace(lower[1], lower[1] + delta * Z.shape[1], Z.shape[1])
                for (j, latitude) in enumerate(self.y):
                    for (i, longitude) in enumerate(self.x):
                        outfile.write("%s %s %s\n" % (longitude, latitude, self.Z[i,j]))

            elif topotype == 2 or topotype == 3:

                if topotype == 2:
                    raise NotImplemented("Writing topo type = 2 is not implemented.")

                # Write out header
                outfile.write('%s ncols\n' % self.Z.shape[1])
                outfile.write('%s nrows\n' % self.Z.shape[0])
                outfile.write('%s xll\n' % self.x[0])
                outfile.write('%s yll\n' % self.y[0])
                outfile.write('%s cellsize\n' % self.delta)
                outfile.write('%s nodata_value\n' % no_data_value)

                # Write out bathy data
                # We flip the output data here since we write from the upper left corner
                # to lower right and the data is ordered from lower left to upper right
                Z_flipped = numpy.flipud(self.Z)
                for i in xrange(self.Z.shape[0]):
                    for j in xrange(self.Z.shape[1]):
                        outfile.write("%s   " % (Z_flipped[i,j]))
                    outfile.write("\n")

            else:
                raise NotImplemented("Output type %s not implemented." % topotype)


    def project_unstructured(self, X_fill, Y_fill, Z_fill, extent=None,
                                   method='nearest', delta_limit=20.0, 
                                   no_data_value=-9999, buffer_length=100.0,
                                   proximity_radius=100.0, 
                                   resolution_limit=2000):
        r"""
        Project unstructured data on to regular grid

        Function to project the unstructured data in the topo object onto a 
        structured grid.  Utilized a bounding box plus a buffer of size 
        *buffer_length* (meters) containing all data unless *extent* is not 
        None.  Then uses the fill data provided (*X_fill*, *Y_fill* and 
        *Z_fill*) to fill in the gaps in the unstructured data.  By default this
        is done by masking the fill data with the extents, the value 
        *no_data_value* and if *proximity_radius* (meters) is not 0, by a radius
        of *proximity_radius* from all grid points in the object.  Stores the 
        result in the *self.X*, *self.Y* and *self.Z* object attributes.  The
        resolution of the final grid is determined by calculating the minimum
        distance between all *x* and *y* data with a hard lower limit of 
        *delta_limit* (meters).

        :Input:
         - *extent* (tuple) - A tuple defining the rectangle of the sub-section.  
           Must be in the form (x lower,x upper,y lower, y upper).
         - *method* (string) - Method used for interpolation, valid methods are
           found in *scipy.interpolate.griddata*.  Default is *nearest*.
         - *delta_limit* (float) - Limit of finest horizontal resolution, 
           default is 20 meters.
         - *no_data_value* (float) - Value to use if no data was found to fill in a 
           missing value, ignored if `method = 'nearest'`. Default is `-9999`.
         - *buffer_length* (float) - Buffer around bounding box, only applicable
           when *extent* is None.  Default is *100.0* meters.
         - *proximity_radius* (float) - Radius every unstructured data point
           used to mask the fill data with.  Default is *100.0* meters.
         - *resolution_limit* (int) - Limit the number of grid points in a
           single dimension.  Raises a *ValueError* if the limit is violated.
           Default value is 

        """

        import scipy.interpolate as interpolate

        mean_lat = 0.5 * (self.y[0] + self.y[-1])
            
        # Calculate new grid coordinates
        if extent is None:
            buffer_degrees = topotools.latlong_resolution(buffer_length, 0.0, mean_lat)
            extent = [ numpy.min(self.x) - buffer_degrees, 
                       numpy.max(self.x) + buffer_degrees, 
                       numpy.min(self.y) - buffer_degrees, 
                       numpy.max(self.y) + buffer_degrees ]
        delta = max( min(numpy.min(numpy.abs(self.x[1:] - self.x[:-1])), 
                         numpy.min(numpy.abs(self.y[1:] - self.y[:-1])) ), 
                     topotools.latlong_resolution(delta_limit, 0.0, mean_lat) )
        N = ( numpy.ceil((extent[1] - extent[0]) / delta),
              numpy.ceil((extent[3] - extent[2]) / delta) )
        assert numpy.all(N[:] < numpy.ones((2)) * resolution_limit), \
               ValueError("Calculated resolution too high, N=%s!" % str(N))
        self._X, self._Y = numpy.meshgrid( 
                                     numpy.linspace(extent[0], extent[1], N[0]),
                                     numpy.linspace(extent[2], extent[3], N[1]))

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
            proximity_radius_deg = \
              topotools.latlong_resolution(proximity_radius, 0.0, mean_lat)
        
            indices = (~all_mask).nonzero()
            for n in xrange(indices[0].shape[0]):
                i = indices[0][n]
                j = indices[1][n]
                all_mask[i,j] = numpy.any(numpy.sqrt((self.x - X_fill[i,j])**2 
                                                   + (self.y - Y_fill[i,j])**2)
                                             < proximity_radius_deg)

        X_fill_masked = numpy.ma.masked_where(all_mask, X_fill)
        Y_fill_masked = numpy.ma.masked_where(all_mask, Y_fill)
        Z_fill_masked = numpy.ma.masked_where(all_mask, Z_fill)    

        # Stick both the input data and fill data into arrays
        fill_points = numpy.column_stack((X_fill_masked.compressed(),
                                          Y_fill_masked.compressed()))
        points = numpy.concatenate((numpy.array([self.x, self.y]).transpose(), 
                                    fill_points))
        values = numpy.concatenate((self.z, Z_fill_masked.compressed()))

        # Use nearest-neighbor interpolation
        self._Z = interpolate.griddata(points, values, (self.X, self.Y), 
                                                                  method=method)

        self._extent = extent
        self._delta = delta


class TimeDependentTography(Topography):
    r"""
    Class represents a patch of time dependent topography (dtopo)

    """

    def __init__(self):
        r"""
        TimeDepedentTopography initialization routine
        
        See :class:`TimeDependentTography` for more info.
        """

        super(TimeDependentTography, self).__init__()



if __name__ == "__main__":
    pass