#!/usr/bin/env python
# encoding: utf-8

r"""
GeoClaw dtopotools Module  `$CLAW/geoclaw/src/python/geoclaw/dtopotools.py`

Module provides several functions for dealing with changes to topography (usually
due to earthquakes) including reading sub-fault specifications, writing out 
dtopo files, and calculating Okada based deformations.

:Classes:

  - DTopography
  - SubFault
  - Fault
  - UCSBFault
  - CSVFault
  - SiftFault
  - SegmentedPlaneFault
    
:Functions:
  - convert_units
  - plot_dz_contours
  - plot_dz_colors
  - Mw
  - strike_direction
  - rise_fraction

"""

from __future__ import absolute_import
from __future__ import print_function
import os
import sys
import re

import numpy

import clawpack.geoclaw.topotools as topotools
import clawpack.geoclaw.util as util

# ==============================================================================
#  Constants
# ==============================================================================
from clawpack.geoclaw.data import DEG2RAD, LAT2METER
import six
from six.moves import range

# Poisson ratio for Okada 
poisson = 0.25

# ==============================================================================
#  Units dictionaries
# ==============================================================================

# Dictionary for standard units to be used for all subfault models.
# The data might be read in from a file where different units are used,
# in which case the *input_units* argument of the *read* method can be used
# to indicate these units.  The *read* function should then convert to these
# standard units:
standard_units = {}
standard_units['length'] = 'm'
standard_units['width'] = 'm'
standard_units['depth'] = 'm'
standard_units['slip'] = 'm'
standard_units['mu'] = 'Pa'

# Dictionary for converting input_units specified by user to or from
# standard units used internally:
# (Conversion is performed by the module function *convert_units*, which
# is called by *SubFault.convert_to_standard_units*)

unit_conversion_factor = {}  
# for length, width, depth, slip:  (standard units = 'm')
unit_conversion_factor['m'] = 1.
unit_conversion_factor['cm'] = 0.01
unit_conversion_factor['km'] = 1000.
unit_conversion_factor['nm'] = 1852.0  # nautical miles
# for rigidity (shear modulus) mu:  (standard units = 'Pa')
unit_conversion_factor['Pa'] = 1.
unit_conversion_factor['GPa'] = 1.e9
unit_conversion_factor['dyne/cm^2'] = 0.1
unit_conversion_factor['dyne/m^2'] = 1.e-5
# for seismic moment Mo:  (standard units = 'N-m', Newton-meters)
unit_conversion_factor['N-m'] = 1.
unit_conversion_factor['dyne-cm'] = 1.e-7

# Check that these are consistent:
check = [unit_conversion_factor[standard_units[param]] == 1. for param in \
         standard_units.keys()]
if not numpy.alltrue(check):
    raise ValueError("Conversion factors should be 1 for all standard_units")


# ==============================================================================
#  General utility functions
# ==============================================================================
def convert_units(value, io_units, direction=1, verbose=False):
    r"""
    convert *value* to standard units from *io_units* or vice versa.
    *io_units* (str) refers to the units used in the subfault file read or to be
    written.  The standard units are those used internally in this module.
    See the comments below for the standard units.
    If *direction==1*, *value* is in *io_units* and convert to standard.
    If *direction==2*, *value* is in standard units and convert to *io_units*.

    """
    try:
        factor = unit_conversion_factor[io_units]
    except:
        factor = 1.
        print("*** Warning: unrecoginized units in convert_units, not converting")
        #raise ValueError("Unrecognized io_units %s, must be one of %s" \
        #      % (io_units, unit_conversion_factor.keys()))
    if direction == 1:
        converted_value = value * factor
    elif direction == 2:
        converted_value = value / factor
    else:
        raise ValueError("Unrecognized direction, must be 1 or 2")

    return converted_value


def plot_dZ_contours(x, y, dZ, axes=None, dZ_interval=0.5, verbose=False,
                               fig_kwargs={}):
    r"""For plotting seafloor deformation dZ"""
    import matplotlib.pyplot as plt

    dZ_max = max(dZ.max(), -dZ.min()) + dZ_interval
    clines1 = numpy.arange(dZ_interval, dZ_max, dZ_interval)
    clines = list(-numpy.flipud(clines1)) + list(clines1)

    # Create axes if needed
    if axes is None:
        fig = plt.figure(**fig_kwargs)
        axes = fig.add_subplot(111)

    if len(clines) > 0:
        if verbose:
            print("Plotting contour lines at: ",clines)
        axes.contour(x, y, dZ, clines, colors='k')
    else:   
        print("No contours to plot")

    return axes


def plot_dZ_colors(x, y, dZ, axes=None, cmax_dZ=None, dZ_interval=None,
                   add_colorbar=True, verbose=False, fig_kwargs={}):
    r"""
    Plot sea floor deformation dZ as colormap with contours
    """

    from clawpack.visclaw import colormaps
    import matplotlib.pyplot as plt

    if axes is None:
        fig = plt.figure(**fig_kwargs)
        axes = fig.add_subplot(1, 1, 1)
    #print "+++ in plot_dZ_colors, axes = ",axes
    #print "+++ in plot_dZ_colors, id(axes) = ",id(axes)

    dZmax = numpy.abs(dZ).max()
    if cmax_dZ is None:
        if dZmax < 1.e-12:
            cmax_dZ = 0.1
        else:
            cmax_dZ = dZmax
    cmap = colormaps.blue_white_red
    extent = [x.min(), x.max(), y.min(), y.max()]
    im = axes.imshow(dZ, extent=extent, cmap=cmap, origin='lower')
    im.set_clim(-cmax_dZ,cmax_dZ)
    if add_colorbar:
        cbar = plt.colorbar(im, ax=axes)
        cbar.set_label("Deformation (m)")
    

    if dZ_interval is None:
        dZ_interval = cmax_dZ/10.
    clines1 = numpy.arange(dZ_interval, dZmax + dZ_interval, dZ_interval)
    clines = list(-numpy.flipud(clines1)) + list(clines1)
    if len(clines) > 0:
        if verbose:
            print("Plotting contour lines at: ",clines)
        axes.contour(x,y,dZ,clines,colors='k',linestyles='solid')
    elif verbose:
        print("No contours to plot")

    y_ave = 0.5 * (y.min() + y.max())
    axes.set_aspect(1. / numpy.cos(y_ave * numpy.pi / 180.))
    axes.ticklabel_format(format='plain', useOffset=False)
    axes.set_title('Seafloor deformation')
    for label in axes.get_xticklabels():
        label.set_rotation(20)

    return axes


def Mw(Mo, units="N-m"):
    """ 
    Calculate moment magnitude based on seismic moment Mo.
    Follows USGS recommended definition from 

        http://earthquake.usgs.gov/aboutus/docs/020204mag_policy.php

    The SubFault and Fault classes each have a function Mo to compute
    the seismic moment for a single subfault or collection respectively.
    """

    if units == "N-m":
        Mw = 2/3.0 * (numpy.log10(Mo) - 9.05)
    elif units == "dyne-cm":
        Mw = 2/3.0 * numpy.log10(Mo) - 10.7
        #  = 2/3.0 * (numpy.log10(1e-7 * Mo) - 9.05)
    else:
        raise ValueError("Unknown unit for Mo: %s." % units)

    return Mw

    
def strike_direction(x1, y1, x2, y2):
    """
    Calculate strike direction between two points.
    Actually calculates "initial bearing" from (x1,y1) in direction
    towards (x2,y2), following
    http://www.movable-type.co.uk/scripts/latlong.html
    """

    x1 = x1*numpy.pi/180.
    y1 = y1*numpy.pi/180.
    x2 = x2*numpy.pi/180.
    y2 = y2*numpy.pi/180.
    dx = x2-x1
    theta = numpy.arctan2(numpy.sin(dx)*numpy.cos(y2), \
            numpy.cos(y1)*numpy.sin(y2) \
            - numpy.sin(y1)*numpy.cos(y2)*numpy.cos(dx))
    s = theta*180./numpy.pi
    if s<0:
        s = 360+s
    return s


def rise_fraction(t, t0, t_rise, t_rise_ending=None):
    """
    A continuously differentiable piecewise quadratic function of t that is 

     *  0 for t <= t0, 
     *  1 for t >= t0 + t_rise + t_rise_ending 

    with maximum slope at t0 + t_rise.
    For specifying dynamic fault ruptures:  Subfault files often contain these
    parameters for each subfault for an earthquake event.

    *t* can be a scalar or a numpy array of times and the returned result
    will have the same type.  A list or tuple of times returns a numpy array.
    """

    scalar = (type(t) in [float,int])
    t = numpy.array(t)

    if t_rise_ending is None: 
        t_rise_ending = t_rise

    t1 = t0+t_rise
    t2 = t1+t_rise_ending

    rf = numpy.where(t<=t0, 0., 1.)
    if t2 != t0:

        t20 = float(t2-t0)
        t10 = float(t1-t0)
        t21 = float(t2-t1)

        c1 = t21 / (t20*t10*t21) 
        c2 = t10 / (t20*t10*t21) 

        rf = numpy.where((t>t0) & (t<=t1), c1*(t-t0)**2, rf)
        rf = numpy.where((t>t1) & (t<=t2), 1. - c2*(t-t2)**2, rf)

    if scalar:
        rf = float(rf)   # return a scalar if input t is scalar

    return rf



# ==============================================================================
#  DTopography Base Class
# ==============================================================================
class DTopography(object):
    r"""Basic object representing moving topography



    """


    def __init__(self, path=None, dtopo_type=None):
        r"""DTopography initialization routine.
        
        See :class:`DTopography` for more info.

        """

        self.dZ = None
        self.times = []
        self.x = None
        self.y = None
        self.X = None
        self.Y = None
        self.delta = None
        self.path = path
        if path:
            self.read(path, dtopo_type)


    def read(self, path=None, dtopo_type=None, verbose=False):
        r"""
        Read in a dtopo file and use to set attributes of this object.

        :input:
        
         - *path* (path) - Path to existing dtopo file to read in.
         - *dtopo_type* (int) - Type of topography file to read.  Default is 3
            if not specified or apparent from file extension.
        """

        if path is not None:
            self.path = path
        else:
            if self.path is None:
                raise ValueError("Need to specify a path to a file.")
            else:
                path = self.path

        if dtopo_type is None:
            dtopo_type = topotools.determine_topo_type(path, default=3)

        if dtopo_type == 1:
            data = numpy.loadtxt(path)
            if verbose:
                print("Loaded file %s with %s lines" %(path,data.shape[0]))
            t = list(set(data[:,0]))
            t.sort()
            if verbose:
                print("times found: ",t)
            ntimes = len(t)
            tlast = t[-1]
            lastlines = data[data[:,0]==tlast]
            xvals = list(set(lastlines[:,1]))
            xvals.sort()
            mx = len(xvals)
            my = len(lastlines) / mx
            if verbose:
                print("Read dtopo: mx=%s and my=%s, at %s times" % (mx,my,ntimes))
            X = numpy.reshape(lastlines[:,1],(my,mx))
            Y = numpy.reshape(lastlines[:,2],(my,mx))
            Y = numpy.flipud(Y)
            if verbose:
                print("Returning dZ as a list of mx*my arrays")
            dZ = None
            for n in range(ntimes):
                i1 = n*mx*my
                i2 = (n+1)*mx*my
                dzt = numpy.reshape(data[i1:i2,3],(my,mx))
                dzt = numpy.flipud(dzt)
                dzt = numpy.array(dzt, ndmin=3)  # convert to 3d array
                if dZ is None:
                    dZ = dzt.copy()
                else:
                    dZ = numpy.append(dZ, dzt, axis=0)
            self.X = X
            self.Y = Y
            self.x = X[0,:]
            self.y = Y[:,0]
            self.times = t
            self.dZ = dZ

        elif dtopo_type == 2 or dtopo_type == 3:
            fid = open(path)
            mx = int(fid.readline().split()[0])
            my = int(fid.readline().split()[0])
            mt = int(fid.readline().split()[0])
            xlower = float(fid.readline().split()[0])
            ylower = float(fid.readline().split()[0])
            t0 = float(fid.readline().split()[0])
            dx = float(fid.readline().split()[0])
            dy = float(fid.readline().split()[0])
            dt = float(fid.readline().split()[0])
            fid.close()
    
            xupper = xlower + (mx-1)*dx
            yupper = ylower + (my-1)*dy
            x=numpy.linspace(xlower,xupper,mx)
            y=numpy.linspace(ylower,yupper,my)
            times = numpy.linspace(t0, t0+(mt-1)*dt, mt)
    
            dZvals = numpy.loadtxt(path, skiprows=9)
            if dtopo_type==3:
                # my lines with mx values on each
                for k,t in enumerate(times):
                    dZk = numpy.reshape(dZvals[k*my:(k+1)*my, :], (my,mx))
                    dZk = numpy.flipud(dZk)
                    dZk = numpy.array(dZk, ndmin=3)  # convert to 3d array
                    if k==0:
                        dZ = dZk.copy()
                    else:
                        dZ = numpy.append(dZ, dZk, axis=0)
            else:
                # dtopo_type==2 ==> mx*my lines with 1 values on each
                for k,t in enumerate(times):
                    dZk = numpy.reshape(dZvals[k*mx*my:(k+1)*mx*my], (my,mx))
                    dZk = numpy.flipud(dZk)
                    dZk = numpy.array(dZk, ndmin=3)  # convert to 3d array
                    if k==0:
                        dZ = dZk.copy()
                    else:
                        dZ = numpy.append(dZ, dZk, axis=0)
                    
            self.x = x
            self.y = y
            self.X, self.Y = numpy.meshgrid(x,y)
            self.times = times
            self.dZ = dZ

        else:
            raise ValueError("Only topography types 1, 2, and 3 are supported,",
                             " given %s." % dtopo_type)


    def write(self, path=None, dtopo_type=None):
        r"""Write out subfault resulting dtopo to file at *path*.

        :input:
        
         - *path* (path) - Path to the output file to written to.
         - *dtopo_type* (int) - Type of topography file to write out.  Default
           is 3.

        """

        if path is not None:
            self.path = path
        if self.path is None:
            raise IOError("*** need to specify path to file for writing")
        path = self.path

        if dtopo_type is None:
            dtopo_type = topotools.determine_topo_type(path, default=3)

        x = self.X[0,:]
        y = self.Y[:,0]
        dx = x[1] - x[0]
        dy = y[1] - y[0]

        # This is no longer required in GeoClaw...
        #if abs(dx - dy) >= 1e-12:
        #    raise ValueError("dx = %g not equal to dy = %g" % (dx,dy))

        # Construct each interpolating function and evaluate at new grid
        ## Shouldn't need to interpolate in time.
        with open(path, 'w') as data_file:

            if dtopo_type == 0:
                # Topography file with 3 columns, x, y, dz written from the
                # upper left corner of the region. Only final time.
                Y_flipped = numpy.flipud(self.Y)
                dZ_flipped = numpy.flipud(self.dZ[-1,:,:]) 

                for j in range(self.Y.shape[0]):
                    for i in range(self.X.shape[1]):
                        data_file.write("%s %s %s\n" % self.X[j,i], 
                            Y_flipped[j,i], dZ_flipped[j,i])

            elif dtopo_type == 1:
                # Topography file with 4 columns, t, x, y, dz written from the
                # upper
                # left corner of the region
                Y_flipped = numpy.flipud(self.Y)
                for (n, time) in enumerate(self.times):
                    #alpha = (time - self.t[0]) / self.t[-1]
                    #dZ_flipped = numpy.flipud(alpha * self.dZ[:,:])
                    dZ_flipped = numpy.flipud(self.dZ[n,:,:])  

                    for j in range(self.Y.shape[0]):
                        for i in range(self.X.shape[1]):
                            data_file.write("%s %s %s %s\n" % (self.times[n],
                                self.X[j,i], Y_flipped[j,i], dZ_flipped[j,i]))
        
            elif dtopo_type == 2 or dtopo_type == 3:
                if len(self.times) == 1:
                    dt = 0.
                else:
                    dt = float(self.times[1] - self.times[0])
                # Write out header
                data_file.write("%7i       mx \n" % x.shape[0])
                data_file.write("%7i       my \n" % y.shape[0])
                data_file.write("%7i       mt \n" % len(self.times))
                data_file.write("%20.14e   xlower\n" % x[0])
                data_file.write("%20.14e   ylower\n" % y[0])
                data_file.write("%20.14e   t0\n" % self.times[0])
                data_file.write("%20.14e   dx\n" % dx)
                data_file.write("%20.14e   dy\n" % dy)
                data_file.write("%20.14e   dt\n" % dt)

                if dtopo_type == 2:
                    raise ValueError("Topography type 2 is not yet supported.")
                elif dtopo_type == 3:
                    for (n, time) in enumerate(self.times):
                        #alpha = (time - self.t[0]) / (self.t[-1])
                        for j in range(self.Y.shape[0]-1, -1, -1):
                            data_file.write(self.X.shape[1] * '%012.6e  ' 
                                                  % tuple(self.dZ[n,j,:]))
                            data_file.write("\n")

            else:
                raise ValueError("Only topography types 1, 2, and 3 are ",
                                 "supported, given %s." % dtopo_type)


    def dZ_at_t(self, t):
        """
        Interpolate dZ to specified time t and return deformation.
        """
        from matplotlib.mlab import find
        if t <= self.times[0]:
            return self.dZ[0,:,:]
        elif t >= self.times[-1]:
            return self.dZ[-1,:,:]
        else:
            n = max(find(self.times <= t))
            t1 = self.times[n]
            t2 = self.times[n+1]
            dz = (t2-t)/(t2-t1) * self.dZ[n,:,:] + \
                 (t-t1)/(t2-t1) * self.dZ[n+1,:,:]
            return dz


    def dZ_max(self):
        r"""Return max(abs(dZ)) over all dz in self.dZ, the maximum
        surface deformation for this dtopo.
        DEPRECATE? -- it's now a 1-liner
        """

        return abs(self.dZ).max()


    def plot_dZ_colors(self, t, axes=None, cmax_dZ=None, dZ_interval=None, 
                                fig_kwargs={}):
        """
        Interpolate self.dZ to specified time t and then call module function
        plot_dZ_colors.
        """
        axes = plot_dZ_colors(self.X, self.Y, self.dZ_at_t(t), axes=axes,
                              cmax_dZ=cmax_dZ, dZ_interval=dZ_interval,
                              fig_kwargs=fig_kwargs)
        return axes


    def plot_dZ_contours(self, t, dZ_interval=0.5, axes=None, fig_kwargs={}):
        """
        Interpolate self.dZ to specified time t and then call module function
        plot_dZ_contours.
        """
        axes = plot_dZ_contours(self.X, self.Y, self.dZ_at_t(t), axes=axes,
                                dZ_interval=dZ_interval)
        return axes


                

# ==============================================================================
#  Generic Fault Class
# ==============================================================================
class Fault(object):
    
    r"""Base Fault class

    A class describing a fault possibly composed of subfaults.

    :Properties:

    :Initialization:

    :Examples:

    """

    def __init__(self, subfaults=None, input_units={},
                 coordinate_specification=None):
        r"""Fault initialization routine.
        
        See :class:`Fault` for more info.

        """

        # Parameters for subfault specification
        self.rupture_type = 'static' # 'static' or 'dynamic'
        #self.times = numpy.array([0., 1.])   # or just [0.] ??
        self.dtopo = None

        # Default units of each parameter type
        self.input_units = standard_units.copy()
        self.input_units.update(input_units)

        # Set the coordinate specification, e.g. 'top center':
        self.coordinate_specification = coordinate_specification
        
        if subfaults is not None:
            if not isinstance(subfaults, list):
                raise ValueError("Input parameter subfaults must be a list.")
            self.subfaults = subfaults
            for subfault in self.subfaults:
                subfault.convert_to_standard_units(self.input_units)
                if subfault.coordinate_specification is None:
                    subfault.coordinate_specification = coordinate_specification
                if subfault.coordinate_specification is None:
                    raise ValueError("Must specify coordinate_specification, " + \
                            "either for fault or for each subfault")


    def read(self, path, column_map, coordinate_specification="centroid",
                                     rupture_type="static", skiprows=0, 
                                     delimiter=None, input_units={}, defaults=None):
        r"""Read in subfault specification at *path*.

        Creates a list of subfaults from the subfault specification file at
        *path*.
    
        :Inputs:

          - *path* (str) file to read in, should contain subfaults, one per line
          - *column_map* (dict) specifies mapping from parameter to the column
            of the input file that contains values for this parameter, e.g.

                column_map = {"latitude":0, "longitude":1, "depth":2, "slip":3,
                "rake":4, "strike":5, "dip":6}

          - *coordinate_specification* (str) specifies the location on each
            subfault that corresponds to the (longitude,latitude) and depth 
            of the subfault.  See the documentation for 
            *SubFault.calculate_geometry*.
          - *rupture_type* (str) either "static" or "dynamic"
          - *skiprows* (int) number of header lines to skip before data
          - *delimiter* (str) e.g. ',' for csv files
          - *input_units* (dict) indicating units for length, width, slip, depth,
                           and for rigidity mu as specified in file.  These
                           will be converted to "standard units".
          - *defaults* (dict) default values for all subfaults, for values not
                       included in subfault file on each line.

        """

        # Read in rest of data
        # (Use genfromtxt to deal with files containing strings, e.g. unit
        # source name, in some column)
        data = numpy.genfromtxt(path, skip_header=skiprows, delimiter=delimiter)
        if len(data.shape) == 1:
            data = numpy.array([data])

        self.coordinate_specification = coordinate_specification
        self.input_units = standard_units.copy()
        self.input_units.update(input_units)
        self.subfaults = []
        for n in range(data.shape[0]):

            new_subfault = SubFault()
            new_subfault.coordinate_specification = coordinate_specification
            
            for (var, column) in six.iteritems(column_map):
                if isinstance(column, tuple) or isinstance(column, list):
                    setattr(new_subfault, var, [None for k in column])
                    for (k, index) in enumerate(column):
                        getattr(new_subfault, var)[k] = data[n, index]
                else:
                    setattr(new_subfault, var, data[n, column])

            if defaults is not None:
                for param in six.iterkeys(defaults):
                    setattr(new_subfault, param, defaults[param]) 

            new_subfault.convert_to_standard_units(self.input_units)
            self.subfaults.append(new_subfault)


    def write(self, path, style=None, column_list=None, output_units={}, 
                    delimiter='  '):
        r"""
        Write subfault format file with one line for each subfault.
        Can either specify a *style* that determines the columns, 
        or a *column_list*.  Must specify one but not both.  See below for
        details.

        Inputs:
          - *path* (str) file to write to. 
          - *style* (str) to write in a style that matches standard styles
            adopted by various groups.  One of the following:

              - "usgs"  (Not implemented)
              - "noaa sift"  (Not implemented)
              - "ucsb"  (Not implemented)

          - *column_list* (list) specifies what order the parameters should
            be written in the output file, e.g.

                column_list = ['longitude','latitude','length','width',
                'depth','strike','rake','dip','slip']

          - *output_units* (dict) specifies units to convert to before writing.
            Defaults to "standard units".
          - *delimiter* (str) specifies delimiter between columns, e.g.
            "," to create a csv file.  Defaults to "  ".

        """

        self.output_units = standard_units.copy()
        self.output_units.update(output_units)

        if style is not None:
            msg =  "style option not yet implemented, use column_list"
            raise NotImplementedError(msg)

        if column_list is None:
            raise Exception("Must specify column_list")
        
        format = {}
        format['longitude'] = '%15.5f'
        format['latitude'] = '%15.5f'
        format['strike'] = '%15.5f'
        format['rake'] = '%15.5f'
        format['dip'] = '%15.5f'
        format['depth'] = '%15.8e'
        format['length'] = '%15.8e'
        format['width'] = '%15.8e'
        format['slip'] = '%15.8e'

        with open(path, 'w') as data_file:
            c_s_list = set([s.coordinate_specification for s in self.subfaults])
            if (len(c_s_list) >= 1) and \
                    (c_s_list.pop() != self.coordinate_specification):
                raise ValueError("Subfaults do not have common " +
                    "coordinate_specification that agrees with fault attribute")
            # write header:
            data_file.write('Subfaults file with coordinate_specification:  ')
            data_file.write('%s, \n' % self.coordinate_specification)
            data_file.write('Units: %s, \n' % str(output_units))
            s = ""
            for param in column_list:
                s = s + param.rjust(15) + delimiter
            data_file.write(s + '\n')
            for subfault in self.subfaults:
                s = ""
                for param in column_list:
                    value = getattr(subfault,param)
                    if param in output_units:
                        converted_value = convert_units(value, 
                                    self.output_units[param], direction=2)
                    s = s + format[param] % value + delimiter 
                data_file.write(s + '\n')
                


    def Mo(self):
        r""" 
        Calculate the seismic moment for a fault composed of subfaults,
        in units N-m.
        """

        total_Mo = 0.0
        for subfault in self.subfaults:
            total_Mo += subfault.Mo()
        return total_Mo


    def Mw(self):
        r"""Calculate the moment magnitude for a fault composed of subfaults."""
        return Mw(self.Mo())

    
    def create_dtopography(self, x, y, times=[0., 1.], verbose=False):
        r"""Compute change in topography and construct a dtopography object.

        Use subfaults' `okada` routine and add all 
        deformations together.

        Raises a ValueError exception if the *rupture_type* is an unknown type.

        returns a :class`DTopography` object.
        """

        dtopo = DTopography()
        dtopo.x = x
        dtopo.y = y
        X, Y = numpy.meshgrid(x,y)
        dtopo.X = X
        dtopo.Y = Y
        dtopo.times = times

        if verbose:
            print("Making Okada dz for each of %s subfaults" \
                  % len(self.subfaults))

        for k,subfault in enumerate(self.subfaults):
            if verbose:
                sys.stdout.write("%s.." % k)
                sys.stdout.flush()
            subfault.okada(x,y)  # sets subfault.dtopo with times=[0]
                                 # and subfault.dtopo.dZ.shape[0] == 1
        if verbose:
            sys.stdout.write("\nDone\n")

        if self.rupture_type == 'static':
            if len(times) > 2:
                raise ValueError("For static deformation, need len(times) <= 2")
            dz = numpy.zeros(X.shape)
            for subfault in self.subfaults:
                dz += subfault.dtopo.dZ[0,:,:]

            if len(times) == 1:
                # only final deformation stored:
                dtopo.dZ = numpy.array(dz, ndmin=3) 
            elif len(times) == 2:
                # store 0 at first time and final deformation at second:
                dz0 = numpy.zeros(X.shape)
                dtopo.dZ = numpy.array([dz0, dz])
                if dtopo.dZ.shape != (2, dz.shape[0], dz.shape[1]):
                    raise ValueError("dtopo.dZ does not have expected shape")

        elif self.rupture_type in ['dynamic','kinematic']:

            t_prev = -1.e99
            dzt = numpy.zeros(X.shape)
            dZ = None
            for t in times:
                for k,subfault in enumerate(self.subfaults):
                    t0 = getattr(subfault,'rupture_time',0)
                    t1 = getattr(subfault,'rise_time',0.5)
                    t2 = getattr(subfault,'rise_time_ending',None)
                    rf = rise_fraction([t_prev,t],t0,t1,t2)
                    dfrac = rf[1] - rf[0]
                    if dfrac > 0.:
                        dzt = dzt + dfrac * subfault.dtopo.dZ[0,:,:]
                dzt = numpy.array(dzt, ndmin=3)  # convert to 3d array
                if dZ is None:
                    dZ = dzt.copy()
                else:
                    dZ = numpy.append(dZ, dzt, axis=0)
                t_prev = t
            dtopo.dZ = dZ

        else:   
            raise ValueError("Unrecognized rupture_type: %s" % self.rupture_type)

        # Store for user
        self.dtopo = dtopo

        return dtopo


    
    def plot_subfaults(self, axes=None, plot_centerline=False, slip_color=False,
                             cmap_slip=None, cmin_slip=None, cmax_slip=None,
                             slip_time=None, plot_rake=False, xylim=None, 
                             plot_box=True, colorbar_shrink=1, verbose=False):
        """
        Plot each subfault projected onto the surface.

        *axes* can be passed in to specify the *matplotlib.axes.AxesSubplot*
        on which to add this plot.  If *axes == None*, a new figure window
        will be opened.  The *axes* on which it is plotted is the return
        value of this call.

        If *plot_centerline == True*, plot a line from the centroid to the
        top center of each subfault to show what direction is up-dip.

        If *slip_color == True* then use the color map *cmap_slip* 
        (which defaults to *matplotlib.cm.jet*) to color the subplots based
        on the magnitude of slip, scaled between *cmin_slip* and *cmax_slip*.  
        (If these are *None* then scaled automatically based on range of slip.)
        If *slip_time == None* then colors are based on the final slip.
        For dynamic faults, *slip_time* can be set to a time and the 
        dynamic timing of each subfault will be used to compute and 
        plot the slip at this time.

        If *plot_rake == True*, plot a line from the centroid pointing in
        the direction of the rake (the direction in which the top block is
        moving relative to the lower block.  The distance it moves is given
        by the *slip*.)

        *xylim* can be set to a list or tuple of length 4 of the form
        [x1,x2,y1,y2] to specify the x- and y-axis limits.

        If *plot_box == True*, a box will be drawn around each subfault.
        """
    
        import matplotlib
        import matplotlib.pyplot as plt

        if (slip_time is not None) and (self.rupture_type == 'static'):
            raise Exception("slip_time can only be specified for dynamic faults")
        if axes is None:
            fig = plt.figure()
            axes = fig.add_subplot(1, 1, 1)
    
        max_slip = 0.
        min_slip = 0.
        for subfault in self.subfaults:
            slip = subfault.slip
            max_slip = max(abs(slip), max_slip)
            min_slip = min(abs(slip), min_slip)
        if verbose:
            print("Max slip, Min slip: ",max_slip, min_slip)
    
        if slip_color:
            if cmap_slip is None:
                cmap_slip = matplotlib.cm.jet
                #white_purple = colormaps.make_colormap({0.:'w', 1.:[.6,0.2,.6]})
                #cmap_slip = white_purple
            if cmax_slip is None:
                cmax_slip = max_slip
            if cmin_slip is None:
                cmin_slip = 0.
            
        y_ave = 0.
        for subfault in self.subfaults:

            x_top = subfault.centers[0][0]
            y_top = subfault.centers[0][1]
            x_centroid = subfault.centers[1][0]
            y_centroid = subfault.centers[1][1]
            x_corners = [subfault.corners[2][0],
                         subfault.corners[3][0],
                         subfault.corners[0][0],
                         subfault.corners[1][0],
                         subfault.corners[2][0]]
            y_corners = [subfault.corners[2][1],
                         subfault.corners[3][1],
                         subfault.corners[0][1],
                         subfault.corners[1][1],
                         subfault.corners[2][1]]
    
            y_ave += y_centroid
    
    
            # Plot projection of planes to x-y surface:
            if plot_centerline:
                axes.plot([x_top],[y_top],'bo',label="Top center")
                axes.plot([x_centroid],[y_centroid],'ro',label="Centroid")
                axes.plot([x_top,x_centroid],[y_top,y_centroid],'r-')
            if plot_rake:
                tau = (subfault.rake - 90) * numpy.pi/180.
                axes.plot([x_centroid],[y_centroid],'go',markersize=5,label="Centroid")
                dxr = x_top - x_centroid
                dyr = y_top - y_centroid
                x_rake = x_centroid + numpy.cos(tau)*dxr - numpy.sin(tau)*dyr
                y_rake = y_centroid + numpy.sin(tau)*dxr + numpy.cos(tau)*dyr
                axes.plot([x_rake,x_centroid],[y_rake,y_centroid],'g-',linewidth=1)
            if slip_color:
                if slip_time is not None:
                    slip = subfault.dynamic_slip(slip_time)
                else:
                    slip = subfault.slip
                s = min(1, max(0, (slip-cmin_slip)/(cmax_slip-cmin_slip)))
                c = cmap_slip(s*.99)  # since 1 does not map properly with jet
                axes.fill(x_corners,y_corners,color=c,edgecolor='none')
            if plot_box:
                axes.plot(x_corners, y_corners, 'k-')
    
        slipax = axes
            
        y_ave = y_ave / len(self.subfaults)
        slipax.set_aspect(1./numpy.cos(y_ave*numpy.pi/180.))

        if xylim is not None:
            slipax.set_xlim(xylim[:2])
            slipax.set_ylim(xylim[2:])
        if slip_color:
            if slip_time is None:
                slipax.set_title('Slip on fault')
            else:
                slipax.set_title('Slip on fault at time %6.1fs' % slip_time)
        else:
            slipax.set_title('Fault planes')

        slipax.ticklabel_format(format='plain', useOffset=False)
        for label in slipax.get_xticklabels():
            label.set_rotation(20)
        

        if slip_color and (colorbar_shrink > 0):
            cax,kw = matplotlib.colorbar.make_axes(slipax, 
                     shrink=colorbar_shrink)
            norm = matplotlib.colors.Normalize(vmin=cmin_slip,vmax=cmax_slip)
            cb1 = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap_slip, norm=norm)
            cb1.set_label("Slip (m)")
        plt.sca(slipax) # reset the current axis to the main figure

        return slipax
    


    def plot_subfaults_depth(self, axes=None):
        """
        Plot the depth of each subfault vs. x and vs. y in a second plot.
        """
    
        import matplotlib.pyplot as plt

        if axes is None:
            fig, axes = plt.subplots(nrows=2, ncols=1)
        else:
            if len(axes) != 2:
                raise ValueError("The *axes* argument should be a list of ",
                                 "axes objects of length == 2.")
    
        for subfault in self.subfaults:

            x_top = subfault.centers[0][0]
            y_top = subfault.centers[0][1]
            depth_top = subfault.centers[0][2]
            x_bottom = subfault.centers[2][0]
            y_bottom = subfault.centers[2][1]
            depth_bottom = subfault.centers[2][2]
    
            # Plot planes in x-z and y-z to see depths:
            axes[0].plot([x_top, x_bottom], [-depth_top, -depth_bottom])
            axes[1].plot([y_top, y_bottom], [-depth_top, -depth_bottom])
    
        axes[0].set_title('depth vs. x')
        axes[1].set_title('depth vs. y')

        return axes
    

    def containing_rect(self):
        r"""Find containing rectangle of fault in x-y plane.

        Returns tuple of x-limits and y-limits.

        """

        extent = [numpy.infty, -numpy.infty, numpy.infty, -numpy.infty]
        for subfault in self.subfaults:
            for corner in subfault.corners:
                extent[0] = min(corner[0], extent[0])
                extent[1] = max(corner[0], extent[1])
                extent[2] = min(corner[1], extent[2])
                extent[3] = max(corner[1], extent[3])

        return extent

    
    def create_dtopo_xy(self, rect=None, dx=1/60., buffer_size=0.5):
        r"""Create coordinate arrays containing fault with a buffer.

        :Input:
        
         - *rect* - if None, use self.containing_rect
            Otherwise a list [x1,x2,y1,y2]
         - *dx* (int) - Spatial resolution. Defaults to 1" resolution.
         - *buffer_size* (float) - Buffer distance around edge of fault in
           degrees, defaults to 0.5 degrees.

        :Output:
        
         - *x,y* 1-dimensional arrays that cover the desired rect.
           They start at (x1,y1) and may go a bit beyond (x2,y2) depending on dx

        """

        if rect is None:
            rect = self.containing_rect()
        
        rect[0] -= buffer_size
        rect[1] += buffer_size
        rect[2] -= buffer_size
        rect[3] += buffer_size

        mx = int(numpy.ceil(rect[1] - rect[0]) / dx) + 1
        x1 = rect[0]
        x2 = x1 + (mx-1)*dx
        my = int(numpy.ceil(rect[3] - rect[2]) / dx) + 1
        y1 = rect[2]
        y2 = y1 + (my-1)*dx   # note dy==dx

        x = numpy.linspace(x1,x2,mx)
        y = numpy.linspace(y1,y2,my)

        return x,y


    def set_dynamic_slip(self, t):
        r"""
        Set *slip_at_dynamic_t* attribute of all subfaults to slip at the
        requested time *t*.

        :Input:
         - *t* (float) -

        Raises a ValueError exception if this object's rupture_type attribute
        is set to static.
        """

        if self.rupture_type is 'static':
            raise ValueError("Rupture type is set to static.")

        self.dynamic_t = t
        for subfault in self.subfaults:
            subfault.slip_at_dynamic_t = subfault.dynamic_slip(t)




# ==============================================================================
#  Sub-Fault Class
# ==============================================================================
class SubFault(object):
    r"""Basic sub-fault specification.

    Locate fault plane in 3D space
    Note that the coodinate specification is in reference to the fault 
        

    :Coordinates of Fault Plane:

    The attributes *centers* and *corners* are described by the figure below.

    *centers[0,1,2]* refer to the points labeled 0,1,2 below.
      In particular the centroid is given by *centers[1]*.
      Each will be a tuple *(x, y, depth)*.

    *corners[0,1,2,3]* refer to the points labeled a,b,c,d resp. below.
      Each will be a tuple *(x, y, depth)*.
    

    Top edge    Bottom edge
      a ----------- b          ^ 
      |             |          |         ^
      |             |          |         |
      |             |          |         | along-strike direction
      |             |          |         |
      0------1------2          | length  |
      |             |          |
      |             |          |
      |             |          |
      |             |          |
      d ----------- c          v
      <------------->
           width

      <-- up dip direction



    """

    @property
    def corners(self):
        r"""Coordinates of the corners of the fault plane."""
        if self._corners is None:
            self.calculate_geometry()
        return self._corners

    @property
    def centers(self):
        r"""Coordinates along the center-line of the fault plane."""
        if self._centers is None:
            self.calculate_geometry()
        return self._centers

    def __init__(self):
        r"""SubFault initialization routine.
        
        See :class:`SubFault` for more info.

        """
        
        super(SubFault, self).__init__()

        self.strike = None
        r"""Strike direction of subfault in degrees."""
        self.length = None
        r"""Length of subfault in meters."""
        self.width = None
        r"""Width of subfault in meters."""
        self.depth = None
        r"""Depth of subfault based on *coordinate_specification* in meters."""
        self.slip = None
        r"""Slip on subfault in strike direction in meters."""
        self.rake = None
        r"""Rake of subfault movement in degrees."""
        self.dip = None
        r"""Subfault's angle of dip"""
        self.latitude = None
        r"""Latitutde of the subfault based on *coordinate_specification*."""
        self.longitude = None
        r"""Longitude of the subfault based on *coordinate_specification*."""
        self.coordinate_specification = None
        r"""Specifies where the latitude, longitude and depth are measured from."""

        # default value for rigidity = shear modulus
        # Note that standard units for mu is now Pascals.
        # Multiply by 10 to get dyne/cm^2 value.
        self.mu = 4e10
        r"""Rigidity (== shear modulus) in Pascals."""

        self._centers = None
        self._corners = None


    def convert_to_standard_units(self, input_units, verbose=False):
        r"""
        Convert parameters from the units used for input into the standard
        units used in this module.
        """
        params = list(input_units.keys())
        for param in params:
            value = getattr(self, param)
            converted_value = convert_units(value, input_units[param], 1)
            setattr(self,param,converted_value)
            if verbose:
                print("%s %s %s converted to %s %s" \
                    % (param, value, input_units[param], converted_value, \
                       standard_units[param]))


    def Mo(self):
        r"""Calculate the seismic moment for a single subfault

        Returns in units of N-m and assumes mu is in Pascals. 
        """

        total_slip = self.length * self.width * abs(self.slip)
        Mo = self.mu * total_slip
        return Mo


    def __str__(self):
        output = "Subfault Characteristics:\n"
        output += "  Coordinates: (%s, %s) (%s)\n" % (self.longitude,  
                                                      self.latitude, 
                                                  self.coordinate_specification)
        output += "  Dimensions (L,W): (%s, %s) m\n" % (self.length, self.width)
        output += "  Depth: %s m\n" % (self.depth)
        output += "  Rake, Strike, Dip: %s, %s, %s\n" % (self.rake, self.strike,
                                                         self.dip)
        output += "  Slip, Moment: %s m, %s N-m\n" % (self.slip, self.Mo())
        output += "  Fault Centroid: %s\n" % self.centers[1]
        return output


    def calculate_geometry(self):
        r"""Calculate the fault geometry.

        Routine calculates the class attributes *corners* and 
        *centers* which are the corners of the fault plane and 
        points along the centerline respecitvely in 3D space.

        **Note:** *self.coordinate_specification*  specifies the location on each
        subfault that corresponds to the (longitude,latitude) and depth 
        of the subfault.
        Currently must be one of these strings:

          - "bottom center": (longitude,latitude) and depth at bottom center
          - "top center": (longitude,latitude) and depth at top center
          - "centroid": (longitude,latitude) and depth at centroid of plane
          - "noaa sift": (longitude,latitude) at bottom center, depth at top,  
                         This mixed convention is used by the NOAA SIFT
                         database and "unit sources", see:
                         http://nctr.pmel.noaa.gov/propagation-database.html

        The Okada model is expressed assuming (longitude,latitude) and depth
        are at the bottom center of the fault plane, so values must be
        shifted or other specifications.
        """

        # Simple conversion factor of latitude to meters
        lat2meter = util.dist_latlong2meters(0.0, 1.0)[1]

        # Setup coordinate arrays
        self._corners = [[None, None, None], # a 
                                     [None, None, None], # b
                                     [None, None, None], # c
                                     [None, None, None]] # d
        self._centers = [[None, None, None], # 1
                                     [None, None, None], # 2 
                                     [None, None, None]] # 3

        # Set depths
        if self.coordinate_specification == 'top center' or \
           self.coordinate_specification == 'noaa sift':

            self._centers[0][2] = self.depth
            self._centers[1][2] = self.depth            \
                              + 0.5 * self.width * numpy.sin(self.dip * DEG2RAD)
            self._centers[2][2] = self.depth            \
                                    + self.width * numpy.sin(self.dip * DEG2RAD)

        elif self.coordinate_specification == 'centroid':
            self._centers[0][2] = self.depth            \
                              - 0.5 * self.width * numpy.sin(self.dip * DEG2RAD)
            self._centers[1][2] = self.depth
            self._centers[2][2] = self.depth            \
                              + 0.5 * self.width * numpy.sin(self.dip * DEG2RAD)

        elif self.coordinate_specification == 'bottom center':
            self._centers[0][2] = self.depth            \
                              - self.width * numpy.sin(self.dip * DEG2RAD)
            self._centers[1][2] = self.depth            \
                              - 0.5 * self.width * numpy.sin(self.dip * DEG2RAD)
            self._centers[2][2] = self.depth
            
        else:
            raise ValueError("Invalid coordinate specification %s." \
                                                % self.coordinate_specification)

        self._corners[0][2] = self._centers[0][2]
        self._corners[3][2] = self._centers[0][2]
        self._corners[1][2] = self._centers[2][2]
        self._corners[2][2] = self._centers[2][2]
        
        # Locate fault plane in 3D space:  
        # See the class docstring for a guide to labeling of corners/centers.

        # Vector *up_dip* goes from bottom edge to top edge, in meters,
        # from point 2 to point 0 in the figure in the class docstring.


        up_dip = (-self.width * numpy.cos(self.dip * DEG2RAD)           \
                              * numpy.cos(self.strike * DEG2RAD)        \
                              / (LAT2METER * numpy.cos(self.latitude * DEG2RAD)),
                   self.width * numpy.cos(self.dip * DEG2RAD)           \
                              * numpy.sin(self.strike * DEG2RAD) / LAT2METER)
        if self.coordinate_specification == 'top center':

            self._centers[0][:2] = (self.longitude, self.latitude)
            self._centers[1][:2] = (self.longitude - 0.5 * up_dip[0],
                                                self.latitude - 0.5 * up_dip[1])
            self._centers[2][:2] = (self.longitude - up_dip[0],
                                                self.latitude - up_dip[1])

        elif self.coordinate_specification == 'centroid':

            self._centers[0][:2] = (self.longitude + 0.5 * up_dip[0],
                                                self.latitude + 0.5 * up_dip[1])
            self._centers[1][:2] = (self.longitude, self.latitude)
            self._centers[2][:2] = (self.longitude - 0.5 * up_dip[0],
                                                self.latitude  - 0.5 * up_dip[1])

        elif self.coordinate_specification == 'bottom center' or \
             self.coordinate_specification == 'noaa sift':

            # Non-rotated lcoations of center-line coordinates
            self._centers[0][:2] = (self.longitude + up_dip[0],
                                                self.latitude + up_dip[1])
            self._centers[1][:2] = (self.longitude + 0.5 * up_dip[0],
                                                self.latitude + 0.5 * up_dip[1])
            self._centers[2][:2] = (self.longitude, self.latitude)

        else:
            raise ValueError("Unknown coordinate specification '%s'."       \
                                                % self.coordinate_specification)

        # Calculate coordinates of corners:
        # Vector *strike* goes along the top edge from point 1 to point a
        # in the figure in the class docstring.

        up_strike = (0.5 * self.length * numpy.sin(self.strike * DEG2RAD) \
           / (lat2meter * numpy.cos(self._centers[2][1] * DEG2RAD)),
                     0.5 * self.length * numpy.cos(self.strike * DEG2RAD) \
           / lat2meter)
        
        self._corners[0][:2] = (self._centers[0][0] 
                                                                 + up_strike[0],
                                            self._centers[0][1] 
                                                                 + up_strike[1])
        self._corners[1][:2] = (self._centers[2][0] 
                                                                 + up_strike[0], 
                                            self._centers[2][1] 
                                                                 + up_strike[1])
        self._corners[2][:2] = (self._centers[2][0] 
                                                                 - up_strike[0],
                                            self._centers[2][1] 
                                                                 - up_strike[1])
        self._corners[3][:2] = (self._centers[0][0] 
                                                                 - up_strike[0],
                                            self._centers[0][1] 
                                                                 - up_strike[1])

    
    def okada(self, x, y):
        r"""
        Apply Okada to this subfault and return a DTopography object.

        :Input:
          - x,y are 1d arrays
        :Output:
          - DTopography object with dZ array of shape (1,len(x),len(y))
                with single static displacement and times = [0.].

        Currently only calculates the vertical displacement.

        Okada model is a mapping from several fault parameters
        to a surface deformation.
        See Okada 1985 [Okada85]_, or Okada 1992, Bull. Seism. Soc. Am.
        
        okadamap function riginally written in Python by Dave George for
        Clawpack 4.6 okada.py routine, with some routines adapted
        from fortran routines written by Xiaoming Wang.

        Rewritten and made more flexible by Randy LeVeque

        **Note:** *self.coordinate_specification* (str) specifies the location on 
        each subfault that corresponds to the (longitude,latitude) and depth 
        subfault.

        See the documentation for *SubFault.calculate_geometry* for dicussion of the 
        possible values *self.coordinate_specification* can take.

        """

        # Okada model assumes x,y are at bottom center:
        x_bottom = self.centers[2][0]
        y_bottom = self.centers[2][1]
        depth_bottom = self.centers[2][2]

        length = self.length
        width = self.width
        depth = self.depth
        slip = self.slip

        halfL = 0.5*length
        w  =  width

        # convert angles to radians:
        ang_dip = DEG2RAD * self.dip
        ang_rake = DEG2RAD * self.rake
        ang_strike = DEG2RAD * self.strike
    
        X,Y = numpy.meshgrid(x, y)   # use convention of upper case for 2d
    
        # Convert distance from (X,Y) to (x_bottom,y_bottom) from degrees to
        # meters:
        xx = LAT2METER * numpy.cos(DEG2RAD * Y) * (X - x_bottom)   
        yy = LAT2METER * (Y - y_bottom)
    
    
        # Convert to distance along strike (x1) and dip (x2):
        x1 = xx * numpy.sin(ang_strike) + yy * numpy.cos(ang_strike) 
        x2 = xx * numpy.cos(ang_strike) - yy * numpy.sin(ang_strike) 
    
        # In Okada's paper, x2 is distance up the fault plane, not down dip:
        x2 = -x2
    
        p = x2 * numpy.cos(ang_dip) + depth_bottom * numpy.sin(ang_dip)
        q = x2 * numpy.sin(ang_dip) - depth_bottom * numpy.cos(ang_dip)
    
        f1 = self._strike_slip(x1 + halfL, p,     ang_dip, q)
        f2 = self._strike_slip(x1 + halfL, p - w, ang_dip, q)
        f3 = self._strike_slip(x1 - halfL, p,     ang_dip, q)
        f4 = self._strike_slip(x1 - halfL, p - w, ang_dip, q)
    
        g1=self._dip_slip(x1 + halfL, p,     ang_dip, q)
        g2=self._dip_slip(x1 + halfL, p - w, ang_dip, q)
        g3=self._dip_slip(x1 - halfL, p,     ang_dip, q)
        g4=self._dip_slip(x1 - halfL, p - w, ang_dip, q)
    
        # Displacement in direction of strike and dip:
        ds = slip * numpy.cos(ang_rake)
        dd = slip * numpy.sin(ang_rake)
    
        us = (f1 - f2 - f3 + f4) * ds
        ud = (g1 - g2 - g3 + g4) * dd
    
        dz = (us+ud)

        dtopo = DTopography()
        dtopo.X = X
        dtopo.Y = Y
        dtopo.dZ = numpy.array(dz, ndmin=3)
        dtopo.times = [0.]
        self.dtopo = dtopo
        return dtopo

    # Utility functions for okada:

    def _strike_slip(self, y1, y2, ang_dip, q):
        """
        Used for Okada's model
        Methods from Yoshimitsu Okada (1985)
        """
        sn = numpy.sin(ang_dip)
        cs = numpy.cos(ang_dip)
        d_bar = y2*sn - q*cs
        r = numpy.sqrt(y1**2 + y2**2 + q**2)
        xx = numpy.sqrt(y1**2 + q**2)
        a4 = 2.0*poisson/cs*(numpy.log(r+d_bar) - sn*numpy.log(r+y2))
        f = -(d_bar*q/r/(r+y2) + q*sn/(r+y2) + a4*sn)/(2.0*3.14159)
    
        return f
    
    
    def _dip_slip(self, y1, y2, ang_dip, q):
        """
        Based on Okada's paper (1985)
        Added by Xiaoming Wang
        """
        sn = numpy.sin(ang_dip)
        cs = numpy.cos(ang_dip)
    
        d_bar = y2*sn - q*cs;
        r = numpy.sqrt(y1**2 + y2**2 + q**2)
        xx = numpy.sqrt(y1**2 + q**2)
        a5 = 4.*poisson/cs*numpy.arctan((y2*(xx+q*cs)+xx*(r+xx)*sn)/y1/(r+xx)/cs)
        f = -(d_bar*q/r/(r+y1) + sn*numpy.arctan(y1*y2/q/r) - a5*sn*cs)/(2.0*3.14159)
    
        return f

    
    def dynamic_slip(self, t):
        r"""
        For a dynamic fault, compute the slip at time t.
        Assumes the following attributes are set:

          - *rupture_time*
          - *rise_time*
          - *rise_time_ending*: optional, defaults to *rise_time*

        """
        if (self.rupture_time is None) or (self.rise_time is None):
            raise ValueError("Computing a dynamic slip only works for dynamic ",
                             "rupture types")

        t0 = self.rupture_time
        t1 = self.rise_time
        t2 = getattr(self,'rise_time_ending',None)
        rf = rise_fraction(t,t0,t1,t2)
        return rf * self.slip
            



# ==============================================================================
#  UCSB sub-class of Fault
# ==============================================================================
class UCSBFault(Fault):

    r"""Fault subclass for reading in subfault format models from UCSB

    Read in subfault format models produced by Chen Ji's group at UCSB,
    downloadable from:  

        http://www.geol.ucsb.edu/faculty/ji/big_earthquakes/home.html

    """

    def __init__(self, path=None, **kwargs):
        r"""UCSBFault initialization routine.
        
        See :class:`UCSBFault` for more info.

        """

        self.num_cells = [None, None]   # RJL: Why needed??
                                        # Do not really but the specification
                                        # has it.

        super(UCSBFault, self).__init__()

        if path is not None:
            self.read(path, **kwargs)


    def read(self, path, rupture_type='static'):
        r"""Read in subfault specification at *path*.

        Creates a list of subfaults from the subfault specification file at
        *path*.

        Subfault format contains info for dynamic rupture, so can specify 
        rupture_type = 'static' or 'dynamic'

        """

        self.rupture_type = rupture_type

        # Read header of file
        regexp_dx = re.compile(r"Dx=[ ]*(?P<dx>[^k]*)")
        regexp_dy = re.compile(r"Dy=[ ]*(?P<dy>[^k]*)")
        regexp_nx = re.compile(r"nx[^=]*=[ ]*(?P<nx>[^D]*)")
        regexp_ny = re.compile(r"ny[^=]*=[ ]*(?P<ny>[^D]*)")
        found_subfault_discretization = False
        found_subfault_boundary = False
        header_lines = 0
        with open(path, 'r') as subfault_file:
            # Find fault secgment discretization
            for (n,line) in enumerate(subfault_file):
                result_dx = regexp_dx.search(line)
                result_dy = regexp_dy.search(line)
                result_nx = regexp_nx.search(line)
                result_ny = regexp_ny.search(line)

                if result_dx and result_dy:
                    dx = float(result_dx.group('dx'))
                    dy = float(result_dy.group('dy'))
                    self.num_cells[0] = int(result_nx.group('nx'))
                    self.num_cells[1] = int(result_ny.group('ny'))
                    found_subfault_discretization = True
                    break
            header_lines += n

            # Parse boundary
            in_boundary_block = False
            boundary_data = []
            for (n,line) in enumerate(subfault_file):
                if line[0].strip() == "#":
                    if in_boundary_block and len(boundary_data) == 5:
                        found_subfault_boundary = True
                        break
                else:
                    in_boundary_block = True
                    boundary_data.append([float(value) for value in line.split()])

        # Assume that there is a column label right underneath the boundary
        # specification
        header_lines += n + 2

        # Check to make sure last boundary point matches, then throw away
        if boundary_data[0] != boundary_data[4]:
            raise ValueError("Boundary specified incomplete: ",
                             "%s" % boundary_data)

        # Locate fault plane in 3D space - see SubFault `calculate_geometry`
        # for
        # These are only useful in sub-faults
        # a schematic of where these points are
        # self._corners = [None, # a 
        #                              None, # b
        #                              None, # c
        #                              None] # d
        # self._centers = [[0.0, 0.0, 0.0], # 1
        #                              [0.0, 0.0, 0.0], # 2 
        #                              [0.0, 0.0, 0.0]] # 3
        # # :TODO: Is the order of this a good assumption?
        # self._corners[0] = boundary_data[0]
        # self._corners[3] = boundary_data[1]
        # self._corners[2] = boundary_data[2]
        # self._corners[1] = boundary_data[3]
        
        # # Calculate center by averaging position of appropriate corners
        # for (n, corner) in enumerate(self._corners):
        #     for i in xrange(3):
        #         self._centers[1][i] += corner[i] / 4
        #     if n == 0 or n == 4:
        #         for i in xrange(3):
        #             self._centers[0][i] += corner[i] / 2
        #     else:
        #         for i in xrange(3):
        #             self._centers[2][i] += corner[i] / 2


        if not (found_subfault_boundary and found_subfault_discretization):
            raise ValueError("Could not find base fault characteristics in ",
                             "subfault specification file at %s." % path)

        # Calculate center of fault

        column_map = {"latitude":0, "longitude":1, "depth":2, "slip":3,
                       "rake":4, "strike":5, "dip":6, "rupture_time":7,
                       "rise_time":8, "rise_time_ending":9, "mu":10}
        defaults = {"length":dx, "width":dy}
        input_units = {"slip":"cm", "depth":"km", 'mu':"dyne/cm^2",
                                   "length":"km", "width":"km"}

        super(UCSBFault, self).read(path, column_map, skiprows=header_lines,
                                coordinate_specification="centroid",
                                input_units=input_units, defaults=defaults)



# ==============================================================================
#  CSV sub-class of Fault
# ==============================================================================
class CSVFault(Fault):

    r"""Fault subclass for reading in CSV formatted files

    Assumes that the first row gives the column headings
    """

    def read(self, path, input_units={}, coordinate_specification="top center",
                         rupture_type="static", verbose=False):
        r"""Read in subfault specification at *path*.

        Creates a list of subfaults from the subfault specification file at
        *path*.

        """

        possible_column_names = """longitude latitude length width depth strike dip
                          rake slip mu rupture_time rise_time rise_time_ending""".split()
        param = {}
        for n in possible_column_names:
            param[n] = n
        # alternative names that might appear in csv file:
        param["rigidity"] = "mu"
        param["rupture time"] = "rupture_time"
        param["rise time"] = "rise_time"

        # Read header of file
        with open(path, 'r') as subfault_file:
            header_line = subfault_file.readline().split(",")
            column_map = {}
            for (n,column_heading) in enumerate(header_line):
                if "(" in column_heading:
                    # Strip out units if present
                    unit_start = column_heading.find("(")
                    unit_end = column_heading.find(")")
                    column_name = column_heading[:unit_start].lower()
                    units = column_heading[unit_start+1:unit_end]
                    if verbose and input_units.get(column_name,units) != units:
                        print("*** Warning: input_units[%s] reset to %s" \
                              % (column_name, units))
                        print("    based on file header")
                        input_units[column_name] = units

                else:
                    column_name = column_heading.lower()
                column_name = column_name.strip()
                if column_name in list(param.keys()):
                    column_key = param[column_name]
                    column_map[column_key] = n
                else:
                    print("*** Warning: column name not recognized: %s" \
                        % column_name)

        super(CSVFault, self).read(path, column_map=column_map, skiprows=1,
                                delimiter=",", input_units=input_units,
                                coordinate_specification=coordinate_specification,
                                rupture_type=rupture_type)



# ==============================================================================
#  Sift sub-class of Fault
# ==============================================================================
class SiftFault(Fault):

    r"""
    Define a fault by specifying the slip on a subset of the SIFT unit sources.
    The database is read in by load_sift_unit_sources.
    See http://www.pmel.noaa.gov/pubs/PDF/gica2937/gica2937.pdf
    for a discussion of these unit sources, although the database used
    is more recent than what is reported in that paper and uses different
    notation for the subfault names.
    The subfault database used was downloaded from

        http://sift.pmel.noaa.gov/ComMIT/compressed/info_sz.dat

    Example:

        >>> sift_slip = {'acsza1':2, 'acszb1':3}
        >>> fault = SiftFault(sift_slip)

    results in a fault with two specified subfaults with slip of 2 and 3 meters.
    """

    def __init__(self, sift_slip=None):
        r"""SiftFault initialization routine.
        
        See :class:`SiftFault` for more info.

        """
        
        super(SiftFault, self).__init__()
        self._load_sift_unit_sources()
        if sift_slip is not None:
            self.set_subfaults(sift_slip)
            

    def set_subfaults(self,sift_slip):
        r"""
        *sift_slip* (dict) is a dictionary with key = name of unit source
                    and value = magnitude of slip to assign (in meters).
        """
        self.subfaults = []
        for k,v in six.iteritems(sift_slip):
            subfault = self.sift_subfaults[k]
            subfault.slip = v
            self.subfaults.append(subfault)


    def _load_sift_unit_sources(self):
        r"""
        Load SIFT unit source subfault data base. 
        File was downloaded from
            http://sift.pmel.noaa.gov/ComMIT/compressed/info_sz.dat
        """

        unit_source_file = os.path.join(os.path.dirname(__file__), 'data', 
                                        'info_sz.dat.txt')
        self.input_units = {'length':'km', 'width':'km', 'depth':'km', 'slip':'m',
                 'mu':"dyne/cm^2"}

        self.sift_subfaults = {}

        with open(unit_source_file, 'r') as sift_file:
            # Skip first two lines
            sift_file.readline(); sift_file.readline()
            for line in sift_file:
                tokens = line.split(',')
                name = tokens[0]
                # url = tokens[1]
                subfault = SubFault()
                subfault.longitude = float(tokens[2])
                subfault.latitude = float(tokens[3])
                subfault.slip = float(tokens[4])
                subfault.strike = float(tokens[5])
                subfault.dip = float(tokens[6])
                subfault.depth = float(tokens[7])
                subfault.length = float(tokens[8])
                subfault.width = float(tokens[9])
                subfault.rake = float(tokens[10])
                subfault.coordinate_specification = "noaa sift"
                # subfault.mu = ??  ## currently using SubFault default
                subfault.convert_to_standard_units(self.input_units)
                self.sift_subfaults[name] = subfault


# ==============================================================================
#  Subdivided plane sub-class of Fault
# ==============================================================================
class SubdividedPlaneFault(Fault):

    r"""
    Define a fault by starting with a single fault plane (specified as 
    *base_subfault* of class *SubFault*) and subdividing the fault plane
    into a rectangular array of *nstrike* by *ndip* equally sized subfaults.

    By default,the slip on each subfault will be initialized to
    *base_subfault.slip* so that the slip is uniform over the original plane
    and the seismic moment is independent of the number of subdivisions.  

    Alternatively, the slip distribution can be specified by providing a
    function *slip_distribution*, which should be a function of *(xi,eta)*
    with each variable ranging from 0 to 1.  *xi* varies from 0 at the top
    of the fault to 1 at the bottom in the down-dip direction. 
    *eta* varies from one edge of the fault to the other moving in the
    strike direction.  This function will be evaluated at the centroid of
    each subfault to set the slip.

    Can also specify a desired seismic moment Mo in which case the slips will
    be rescaled at the end so the total seismic moment is Mo.  In this case
    the *slip_distribution* function only indicates the relative slip
    between subfaults.

    """


    def __init__(self, base_subfault, nstrike=1, ndip=1,
                       slip_function=None, Mo=None):
        r"""SubdivdedPlaneFault initialization routine.
        
        See :class:`SubdivdedPlaneFault` for more info.

        """

        super(SubdividedPlaneFault, self).__init__()

        self.base_subfault = base_subfault
        self.nstrike = nstrike
        self.ndip = ndip

        self.subdivide(nstrike, ndip, slip_function, Mo)


    def subdivide(self, nstrike=1, ndip=1, slip_function=None, Mo=None):
        r"""Subdivide the fault plane into nstrike * ndip subfaults."""

        # may have changed resolution:
        self.nstrike = nstrike
        self.ndip = ndip

        base_subfault = self.base_subfault
        strike = base_subfault.strike
        dip = base_subfault.dip
        rake = base_subfault.rake
        slip = base_subfault.slip
        length = base_subfault.length
        width = base_subfault.width

        # unpack corners from fault plane geometry:
        x_corners = [base_subfault.corners[2][0],
                     base_subfault.corners[3][0],
                     base_subfault.corners[0][0],
                     base_subfault.corners[1][0],
                     base_subfault.corners[2][0]]
        y_corners = [base_subfault.corners[2][1],
                     base_subfault.corners[3][1],
                     base_subfault.corners[0][1],
                     base_subfault.corners[1][1],
                     base_subfault.corners[2][0]]

        # set depth at corners:
        depth_top = base_subfault.centers[0][2]
        depth_bottom = base_subfault.centers[2][2]
        d_corners = [depth_bottom, depth_top, depth_top, depth_bottom,
                     depth_bottom]

        # coefficients for bilinear interpolants:
        cx = [x_corners[1], 
              x_corners[0] - x_corners[1],
              x_corners[2] - x_corners[1], 
              x_corners[3] + x_corners[1] - x_corners[2] - x_corners[0]]
        cy = [y_corners[1], 
              y_corners[0] - y_corners[1],
              y_corners[2] - y_corners[1], 
              y_corners[3] + y_corners[1] - y_corners[2] - y_corners[0]]
        cd = [d_corners[1], 
              d_corners[0] - d_corners[1],
              d_corners[2] - d_corners[1], 
              d_corners[3] + d_corners[1] - d_corners[2] - d_corners[0]]

        self.subfaults = []

        # determine coordinates for each subfault.
        # note that xi goes from 0 to 1 from top to bottom in dip direction,
        #          eta goes from 0 to 1 in along-strike direction.
        dxi = 1. / ndip
        deta = 1. / nstrike
        for i in range(ndip):
            xi = numpy.array([i, i+0.5, i+1.]) * dxi # xi at top, center, bottom
            for j in range(nstrike):
                eta = (j+0.5)*deta
                # interpolate longitude,latitude,depth from corners:
                x_sf = cx[0] + cx[1]*xi + cx[2]*eta + cx[3]*xi*eta
                y_sf = cy[0] + cy[1]*xi + cy[2]*eta + cy[3]*xi*eta
                d_sf = cd[0] + cd[1]*xi + cd[2]*eta + cd[3]*xi*eta

                subfault = SubFault()

                if base_subfault.coordinate_specification == 'centroid':
                    subfault.longitude = x_sf[1]
                    subfault.latitude = y_sf[1]
                    subfault.depth = d_sf[1]
                elif base_subfault.coordinate_specification == 'top center':
                    subfault.longitude = x_sf[0]
                    subfault.latitude = y_sf[0]
                    subfault.depth = d_sf[0]
                elif base_subfault.coordinate_specification == 'noaa sift':
                    subfault.longitude = x_sf[2]
                    subfault.latitude = y_sf[2]
                    subfault.depth = d_sf[0]
                else:   
                    msg = "Unrecognized coordinate_specification: %s" \
                            % base_subfault.coordinate_specification
                    raise NotImplementedError(msg)

                subfault.dip = dip
                subfault.strike = strike
                subfault.rake = rake
                subfault.length = length / nstrike
                subfault.width = width / ndip
                subfault.slip = slip
                subfault.coordinate_specification = \
                        base_subfault.coordinate_specification
                subfault.mu = base_subfault.mu

                self.subfaults.append(subfault)

        if slip_function is not None:
            self.set_slip(nstrike, ndip, slip_function, Mo)

    def set_slip(self, nstrike, ndip, slip_function, Mo=None):

        self.slip_function = slip_function
        dxi = 1. / ndip
        deta = 1. / nstrike
        Mo_0 = 0.
        k = 0
        for i in range(ndip):
            xi = (i+0.5) * dxi
            for j in range(nstrike):
                eta = (j+0.5) * deta
                subfault = self.subfaults[k]  
                k = k+1
                subfault.slip = slip_function(xi,eta)
                Mo_0 += subfault.Mo()

        if Mo is not None:
            # rescale slip on each subfault to achieve desired seismic moment
            Mo_ratio = Mo / Mo_0
            for k in range(len(self.subfaults)):
                self.subfaults[k].slip *= Mo_ratio



# ==============================================================================
#  Tensor product sub-class of Fault
# ==============================================================================
class TensorProductFault(SubdividedPlaneFault):

    r"""
    Define a fault by starting with a single fault plane (specified as 
    *fault_plane* of class *SubFault*) and subdividing the fault plane
    into a rectangular array of *nstrike* by *ndip* equally sized subfaults.
    
    Then define the slip on each subfault via
    two one-dimensional functions *slip_along_strike* and
    *slip_down_dip* that specify the slip as a function of fractional
    distance in the along-strike and down-dip direction respectively
    (i.e. the argument of each goes from 0 to 1).

    Setting either to None defaults to constant function 1.

    The slip is set by evaluating the tensor product at the centroid of
    each subfault.

    Can specify a desired seismic moment Mo in which case the slips will
    be rescaled at the end.

    """

    def __init__(self, fault_plane, slip_along_strike=None, slip_down_dip=None,
                      nstrike=1, ndip=1, Mo=None):
        r"""TensorProductFault initialization routine.
        
        See :class:`TensorProductFault` for more info.

        """
        
        # perform the subdivision and set parameters on each subfault:
        super(TensorProductFault, self).__init__(fault_plane, nstrike, ndip)

        if slip_along_strike is None:
            # set to constant in the strike direction if not specified
            slip_along_strike = lambda eta: 1.0
        if slip_down_dip is None:
            # set to constant in the dip direction if not specified
            slip_down_dip = lambda xi: 1.0
