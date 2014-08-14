#!/usr/bin/env python
# encoding: utf-8

r"""
GeoClaw dtopotools Module

Module provides several functions for dealing with changes to topography (usually
due to earthquakes) including reading sub-fault specifications, writing out 
dtopo files, and calculating Okada based deformations.

:Classes:

    DTopography
    SubFault
    Fault
    UCSBFault
    CSVFault
    SiftFault
    SegmentedPlaneFault
    
:Functions:


:TODO:
 - List functions and classes in this docstring?
 - This file contains both the older dtopo functionality and a new class, should
   merge as much functionality into the class as possible ensuring nothing is
   left behind.
 - Refactor Okada functionality
"""

import os
import sys
import re

import numpy

import clawpack.geoclaw.topotools as topotools

# ==============================================================================
#  Constants
# ==============================================================================
# Poisson ratio for Okada 
from clawpack.geoclaw.util import DEG2RAD, LAT2METER
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
check = [unit_conversion_factor[standard_units[param]] is 1. for param in \
         standard_units.keys()]
assert numpy.alltrue(check), \
        """Conversion factors should be 1 for all standard_units"""


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
        raise ValueError("Unrecognized io_units %s, must be one of %s" \
              % (io_units, unit_conversion_factor.keys()))
    if direction == 1:
        converted_value = value * factor
    elif direction == 2:
        converted_value = value / factor
    else:
        raise ValueError("Unrecognized direction, must be 1 or 2")

    return converted_value



def plot_dz_contours(x, y, dz, axes=None, dz_interval=0.5, verbose=False,
                               fig_kwargs={}):
    r"""For plotting seafloor deformation dz"""
    import matplotlib.pyplot as plt

    if axes is None:
        fig = plt.figure(**fig_kwargs)
        axes = fig.add_subplot(1, 1, 1)

    dzmax = max(dz.max(), -dz.min()) + dz_interval
    clines1 = numpy.arange(dz_interval, dzmax, dz_interval)
    clines = list(-numpy.flipud(clines1)) + list(clines1)

    # Create axes if needed
    if axes is None:
        fig = plt.figure(**fig_kwargs)
        axes = fig.add_subplot(111)

    if len(clines) > 0:
        if verbose:
            print "Plotting contour lines at: ",clines
        axes.contour(x, y, dz, clines, colors='k')
    else:   
        print "No contours to plot"

    return axes

def plot_dz_colors(x, y, dz, axes=None, cmax_dz=None, dz_interval=None,
                   add_colorbar=True, verbose=False, fig_kwargs={}):
    r"""
    Plot sea floor deformation dz as colormap with contours
    """

    from clawpack.visclaw import colormaps
    import matplotlib.pyplot as plt

    if axes is None:
        fig = plt.figure(**fig_kwargs)
        axes = fig.add_subplot(1, 1, 1)
    #print "+++ in plot_dz_colors, axes = ",axes
    #print "+++ in plot_dz_colors, id(axes) = ",id(axes)

    dzmax = numpy.abs(dz).max()
    if cmax_dz is None:
        if dzmax < 1.e-12:
            cmax_dz = 0.1
        else:
            cmax_dz = dzmax
    cmap = colormaps.blue_white_red
    extent = [x.min(), x.max(), y.min(), y.max()]
    im = axes.imshow(dz, extent=extent, cmap=cmap, origin='lower')
    im.set_clim(-cmax_dz,cmax_dz)
    if add_colorbar:
        cbar = plt.colorbar(im, ax=axes)
        cbar.set_label("Deformation (m)")
    
    if dz_interval is None:
        dz_interval = cmax_dz/10.
    clines1 = numpy.arange(dz_interval, dzmax + dz_interval, dz_interval)
    clines = list(-numpy.flipud(clines1)) + list(clines1)
    if len(clines) > 0:
        if verbose:
            print "Plotting contour lines at: ",clines
        axes.contour(x,y,dz,clines,colors='k',linestyles='solid')
    elif verbose:
        print "No contours to plot"

    y_ave = 0.5*(y.min() + y.max())
    axes.set_aspect(1./numpy.cos(y_ave*numpy.pi/180.))
    axes.ticklabel_format(format='plain', useOffset=False)
    ## RJL: setting labels like this gives None's as labels:
    #axes.set_xticklabels([label.set_rotation(80) 
    #                                       for label in axes.get_xticklabels()])
    axes.set_title('Seafloor deformation')
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
       0 for t <= t0, 
       1 for t >= t0 + t_rise + t_rise_ending 
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

        self.dz_list = []
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

        input
        -----
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
                print "Loaded file %s with %s lines" %(path,data.shape[0])
            t = list(set(data[:,0]))
            t.sort()
            print "times found: ",t
            ntimes = len(t)
            tlast = t[-1]
            lastlines = data[data[:,0]==tlast]
            xvals = list(set(lastlines[:,1]))
            xvals.sort()
            mx = len(xvals)
            my = len(lastlines) / mx
            print "Read dtopo: mx=%s and my=%s, at %s times" % (mx,my,ntimes)
            X = numpy.reshape(lastlines[:,1],(my,mx))
            Y = numpy.reshape(lastlines[:,2],(my,mx))
            Y = numpy.flipud(Y)
            dz_list = []
            print "Returning dZ as a list of mx*my arrays"
            for n in range(ntimes):
                i1 = n*mx*my
                i2 = (n+1)*mx*my
                dz = numpy.reshape(data[i1:i2,3],(my,mx))
                dz = numpy.flipud(dz)
                dz_list.append(dz)
            self.X = X
            self.Y = Y
            self.x = X[0,:]
            self.y = Y[:,0]
            self.times = t
            self.dz_list = dz_list

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
            dz_list = []
            if dtopo_type==3:
                # my lines with mx values on each
                for k,t in enumerate(times):
                    dZk = numpy.reshape(dZvals[k*my:(k+1)*my, :], (my,mx))
                    dZk = numpy.flipud(dZk)
                    dz_list.append(dZk)
            else:
                # dtopo_type==2 ==> mx*my lines with 1 values on each
                for k,t in enumerate(times):
                    dZk = numpy.reshape(dZvals[k*mx*my:(k+1)*mx*my], (my,mx))
                    dZk = numpy.flipud(dZk)
                    dz_list.append(dZk)
                    
            self.x = x
            self.y = y
            self.X, self.Y = numpy.meshgrid(x,y)
            self.times = times
            self.dz_list = dz_list

        else:
            raise ValueError("Only topography types 1, 2, and 3 are supported,",
                             " given %s." % dtopo_type)


    def write(self, path=None, dtopo_type=None):
        r"""Write out subfault resulting dtopo to file at *path*.

        input
        -----
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
        assert abs(dx-dy) <1e-12, \
            "*** dx = %g not equal to dy = %g" % (dx,dy)

        # Construct each interpolating function and evaluate at new grid
        ## Shouldn't need to interpolate in time.
        with open(path, 'w') as data_file:

            if dtopo_type == 0:
                # Topography file with 3 columns, x, y, dz written from the
                # upper left corner of the region
                Y_flipped = numpy.flipud(self.Y)
                dZ_flipped = numpy.flipud(self.dz_list[0])

                for j in xrange(self.Y.shape[0]):
                    for i in xrange(self.X.shape[1]):
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
                    dZ_flipped = numpy.flipud(self.dz_list[n])

                    for j in xrange(self.Y.shape[0]):
                        for i in xrange(self.X.shape[1]):
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
                                                  % tuple(self.dz_list[n][j,:]))
                            data_file.write("\n")

            else:
                raise ValueError("Only topography types 1, 2, and 3 are ",
                                 "supported, given %s." % dtopo_type)


    def dz(self, t):
        """
        Interpolate dz_list to specified time t and return deformation dz.
        """
        from matplotlib.mlab import find
        if t <= self.times[0]:
            return self.dz_list[0]
        elif t >= self.times[-1]:
            return self.dz_list[-1]
        else:
            n = max(find(self.times <= t))
            t1 = self.times[n]
            t2 = self.times[n+1]
            dz = (t2-t)/(t2-t1) * self.dz_list[n] + (t-t1)/(t2-t1) * self.dz_list[n+1]
            return dz

    def dz_max(self):
        r"""Return max(abs(dz)) over all dz in self.dz_list, the maximum
        surface deformation for this dtopo."""

        dzm = 0.
        for dz in self.dz_list:
            dzm = max(dzm, abs(dz).max())
        return dzm

    def plot_dz_colors(self, t, axes=None, cmax_dz=None, dz_interval=None, 
                                fig_kwargs={}):
        """
        Interpolate dz_list to specified time t and then call module function
        plot_dz_colors.
        """
        axes = plot_dz_colors(self.X, self.Y, self.dz(t), axes=axes,
                              cmax_dz=cmax_dz, dz_interval=dz_interval,
                              fig_kwargs=fig_kwargs)
        return axes


    def plot_dz_contours(self, t, dz_interval=0.5, axes=None, fig_kwargs={}):
        """
        Interpolate dz_list to specified time t and then call module function
        plot_dz_contours.
        """
        axes = plot_dz_contours(self.X, self.Y, self.dz(t), 
                                dz_interval=dz_interval)
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

    def __init__(self, subfaults=None, input_units={}):
        r"""Fault initialization routine.
        
        See :class:`Fault` for more info.

        """

        # Parameters for subfault specification
        self.rupture_type = 'static' # 'static' or 'dynamic'
        #self.times = numpy.array([0., 1.])   # or just [0.] ??
        self.dtopo = None

        # Default units of each parameter type
        self.input_units = standard_units
        self.input_units.update(input_units)
        
        if subfaults is not None:
            if not isinstance(subfaults, list):
                raise ValueError("Input parameter subfaults must be a list.")
            self.subfaults = subfaults
            for subfault in self.subfaults:
                subfault.convert_to_standard_units(input_units)


    def read(self, path, column_map, coordinate_specification="centroid",
                                     rupture_type="static", skiprows=0, 
                                     delimiter=None, input_units={}, defaults=None):
        r"""Read in subfault specification at *path*.

        Creates a list of subfaults from the subfault specification file at
        *path*.
        Inputs:
          - *path* (str) file to read in, should contain subfaults, one per line
          - *column_map* (dict) specifies mapping from parameter to the column
            of the input file that contains values for this parameter, e.g.
                column_map = {"latitude":0, "longitude":1, "depth":2, "slip":3,
                               "rake":4, "strike":5, "dip":6}
          - *coordinate_specification* (str) specifies the location on each
            subfault that corresponds to the (longitude,latitude) and depth 
            of the subfault.  See the documentation for *SubFault.set_geometry*.
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
        data = numpy.genfromtxt(path, skiprows=skiprows, delimiter=delimiter)
        if len(data.shape) == 1:
            data = numpy.array([data])

        self.coordinate_specification = coordinate_specification
        self.input_units = standard_units
        self.input_units.update(input_units)
        self.subfaults = []
        for n in xrange(data.shape[0]):

            new_subfault = SubFault()
            new_subfault.coordinate_specification = coordinate_specification
            
            for (var, column) in column_map.iteritems():
                if isinstance(column, tuple) or isinstance(column, list):
                    setattr(new_subfault, var, [None for k in column])
                    for (k, index) in enumerate(column):
                        getattr(new_subfault, var)[k] = data[n, index]
                else:
                    setattr(new_subfault, var, data[n, column])

            if defaults is not None:
                for param in defaults.iterkeys():
                    setattr(new_subfault, param, defaults[param]) 

            new_subfault.input_units = input_units
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

        self.output_units = standard_units
        self.output_units.update(output_units)

        if style is not None:
            msg =  "style option not yet implemented, use column_map"
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
            # write header:
            data_file.write('Subfaults file with coordinate_specification:  ')
            data_file.write('%s, \n' % self.coordinate_specification)
            data_file.write('Units: %s, \n' % str(output_units))
            s = ""
            for param in column_list:
                s = s + delimiter + param.rjust(15)
            data_file.write(s + '\n')
            for subfault in self.subfaults:
                s = ""
                for param in column_list:
                    value = getattr(subfault,param)
                    if output_units.has_key(param):
                        converted_value = convert_units(value, 
                                    self.output_units[param], direction=2)
                    s = s + delimiter + format[param] % value
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
            print "Making Okada dz for each of %s subfaults" \
                  % len(self.subfaults)

        for k,subfault in enumerate(self.subfaults):
            if verbose:
                sys.stdout.write("%s.." % k)
                sys.stdout.flush()
            subfault.okada(x,y)  # sets subfault.dtopo with times=[0]
                                 # and len(subfault.dtopo.dz_list) == 1
        if verbose:
            sys.stdout.write("\nDone\n")

        if self.rupture_type == 'static':
            if len(times) > 2:
                raise ValueError("For static deformation, need len(times) <= 2")
            dz = numpy.zeros(X.shape)
            for subfault in self.subfaults:
                dz += subfault.dtopo.dz_list[0]

            if len(times) == 1:
                dtopo.dz_list = [dz]   # only final deformation stored
            elif len(times) == 2:
                dz0 = numpy.zeros(X.shape)
                dtopo.dz_list = [dz0, dz]

        elif self.rupture_type in ['dynamic','kinematic']:

            t_prev = -1.e99
            dz_list = []
            dz = numpy.zeros(X.shape)
            for t in times:
                for k,subfault in enumerate(self.subfaults):
                    t0 = getattr(subfault,'rupture_time',0)
                    t1 = getattr(subfault,'rise_time',0.5)
                    t2 = getattr(subfault,'rise_time_ending',None)
                    rf = rise_fraction([t_prev,t],t0,t1,t2)
                    dfrac = rf[1] - rf[0]
                    if dfrac > 0.:
                        dz = dz + dfrac * subfault.dtopo.dz_list[0]
                dz_list.append(dz)
                t_prev = t
            dtopo.dz_list = dz_list

        else:   
            raise ValueError("Unrecognized rupture_type: %s" % self.rupture_type)

        # Store for user
        self.dtopo = dtopo

        return dtopo


    
    def plot_subfaults(self, axes=None, plot_centerline=False, slip_color=False,
                             cmap_slip=None, cmin_slip=None, cmax_slip=None,
                             slip_time=None,
                             plot_rake=False, xylim=None, plot_box=True):
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
        print "Max slip, Min slip: ",max_slip, min_slip
    
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

            # unpack parameters:
            paramlist = """x_top y_top x_bottom y_bottom x_centroid y_centroid
                depth_top depth_bottom x_corners y_corners""".split()
    
            for param in paramlist:
                cmd = "%s = subfault.geometry['%s']" % (param,param)
                exec(cmd)
    
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
            axes.set_xlim(xylim[:2])
            axes.set_ylim(xylim[2:])
        if slip_color:
            if slip_time is None:
                axes.set_title('Slip on fault')
            else:
                axes.set_title('Slip on fault at time %6.1fs' % slip_time)
        else:
            axes.set_title('Fault planes')

        axes.ticklabel_format(format='plain', useOffset=False)
        labels = axes.get_xticks().tolist()
        axes.set_xticklabels(labels, rotation=80)

        if slip_color:
            cax,kw = matplotlib.colorbar.make_axes(slipax)
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
    
            # unpack parameters:
            paramlist = """x_top y_top x_bottom y_bottom x_centroid y_centroid
                depth_top depth_bottom x_corners y_corners""".split()
    
            for param in paramlist:
                cmd = "%s = subfault.geometry['%s']" % (param,param)
                exec(cmd)
    
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

        rect = [numpy.infty, -numpy.infty, numpy.infty, -numpy.infty]
        for subfault in self.subfaults:
            xmin = numpy.array(subfault.geometry['x_corners']).min()
            xmax = numpy.array(subfault.geometry['x_corners']).max()
            ymin = numpy.array(subfault.geometry['y_corners']).min()
            ymax = numpy.array(subfault.geometry['y_corners']).max()
            rect[0] = min(xmin, rect[0])
            rect[1] = max(xmax, rect[1])
            rect[2] = min(ymin, rect[2])
            rect[3] = max(ymax, rect[3])

        return rect

    
    def create_dtopo_xy(self, rect=None, dx=1/60., buffer_size=0.5):
        r"""Create coordinate arrays containing fault with a buffer.

        Input
        -----
         - *rect* - if None, use self.containing_rect
            Otherwise a list [x1,x2,y1,y2]
         - *dx* (int) - Spatial resolution. Defaults to 1" resolution.
         - *buffer_size* (float) - Buffer distance around edge of fault in
           degrees, defaults to 0.5 degrees.

        Output
        ------
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

    """

    def __init__(self):
        r"""SubFault initialization routine.
        
        See :class:`SubFault` for more info.

        """
        
        super(SubFault, self).__init__()

        self.strike = None
        r"""Strike direction of subfault in degrees."""
        self.length = None
        r"""Length of subfault in standard units."""
        self.width = None
        r"""Width of subfault in standard units."""
        self.depth = None
        r"""Depth of subfault based on *coordinate_specification*."""
        self.slip = None
        r"""Slip on subfault in strike direction in standard units."""
        self.rake = None
        r"""Rake of subfault movement in degrees."""
        self.dip = None
        r"""Subfault's angle of dip"""
        self.latitude = None
        r"""Latitutde of the subfault based on *coordinate_specification*."""
        self.longitude = None
        r"""Longitude of the subfault based on *coordinate_specification*."""
        self.coordinate_specification = "top center"
        r"""Specifies where the latitude, longitude and depth are measured from."""
        self.mu = 4e11  # default value for rigidity = shear modulus
        r"""Rigidity of subfault movement == shear modulus."""

        # deprecated:
        #self.units = {'mu':"dyne/cm^2", 'length':'km', 'width':'km', 
        #              'depth':'km'}
        #r"""Dictionary of units for the relevant parameters."""

        #self.units.update(units)


        self._geometry = None

    def convert_to_standard_units(self, input_units, verbose=False):
        r"""
        Convert parameters from the units used for input into the standard
        units used in this module.
        """
        params = standard_units.keys()
        for param in params:
            value = getattr(self, param)
            converted_value = convert_units(value, input_units[param], 1)
            setattr(self,param,converted_value)
            if verbose:
                print "%s %s %s converted to %s %s" \
                    % (param, value, input_units[param], converted_value, \
                       standard_units[param])

    def convert2meters(self, parameters): 
        r"""Convert relevant lengths to correct units.

        Returns converted (length, width, depth, slip) 
        Deprecated?

        """

        ## Note nm = Nautical mile
        conversion_dict = {"km":1e3, "cm":1e-2, "nm":1852.0, "m":1.0}
        
        converted_lengths = []
        for (n, parameter) in enumerate(parameters):
            if isinstance(getattr(self, parameter), list) or \
               isinstance(getattr(self, parameter), tuple):
                converted_lengths.append(list())
                for k in xrange(len(getattr(self, parameter))):
                    converted_lengths[n].append(getattr(self, parameter)[k] * \
                            conversion_dict[self.units[parameter]])
            else:
                converted_lengths.append(getattr(self, parameter) * \
                            conversion_dict[self.units[parameter]])

        if len(converted_lengths) == 1:
            return converted_lengths[0]

        return converted_lengths


    def Mo(self):
        r"""Calculate the seismic moment for a single subfault

        Returns in units of N-m and assumes mu is in Pascals. 
        """

        total_slip = self.length * self.width * self.slip
        Mo = self.mu * total_slip
        return Mo


    @property
    def geometry(self):
        r"""Subfault geometry"""
        if self._geometry is None:
            self.set_geometry()
        return self._geometry

    @geometry.setter
    def geometry(self,value):
        # Do we need this?
        self._geometry = value
    @geometry.deleter
    def geometry(self):
        del self._geometry


    def set_geometry(self):
        r"""
        Set self._geometry, a dictionary containing 
        bottom, top, centroid, and corner values of x,y, and depth at top and
        bottom of fault, based on subfault parameters.  
        Automatically called first time user requests self.geometry.

        Note: *self.coordinate_specification*  specifies the location on each
            subfault that corresponds to the (longitude,latitude) and depth 
            of the subfault.
            Currently must be one of these strings:
                "bottom center": (longitude,latitude) and depth at bottom center
                "top center": (longitude,latitude) and depth at top center
                "centroid": (longitude,latitude) and depth at centroid of plane
                "noaa sift": (longitude,latitude) at bottom center, depth at top,  
                             This mixed convention is used by the NOAA SIFT
                             database and "unit sources", see:
                             http://nctr.pmel.noaa.gov/propagation-database.html
            The Okada model is expressed assuming (longitude,latitude) and depth
            are at the bottom center of the fault plane, so values must be
            shifted or other specifications.
        """

        # Convert to meters if necessary:
        #length, width, depth, slip = self.convert2meters(["length","width", \
        #                            "depth","slip"])
        # Should now already be in meters!

        length = self.length
        width = self.width
        depth = self.depth
        slip = self.slip
        x0 =  self.longitude
        y0 =  self.latitude
        location =  self.coordinate_specification

        halfL = 0.5*length
        w  =  width

        # convert angles to radians:
        ang_dip = DEG2RAD * self.dip
        ang_rake = DEG2RAD * self.rake
        ang_strike = DEG2RAD * self.strike
    
        # vector (dx,dy) goes up-dip from bottom to top:
        dx = -w*numpy.cos(ang_dip)*numpy.cos(ang_strike) / \
                (LAT2METER*numpy.cos(y0*DEG2RAD))
        dy = w*numpy.cos(ang_dip)*numpy.sin(ang_strike) / LAT2METER

        if location == "bottom center":
            depth_bottom = depth
            depth_top = depth - w*numpy.sin(ang_dip)
            x_bottom = x0
            y_bottom = y0
            x_top = x0 + dx
            y_top = y0 + dy
            x_centroid = x_bottom + 0.5*dx
            y_centroid = y_bottom + 0.5*dy

        elif location == "top center":
            depth_top = depth
            depth_bottom = depth + w*numpy.sin(ang_dip)
            x_top = x0
            y_top = y0 
            x_bottom = x0 - dx
            y_bottom = y0 - dy
            x_centroid = x_bottom + 0.5*dx
            y_centroid = y_bottom + 0.5*dy
    
        elif location == "centroid":
            depth_top = depth - 0.5*w*numpy.sin(ang_dip)
            depth_bottom = depth + 0.5*w*numpy.sin(ang_dip)
    
            x_centroid = x0
            y_centroid = y0
            x_top = x0 + 0.5*dx
            y_top = y0 + 0.5*dy
            x_bottom = x0 - 0.5*dx
            y_bottom = y0 - 0.5*dy
    
        elif location == "noaa sift":
            depth_top = depth
            depth_bottom = depth + w*numpy.sin(ang_dip)
            x_bottom = x0
            y_bottom = y0
            x_top = x0 + dx
            y_top = y0 + dy
            x_centroid = x_bottom + 0.5*dx
            y_centroid = y_bottom + 0.5*dy

        else:
            raise ValueError("Unrecognized coordinate_specification" \
                    % coordinate_specification)
        

        # distance along strike from center of an edge to corner:
        dx2 = 0.5*length*numpy.sin(ang_strike) \
                / (LAT2METER*numpy.cos(y_bottom*DEG2RAD))
        dy2 = 0.5*length*numpy.cos(ang_strike) / LAT2METER
        x_corners = [x_bottom-dx2,x_top-dx2,x_top+dx2,x_bottom+dx2,x_bottom-dx2]
        y_corners = [y_bottom-dy2,y_top-dy2,y_top+dy2,y_bottom+dy2,y_bottom-dy2]

        # restore proper units to depth if necessary:
        # deprecated
        #if self.units['depth'] == 'km':
            #depth_top = depth_top / 1000.
            #depth_bottom = depth_bottom / 1000.

        paramlist = """x_top y_top x_bottom y_bottom x_centroid y_centroid
            depth_top depth_bottom x_corners y_corners""".split()

        self._geometry = {}
        for param in paramlist:
            cmd = "self._geometry['%s'] = %s" % (param,eval(param))
            exec(cmd)
    

    def okada(self, x, y):
        r"""
        Apply Okada to this subfault and return a DTopography object.

        Input:
            x,y are 1d arrays
        Output:
            DTopography object with dz_list = [dz] being a list with 
                single static displacement and times = [0.].

        Currently only calculates the vertical displacement.

        Okada model is a mapping from several fault parameters
        to a surface deformation.
        See Okada 1985, or Okada 1992, Bull. Seism. Soc. Am.
        
        okadamap function riginally written in Python by Dave George for
        Clawpack 4.6 okada.py routine, with some routines adapted
        from fortran routines written by Xiaoming Wang.

        Rewritten and made more flexible by Randy LeVeque

        Note: *self.coordinate_specification* (str) specifies the location on each
            subfault that corresponds to the (longitude,latitude) and depth 
            of the subfault.
        See the documentation for *SubFault.set_geometry* for dicussion of the 
        possible values *self.coordinate_specification* can take.

        """

        # Okada model assumes x,y are at bottom center:
        x_bottom = self.geometry['x_bottom']
        y_bottom = self.geometry['y_bottom']
        depth_bottom = self.geometry['depth_bottom']


        # Convert some parameters to meters if necessary:
        # deprecated -- assume they are in meters
        #if self.units['depth'] == 'km':
        #    depth_bottom = depth_bottom * 1000.
        #length, width, depth, slip = self.convert2meters(["length","width", \
        #                            "depth","slip"])

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
    
        X,Y = numpy.meshgrid(x,y)   # use convention of upper case for 2d
    
        # Convert distance from (X,Y) to (x_bottom,y_bottom) from degrees to
        # meters:
        xx = LAT2METER * numpy.cos(DEG2RAD*Y)*(X-x_bottom)   
        yy = LAT2METER * (Y-y_bottom)
    
    
        # Convert to distance along strike (x1) and dip (x2):
        x1 = xx*numpy.sin(ang_strike) + yy*numpy.cos(ang_strike) 
        x2 = xx*numpy.cos(ang_strike) - yy*numpy.sin(ang_strike) 
    
        # In Okada's paper, x2 is distance up the fault plane, not down dip:
        x2 = -x2
    
        p = x2*numpy.cos(ang_dip) + depth_bottom * numpy.sin(ang_dip)
        q = x2*numpy.sin(ang_dip) - depth_bottom * numpy.cos(ang_dip)
    
        f1=self._strike_slip (x1+halfL,p,  ang_dip,q)
        f2=self._strike_slip (x1+halfL,p-w,ang_dip,q)
        f3=self._strike_slip (x1-halfL,p,  ang_dip,q)
        f4=self._strike_slip (x1-halfL,p-w,ang_dip,q)
    
        g1=self._dip_slip (x1+halfL,p,  ang_dip,q)
        g2=self._dip_slip (x1+halfL,p-w,ang_dip,q)
        g3=self._dip_slip (x1-halfL,p,  ang_dip,q)
        g4=self._dip_slip (x1-halfL,p-w,ang_dip,q)
    
        # Displacement in direction of strike and dip:
        ds = slip*numpy.cos(ang_rake)
        dd = slip*numpy.sin(ang_rake)
    
        us = (f1-f2-f3+f4)*ds
        ud = (g1-g2-g3+g4)*dd
    
        dz = (us+ud)

        dtopo = DTopography()
        dtopo.X = X
        dtopo.Y = Y
        dtopo.dz_list = [dz]
        dtopo.times = [0.]
        self.dtopo = dtopo
        return dtopo

    # Utility functions for okada:

    def _strike_slip(self,y1,y2,ang_dip,q):
        """
        !.....Used for Okada's model
        !.. ..Methods from Yoshimitsu Okada (1985)
        !-----------------------------------------------------------------------
        """
        sn = numpy.sin(ang_dip)
        cs = numpy.cos(ang_dip)
        d_bar = y2*sn - q*cs
        r = numpy.sqrt(y1**2 + y2**2 + q**2)
        xx = numpy.sqrt(y1**2 + q**2)
        a4 = 2.0*poisson/cs*(numpy.log(r+d_bar) - sn*numpy.log(r+y2))
        f = -(d_bar*q/r/(r+y2) + q*sn/(r+y2) + a4*sn)/(2.0*3.14159)
    
        return f
    
    
    def _dip_slip(self,y1,y2,ang_dip,q):
        """
        !.....Based on Okada's paper (1985)
        !.....Added by Xiaoming Wang
        !-----------------------------------------------------------------------
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
            *rupture_time*
            *rise_time*
            *rise_time_ending*: optional, defaults to *rise_time*

        """
        if (self.rupture_time is None) or (self.rise_time is None):
            raise Exception("dynamic_slip method only works for dynamic ruptures")

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

    def __init__(self):
        r"""UCSBFault initialization routine.
        
        See :class:`UCSBFault` for more info.

        """

        self.num_cells = [None, None]   # RJL: Why needed??

        super(UCSBFault, self).__init__()


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
        # a schematic of where these points are
        self._fault_plane_corners = [None, # a 
                                     None, # b
                                     None, # c
                                     None] # d
        self._fault_plane_centers = [[0.0, 0.0, 0.0], # 1
                                     [0.0, 0.0, 0.0], # 2 
                                     [0.0, 0.0, 0.0]] # 3
        # :TODO: Is the order of this a good assumption?
        self._fault_plane_corners[0] = boundary_data[0]
        self._fault_plane_corners[3] = boundary_data[1]
        self._fault_plane_corners[2] = boundary_data[2]
        self._fault_plane_corners[1] = boundary_data[3]
        
        # Calculate center by averaging position of appropriate corners
        for (n, corner) in enumerate(self._fault_plane_corners):
            for i in xrange(3):
                self._fault_plane_centers[1][i] += corner[i] / 4
            if n == 0 or n == 4:
                for i in xrange(3):
                    self._fault_plane_centers[0][i] += corner[i] / 2
            else:
                for i in xrange(3):
                    self._fault_plane_centers[2][i] += corner[i] / 2


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
                         rupture_type="static"):
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
                    if input_units.get(column_name,units) != units:
                        print "*** Warning: input_units[%s] reset to %s" \
                              % (column_name, units)
                        print "    based on file header"
                        input_units[column_name] = units

                else:
                    column_name = column_heading.lower()
                column_name = column_name.strip()
                if column_name in param.keys():
                    column_key = param[column_name]
                    column_map[column_key] = n
                else:
                    print "*** Warning: column name not recognized: %s" \
                        % column_name

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
        for k,v in sift_slip.iteritems():
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
    *fault_plane* of class *SubFault*) and subdividing the fault plane
    into a rectangular array of *nstrike* by *ndip* equally sized subfaults.

    By default,the slip on each subfault will be initialized to
    *fault_plane.slip* so that the slip is uniform over the original plane
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


    def __init__(self, fault_plane, nstrike=1, ndip=1,
                       slip_function=None, Mo=None):
        
        from numpy import sin,cos
        super(SubdividedPlaneFault, self).__init__()

        self.fault_plane = fault_plane
        self.nstrike = nstrike
        self.ndip = ndip

        self.subdivide(nstrike, ndip, slip_function, Mo)


    def subdivide(self, nstrike=1, ndip=1, slip_function=None, Mo=None):

        # may have changed resolution:
        self.nstrike = nstrike
        self.ndip = ndip

        fault_plane = self.fault_plane
        strike = fault_plane.strike
        dip = fault_plane.dip
        rake = fault_plane.rake
        slip = fault_plane.slip
        length = fault_plane.length
        width = fault_plane.width

        # unpack corners from fault plane geometry:
        x_corners = fault_plane.geometry['x_corners']
        y_corners = fault_plane.geometry['y_corners']

        # set depth at corners:
        depth_top = fault_plane.geometry['depth_top']
        depth_bottom = fault_plane.geometry['depth_bottom']
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

                if fault_plane.coordinate_specification == 'centroid':
                    subfault.longitude = x_sf[1]
                    subfault.latitude = y_sf[1]
                    subfault.depth = d_sf[1]
                elif fault_plane.coordinate_specification == 'top center':
                    subfault.longitude = x_sf[0]
                    subfault.latitude = y_sf[0]
                    subfault.depth = d_sf[0]
                elif fault_plane.coordinate_specification == 'noaa sift':
                    subfault.longitude = x_sf[2]
                    subfault.latitude = y_sf[2]
                    subfault.depth = d_sf[0]
                else:   
                    msg = "Unrecognized coordinate_specification: %s" \
                            % fault_plane.coordinate_specification
                    raise NotImplementedError(msg)

                subfault.dip = dip
                subfault.strike = strike
                subfault.rake = rake
                subfault.length = length / nstrike
                subfault.width = width / ndip
                subfault.slip = slip
                subfault.coordinate_specification = \
                        fault_plane.coordinate_specification
                subfault.mu = fault_plane.mu

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
        
        # perform the subdivision and set parameters on each subfault:
        super(TensorProductFault, self).__init__(fault_plane, nstrike, ndip)

        if slip_along_strike is None:
            # set to constant in the strike direction if not specified
            slip_along_strike = lambda eta: 1.0
        if slip_down_dip is None:
            # set to constant in the dip direction if not specified
            slip_down_dip = lambda xi: 1.0

        

# =============================================================================
# Below are the first attempt at an implementation for these classes
# #  Fault Class 
# # ==============================================================================
# class Fault(object):


#     r"""Generic class representing a fault.

#     :TODO:
#      - Support something other than lat-long
#      - Provide plots (and other plot types)
#      - Provide detailed documentation


#     Fault Parameters
#     -------------------


#     Attributes
#     ----------
#      - *units* (dict) - Dictionary containing unit specifications for the 

#     Properties
#     ----------
#      :Note: All properties are in meters and do not match the units dictionary.

#     """

#     @property
#     def x(self):
#         r"""Coordinate array (x) for fault."""
#         if self._x is None:
#             self.create_coordinate_arrays()
#         return self._x
#     @x.setter
#     def x(self, value):
#         self._x = value
#     @x.deleter
#     def x(self):
#         del self._x
#     @property
#     def X(self):
#         r"""2d x-coordinate array."""
#         if self._X is None:
#             self.create_2d_coordinate_arrays()
#         return self._X
#     @X.deleter
#     def X(self):
#         del self._X

#     @property
#     def y(self):
#         r"""Coordinate array (y) for fault."""
#         if self._y is None:
#             self.create_coordinate_arrays()
#         return self._y
#     @y.setter
#     def y(self, value):
#         self._y = value
#     @y.deleter
#     def y(self):
#         del self._y
#     @property
#     def Y(self):
#         r"""2d y-coordinate array."""
#         if self._Y is None:
#             self.create_2d_coordinate_arrays()
#         return self._Y
#     @Y.deleter
#     def Y(self):
#         del self._Y

#     @property
#     def dZ(self):
#         r"""Deformation dZ of fault."""
#         if self._dZ is None:
#             self.create_deformation_array()
#         return self._dZ
#     @dZ.setter
#     def dZ(self, value):
#         self._dZ = value
#     @dZ.deleter
#     def dZ(self):
#         del self._dZ

#     # Calculated geometry
#     @property
#     def fault_plane_corners(self):
#         r"""Coordinates of the corners of the fault plane"""
#         if self._fault_plane_corners is None:
#             self.calculate_geometry()
#         return self._fault_plane_corners
#     @property
#     def fault_plane_centers(self):
#         r"""Coordinates along the center-line of the fault plane"""
#         if self._fault_plane_centers is None:
#             self.calculate_geometry()
#         return self._fault_plane_centers

#     # Earthquake calculated properties
#     @property
#     def slip(self):
#         r"""Total slip accumulated over all subfaults."""
#         if self._slip is None:
#             self._slip = 0.0
#             for subfault in self.subfaults:
#                 slip = subfault.convert2meters(['slip'])
#                 self._slip += slip
#             # Convert back to units of this object if available
#             # import pdb; pdb.set_trace()
#             if self.units.has_key('slip'):
#                 self._slip = self.convert2meters(['slip'])
#             else:
#                 self.units['slip'] = self.subfaults[0].units['slip']
#         return self._slip

#     @property
#     def Mw(self):
#         r"""Calculate the effective moment magnitude of all subfaults.

#         :TODO:
#          - Need to check that this does the right thing
#          - Need to check mu units

#         """
#         total_Mo = 0.0
#         for subfault in self.subfaults:
#             dimensions, slip = subfault.convert2meters(["dimensions", "slip"])
#             if subfault.units['mu'] in ['dyne/cm^2', "dyne/cm^2"]:
#                 mu = subfault.mu * 0.1
#                 # mu = subfault.mu * 100.0**2
#             else:
#                 mu = subfault.mu
#             # All units should be terms of meters
#             total_Mo += dimensions[0] * dimensions[1] * slip * mu

#         return 2.0 / 3.0 * (numpy.log10(total_Mo) - 9.05)

#     @property
#     def dimensions(self):
#         if self._dimensions is None:
#             self._dimensions = [0.0, 0.0]
#             for subfault in self.subfaults:
#                 local_dimensions = subfault.convert2meters(["dimensions"])
#                 self._dimensions[0] += local_dimensions[0]
#                 self._dimensions[1] += local_dimensions[1]
#         return self._dimensions


#     def __init__(self, path=None, subfaults=None, units={}):
#         r"""Fault initialization routine.
        
#         See :class:`Fault` for more info.

#         """

#         super(Fault, self).__init__()

#         # Object defered storage
#         self._x = None
#         self._X = None
#         self._y = None
#         self._Y = None
#         self._dZ = None
#         self._slip = None
#         self._fault_plane_corners = None
#         self._fault_plane_centers = None
#         self._dimensions = None

#         # Parameters for subfault specification
#         self.rupture_type = 'static' # 'static', 'dynamic', 'kinetic'
#         self.t = numpy.array([0.0, 5.0, 10.0])

#         # Default units of each parameter type
#         self.units = {}
#         self.units.update(units)
        
#         if path is not None:
#             # Read in file at path assuming it is a subfault specification
#             self.read(path)
#         elif subfaults is not None:
#             if not isinstance(subfaults, list):
#                 raise ValueError("Input parameter subfaults must be a list.")
#             self.subfaults = subfaults


#     def __str__(self):
#         output =  "Fault Characteristics:\n"
#         output += "  Mw = %s\n" % self.Mw
#         max_slip = 0.0
#         min_slip = numpy.infty
#         for subfault in self.subfaults:
#             slip = subfault.convert2meters(["slip"])
#             max_slip = max(max_slip, slip)
#             min_slip = min(min_slip, slip)
#         slip_unit = self.units.get('slip', self.subfaults[0].units['slip'])
#         output += "  Slip (Max, Min, Total) %s = (%s, %s, %s)" % (
#                                                              slip_unit,
#                                                              max_slip, 
#                                                              min_slip,
#                                                              self.slip)
#         return output


#     def read_dtopo(self, path, topo_type=3):
#         r"""Read in dtopo file at *path*.

#         """

#         if topo_type == 1:    
#             # Load raw data
#             data = numpy.loadtxt(path)

#             # Parse data
#             t = data[:,0]
#             x = data[:,1]
#             y = data[:,2]

#             # Initialize extents
#             t0 = t[0]
#             lower = [x[0], y[0]]
#             upper = [None, None]
#             num_cells = [0,0]

#             # Count total x-values
#             for row in xrange(1,data.shape[0]):
#                 if x[row] == x[0]:
#                     num_cells[0] = row
#                     break

#             # Count total y-values
#             for row in xrange(num_cells[0], data.shape[0]):
#                 if t[row] != t0:
#                     num_cells[1] = row / num_cells[0]
#                     num_times = data.shape[0] / row
#                     break

#             # Check extents
#             assert(t[0] != t[num_cells[0] * num_cells[1] + 1])
#             # assert(t[0] == t[num_cells[0] * num_cells[1]])

#             # Fill in rest of pertinent data
#             self.t = data[::num_cells[0] * num_cells[1], 0]
#             self._x = data[:num_cells[0], 1]
#             self._yy = data[:num_cells[0] * num_cells[1]:num_cells[0], 2]
#             upper = [x[-1], y[-1]]
#             self._X, self._Y = numpy.meshgrid(x, y)
#             self._dZ = numpy.empty( (num_times, num_cells[0], num_cells[1]) )

#             for (n,time) in enumerate(t):
#                 self._dZ[n,:,:] = data[num_cells[0] * num_cells[1] * n:
#                                 num_cells[0] * num_cells[1] * (n+1), 3].reshape(
#                                         (num_cells[0], num_cells[1]))

#         elif topo_type == 2 or topo_type == 3:
        
#             # Read header
#             mx = int(dtopo_file.readline().split()[0])
#             my = int(dtopo_file.readline().split()[0])
#             mt = int(dtopo_file.readline().split()[0])
#             xlower = float(dtopo_file.readline().split()[0])
#             ylower = float(dtopo_file.readline().split()[0])
#             t0 = float(dtopo_file.readline().split()[0])
#             dx = float(dtopo_file.readline().split()[0])
#             dy = float(dtopo_file.readline().split()[0])
#             dt = float(dtopo_file.readline().split()[0])

#             # Construct coordinate arrays
#             self._x = numpy.linspace(xlower, xlower + (mx - 1) * dx, mx)
#             self._y = numpy.linspace(ylower, ylower + (my - 1) * dy, my) 
#             self._dZ = numpy.empty((mx,my,mt))
#             times = numpy.linspace(t0, t0+(mt-1)*dt, mt)

#             raise NotImplementedError("Reading dtopo type 2 or 3 files not implemented.")

#             # Read data
#             if topo_type == 2:
#                 pass

#             elif topo_type == 3:
#                 pass

#         else:
#             raise ValueError("Topo type %s is invalid." % topo_type)



#     def read(self, path, column_map, coordinate_specification="centroid",
#                          rupture_type="static", t=[0.0, 0.5, 1.0], skiprows=0,
#                          delimiter=None):
#         r"""Read in subfault specification at *path*.

#         Creates a list of subfaults from the subfault specification file at
#         *path*.

#         column_map = {"coordinates":(1,0), "depth":2, "slip":3, "rake":4, 
#                           "strike":5, "dip":6}

#         """

#         # Read in rest of data
#         data = numpy.loadtxt(path, skiprows=skiprows, delimiter=delimiter)
#         if len(data.shape) == 1:
#             data = numpy.array([data])

#         self.subfaults = []
#         for n in xrange(data.shape[0]):

#             new_subfault = SubFault()
#             new_subfault.coordinate_specification = coordinate_specification
#             new_subfault.rupture_type = rupture_type
#             new_subfault.t = numpy.array(t) 

#             for (var, column) in column_map.iteritems():
#                 if isinstance(column, tuple) or isinstance(column, list):
#                     setattr(new_subfault, var, [None for k in column])
#                     for (k, index) in enumerate(column):
#                         getattr(new_subfault, var)[k] = data[n, index]
#                 else:
#                     setattr(new_subfault, var, data[n, column])

#             self.subfaults.append(new_subfault)


#     def convert2meters(self, parameters): 
#         r"""Convert relevant lengths to correct units.

#         Returns converted (dimensions, depth, slip) 

#         """

#         conversion_dict = {"km":1e3, "cm":1e-2, "nm":1852.0, "m":1.0}
        
#         converted_lengths = []
#         for (n, parameter) in enumerate(parameters):
#             if isinstance(getattr(self, parameter), list) or \
#                isinstance(getattr(self, parameter), tuple):
#                 converted_lengths.append(list())
#                 for k in xrange(len(getattr(self, parameter))):
#                     converted_lengths[n].append(getattr(self, parameter)[k] * conversion_dict[self.units[parameter]])
#             else:
#                 converted_lengths.append(getattr(self, parameter) * conversion_dict[self.units[parameter]])

#         if len(converted_lengths) == 1:
#             return converted_lengths[0]

#         return converted_lengths


#     def containing_rect(self):
#         r"""Find containing rectangle of subfault in x-y plane.

#         Returns tuple of x-limits and y-limits.

#         """

#         extent = [numpy.infty, -numpy.infty, numpy.infty, -numpy.infty]
#         for corner in self.fault_plane_corners:
#             extent[0] = min(corner[0], extent[0])
#             extent[1] = max(corner[0], extent[1])
#             extent[2] = min(corner[1], extent[2])
#             extent[3] = max(corner[1], extent[3])

#         return extent


#     def transform(self, origin, point, theta):
#         r"""Transform to rotated coordinate system (via strike).

#         Input
#         -----
#          - *origin* (tuple) - Point being rotated about, (x,y).
#          - *point* (tuple) - Point to be rotated, (x,y).
#          - *theta* (float) - Angle of rotation in radians.

#         Output
#         ------
#          - (tuple) The result of the transforming *point* about the point
#            *origin*.

#         """
#         return (  (point[0] - origin[0]) * numpy.cos(-theta) 
#                 - (point[1] - origin[1]) * numpy.sin(-theta) 
#                 + origin[0],
#                   (point[0] - origin[0]) * numpy.sin(-theta) 
#                 + (point[1] - origin[1]) * numpy.cos(-theta) 
#                 + origin[1]) 


#     def calculate_geometry(self):
#         r"""Calculate the fault geometry.

#         Routine calculates the class attributes *fault_plane_corners* and 
#         *fault_plane_centers* which are the corners of the fault plane and 
#         points along the centerline respecitvely in 3D space.

#         Note that these are simply the largest containing rectangle and the
#         average location of the centers.  Really this should be overridden in a
#         subclass by a particular type of Fault.

#         """

#         # Initialize corners with first subfault
#         self._fault_plane_corners = self.subfaults[0].fault_plane_corners
#         self._fault_plane_centers = self.subfaults[0].fault_plane_centers

#         # If there's only one subfault we are done, otherwise continue
#         if len(self.subfaults) > 1:
#             raise NotImplementedError("Calculating fault geometry not ",
#                                       "implemented.")
    

#     def create_coordinate_arrays(self, resolution=60, buffer_size=0.5):
#         r"""Create coordinate arrays containing subfault.

#         Input
#         -----
#          - *resolution* (int) - Number of grid points per degree.  Defaults to
#            1" resolution.
#          - *buffer_size* (float) - Buffer distance around edge of fault in 
#            degrees, defaults to 0.5 degrees.

#         """

#         rect = self.containing_rect()
#         rect[0] -= buffer_size
#         rect[1] += buffer_size
#         rect[2] -= buffer_size
#         rect[3] += buffer_size

#         self.delta = float(1.0 / resolution)
#         N = [int((rect[1] - rect[0]) * resolution),
#              int((rect[3] - rect[2]) * resolution)]

#         self._x = numpy.linspace(rect[0],rect[1],N[0])
#         self._y = numpy.linspace(rect[2],rect[3],N[1])


#     def create_2d_coordinate_arrays(self):
#         r"""Create 2d-coodrinate arrays."""
        
#         self._X, self._Y = numpy.meshgrid(self.x, self.y)


#     def create_deformation_array(self):
#         r"""Create deformation array dZ.

#         Use subfaults' `create_deformation_array` routine and add all 
#         deformations together.

#         """

#         self._dZ = numpy.zeros(self.X.shape)
#         for subfault in self.subfaults:
#             subfault_dZ = subfault.create_deformation_array(domain=(self.x, self.y))
#             self._dZ += subfault_dZ

#         return self._dZ

#     def write(self, path, topo_type=None):
#         r"""Write out subfault resulting dtopo to file at *path*.

#         input
#         -----
#          - *path* (path) - Path to the output file to written to.
#          - *topo_type* (int) - Type of topography file to write out.  Default is 1.

#         """

#         if topo_type is None:
#             # Try to look at suffix for type
#             extension = os.path.splitext(path)[1][1:]
#             if extension[:2] == "tt":
#                 topo_type = int(extension[2])
#             elif extension == 'xyz':
#                 topo_type = 1
#             else:
#                 # Default to 3
#                 topo_type = 3

#         if self.rupture_type != "static":
#             raise NotImplemented("Only the 'static' rupture type is supported.")

#         # Construct each interpolating function and evaluate at new grid
#         with open(path, 'w') as data_file:

#             if topo_type == 1:
#                 # Topography file with 4 columns, t, x, y, dz written from the upper
#                 # left corner of the region
#                 Y_flipped = numpy.flipud(self.Y)
#                 for (n, time) in enumerate(self.t):
#                     alpha = (time - self.t[0]) / self.t[-1]
#                     dZ_flipped = numpy.flipud(alpha * self.dZ[:,:])
#                     for j in xrange(self.Y.shape[0]):
#                         for i in xrange(self.X.shape[1]):
#                             data_file.write("%s %s %s %s\n" % (self.t[n], self.X[j,i], Y_flipped[j,i], dZ_flipped[j,i]))
        
#             elif topo_type == 2 or topo_type == 3:
#                 # Write out header
#                 data_file.write("%7i       mx \n" % self.x.shape[0])
#                 data_file.write("%7i       my \n" % self.y.shape[0])
#                 data_file.write("%7i       mt \n" % self.t.shape[0])
#                 data_file.write("%20.14e   xlower\n" % self.x[0])
#                 data_file.write("%20.14e   ylower\n" % self.y[0])
#                 data_file.write("%20.14e   t0\n" % self.t[0])
#                 data_file.write("%20.14e   dx\n" % self.delta)
#                 data_file.write("%20.14e   dy\n" % self.delta)
#                 data_file.write("%20.14e   dt\n" % float(self.t[1] - self.t[0]))

#                 if topo_type == 2:
#                     raise ValueError("Topography type 2 is not yet supported.")
#                 elif topo_type == 3:
#                     for (n, time) in enumerate(self.t):
#                         alpha = (time - self.t[0]) / (self.t[-1])
#                         for j in range(self.Y.shape[0]-1, -1, -1):
#                             data_file.write(self.X.shape[1] * '%012.6e  ' 
#                                                   % tuple(alpha * self.dZ[j,:]))
#                             data_file.write("\n")

#             else:
#                 raise ValueError("Only topography types 1, 2, and 3 are supported.")


#     def plot(self, axes=None, region_extent=None, contours=None, 
#                    coastlines=None, limits=None, cmap=None):
#         r"""Plot subfault deformation.


#         Input
#         -----
#          - *axes* (`matplotlib.axes.Axes`) -
#          - *region_extent* (list) - 
#          - *contours* (list) -
#          - *coastlines* (path) -
#          - *limits* (list) -
#          - *cmap* (`matplotlib.colors.Colormap`) -

#         Output
#         ------
#          - *axes* (`matplotlib.axes.Axes`) - Axes used for plot, either
#            this is created in this method if the input argument *axes* is None
#            or the same object is passed back.

#         """

#         import matplotlib.pyplot as plt
#         import matplotlib.colors as colors
#         import clawpack.visclaw.colormaps as colormaps
        
#         # Create axes object if needed
#         if axes is None:
#             fig = plt.figure()
#             axes = fig.add_subplot(1,1,1)

#         # Calculate plot extent and limits if not provided
#         if region_extent is None:
#             region_extent = ( numpy.min(self.x), numpy.max(self.x),
#                               numpy.min(self.y), numpy.max(self.y) )
#         if limits is None:
#             depth_extent = [numpy.min(self.dZ),numpy.max(self.dZ)]
#         else:
#             depth_extent = limits

#         # Setup axes labels, ticks and aspect
#         axes.ticklabel_format(format="plain", useOffset=False)
#         mean_lat = 0.5 * (region_extent[3] - region_extent[2])
#         axes.set_aspect(1.0 / numpy.cos(numpy.pi / 180.0 * mean_lat))
#         axes.set_title("Subfault Deformation")
#         axes.set_xlabel("Longitude")
#         axes.set_ylabel("Latitude")

#         # Colormap and color norm
#         if cmap is None:
#             if depth_extent[0] >= 0.0:
#                 cmap = colormaps.make_colormap({0.0:'w', 1.0:'r'})
#                 extend = 'top'
#             elif depth_extent[1] <= 0.0:
#                 cmap = colormaps.make_colormap({0.0:'b', 1.0:'w'})
#                 extend = 'bottom'
#             else:
#                 cmap = colormaps.make_colormap({0.0:'b', 0.5:'w', 1.0:'r'})
#                 extend = 'both'
#         # Equalize color extents
#         if depth_extent[0] >= 0.0:
#             depth_extent[0] = 0.0
#         elif depth_extent[1] <= 0.0:
#             depth_extent[1] = 0.0
#         else:
#             depth_extent[1] = max(-depth_extent[0], depth_extent[1])
#             depth_extent[0] = -depth_extent[1]
#         color_norm = colors.Normalize(depth_extent[0], depth_extent[1], clip=True)

#         # Plot data
#         if contours is not None:
#             plot = axes.contourf(self.x, self.y, self.dZ, contours, cmap=cmap,
#                                  extend=extend)
#         else:
#             plot = axes.imshow(numpy.flipud(self.dZ), vmin=depth_extent[0], 
#                                                       vmax=depth_extent[1],
#                                                       extent=region_extent,
#                                                       cmap=cmap,
#                                                       norm=color_norm)

#         cbar = plt.colorbar(plot, ax=axes)
#         cbar.set_label("Deformation (m)")

#         # Plot coastlines
#         if coastlines is not None:
#             coastline_data = topotools.Topography(coastlines)
#             axes.contour(coastline_data.X, coastline_data.Y, 
#                          coastline_data.Z, levels=[0.0],colors='r')

#         axes.set_xlim(region_extent[0:2])
#         axes.set_ylim(region_extent[2:])

#         return axes


#     def plot_fault_rect(self, axes=None, color='r', markerstyle="o", 
#                                          linestyle='-'):
#         r"""Plot fault rectangle.

#         Input
#         -----
#          - *axes* (`matplotlib.axes.Axes`) - 

#         Output
#         ------
#          - (`matplotlib.axes.Axes`) - 

#         """

#         import matplotlib.pyplot as plt
        
#         # Create axes object if needed
#         if axes is None:
#             fig = plt.figure()
#             axes = fig.add_subplot(1,1,1)

#         # Plot corners
#         style = color + markerstyle
#         for (n,corner) in enumerate(self.fault_plane_corners):
#             axes.plot(corner[0], corner[1], style)
#             axes.text(corner[0], corner[1], str(n+1))

#         # Plot edges
#         style = color + linestyle
#         edges = []
#         for edge in xrange(len(self.fault_plane_corners) - 1):
#             edges.append([self.fault_plane_corners[edge][:2], 
#                           self.fault_plane_corners[edge+1][:2]])
#         edges.append([self.fault_plane_corners[-1][:2], 
#                       self.fault_plane_corners[0][:2]])
#         for edge in edges:
#             axes.plot([edge[0][0], edge[1][0]], [edge[0][1], edge[1][1]], style)

#         return axes


#     def plot_contours(self, axes=None, dz_interval=0.5):
#         r"""Plot sea-floor deformation contours

#         """

#         dzmax = numpy.max(self.dZ.max(), -self.dZ.min()) + dz_interval
#         clines1 = numpy.arange(dz_interval, dzmax, dz_interval)
#         clines = list(-numpy.flipud(clines1)) + list(clines1)

#         self.plot(axes=axes, contours=clines)

#     def plot_seafloor(self, axes=None, cmax_dz=None, dz_interval=None):
#         r"""Plot sea floor deformation dtopo as colormap
        
#         """

#         self.plot(axes=axes)



# # ==============================================================================
# #  Sub-Fault Class 
# # ==============================================================================
# class SubFault(Fault):

#     r"""Class representing a single subfault.

#     :TODO:
#      - Support something other than lat-long
#      - Provide plots (and other plot types)
#      - Provide detailed documentation


#     Subfault Parameters
#     -------------------
#      - *dimensions* (list) - Dimensions of the fault plane.
#      - *coordinates* (list) - Longitude-latitude of some point on the fault 
#        plane.  Point is specified by *coordinate_specification*.
#      - *coordinate_specification* (string) - Specifies location relative to the
#        fault of the *coordinates* values.  Valid options include "top center",
#        "bottom center", and "centroid".
#      - *depth* (float) - Depth of the specified point below the sea-floor.
#      - *strike* (float) - Orientation of the top edge, measured in degrees 
#        clockwise from North.  Between 0 and 360. The fault plane dips downward 
#        to the right when moving along the top edge in the strike direction.
#      - *dip* (float) - Angle at which the plane dips downward from the top edge,
#        a positive angle between 0 and 90 degrees.
#      - *rake* (float) - Angle in the fault plane in which the slip occurs,
#        measured in degrees counterclockwise from the strike direction. Between 
#        -180 and 180.
#      - *slip* (float) Positive distance the hanging block moves relative to the 
#        foot block in the direction specified by the rake. The "hanging block" is
#        the one above the dipping fault plane (or to the right if you move in the 
#        strike direction).

#     Attributes
#     ----------
#      - *units* (dict) - Dictionary containing unit specifications for the 
#        *coordinates*, *dimenstions*, *slip*, and *depth* subfault parameters.  
#        Defaults to "lat-long", "m", "m", "m" respectively.

#     Properties
#     ----------
#      :Note: All properties are in meters and do not match the units dictionary.

#     """

#     @property
#     def slip(self):
#         r"""Amount of slip on this subfault."""
#         return self._slip
#     @slip.setter
#     def slip(self, value):
#         self._slip = value
#     @slip.deleter
#     def slip(self):
#         del self._slip

#     @property
#     def Mw(self):
#         r"""Calculate the effective moment magnitude of subfault."""
        
#         if self.units["mu"] == "dyne/cm^2":
#             mu = self.mu
#         elif self.units["mu"] == 'dyne/m^2':
#             mu = self.mu / 1e2**2
#         else:
#             raise ValueError("Unknown unit for rigidity %s." % self.units['mu'])

#         dimensions, slip = self.convert2meters(["dimensions", "slip"])
#         total_slip = dimensions[0] * dimensions[1] * slip
#         Mo = 0.1 * mu * total_slip
#         return 2.0 / 3.0 * (numpy.log10(Mo) - 9.05)

#     @property
#     def dimensions(self):
#         r"""Amount of slip on this subfault."""
#         return self._dimensions
#     @dimensions.setter
#     def dimensions(self, value):
#         self._dimensions = value
#     @dimensions.deleter
#     def dimensions(self):
#         del self._dimensions


#     def __init__(self, path=None, units={}):
#         r"""SubFault initialization routine.
        
#         See :class:`SubFault` for more info.

#         """

#         super(SubFault, self).__init__()

#         # Parameters for subfault specification
#         self.coordinates = [] # longitude, latitude
#         self.coordinate_specification = 'centroid' # 'centroid', 'top center', 'epicenter'
#         self._dimensions = [0.0, 0.0] # [length, width]
#         self.rake = None
#         self.strike = None
#         self.dip = None
#         self.depth = None
#         self._slip = None
#         self.mu = 5e11
#         self.t = numpy.array([0.0, 5.0, 10.0])

#         # Default units of each parameter type
#         self.units = {'coordinates':'lat-long', 'dimensions':'m', 'slip':'m',
#                       'depth':"km", "mu":"dyne/cm^2"}
#         self.units.update(units)

#         # Read in file at path if provided
#         if path is not None:
#             self.read(path)


#     def __str__(self):
#         output = "Subfault Characteristics:\n"
#         output += "  Coordinates: %s (%s)\n" % (self.coordinates, self.coordinate_specification)
#         output += "  Dimensions (L,W): %s %s\n" % (self.dimensions, self.units['dimensions'])
#         output += "  Depth: %s %s\n" % (self.depth, self.units["depth"])
#         output += "  Rake, Strike, Dip: %s, %s, %s\n" % (self.rake, self.strike, self.dip)
#         output += "  Slip, Mw: %s %s, %s\n" % (self.slip, self.units['slip'], self.Mw)
#         return output


#     def read(self, path):
#         r"""Read in subfault specification in from file at *path*

#         """
#         raise NotImplementedError("Reading of subfaults not implemented yet.")
        

#     def calculate_slip(self, Mw):
#         r"""Set slip based on a moment magnitude *Mw*."""
#         if self.units["mu"] == "dyne/cm^2":
#             mu = self.mu * 100.0**2

#         Mo = 10.**(Mw * 3.0 / 2.0 + 9.05)

#         dimensions = self.convert2meters("dimensions")
#         subfault_area = dimensions[0] * dimensions[1]
        
#         self.slip = Mo / (0.1 * mu * subfault_area)

#         # Convert back to requested units
#         if self.units['slip'] == 'cm':
#             self.slip *= 1e2
#         elif self.units["slip"] == 'km':
#             self.slip *= 1e-3


#     def calculate_geometry(self):
#         r"""Calculate the fault geometry.

#         Routine calculates the class attributes *fault_plane_corners* and 
#         *fault_plane_centers* which are the corners of the fault plane and 
#         points along the centerline respecitvely in 3D space.

#         """

#         # Simple conversion factor of latitude to meters
#         lat2meter = topotools.dist_latlong2meters(0.0, 1.0)[1]

#         # Setup coordinate arrays
#         self._fault_plane_corners = [[None, None, None], # a 
#                                      [None, None, None], # b
#                                      [None, None, None], # c
#                                      [None, None, None]] # d
#         self._fault_plane_centers = [[None, None, None], # 1
#                                      [None, None, None], # 2 
#                                      [None, None, None]] # 3

#         # Convert values to meters
#         dimensions, depth = self.convert2meters(["dimensions", "slip"])

#         # Depths (in meters)
#         if self.coordinate_specification == 'top center':
#             self._fault_plane_corners[0][2] = depth
#             self._fault_plane_corners[1][2] = depth     \
#                                  + dimensions[1] * numpy.sin(self.dip * DEG2RAD)
#             self._fault_plane_centers[1][2] = depth      \
#                            + 0.5 * dimensions[1] * numpy.sin(self.dip * DEG2RAD)
#         elif self.coordinate_specification == 'centroid':
#             self._fault_plane_corners[0][2] = depth     \
#                                  - dimensions[1] * numpy.sin(self.dip * DEG2RAD)
#             self._fault_plane_corners[1][2] = depth     \
#                                  + dimensions[1] * numpy.sin(self.dip * DEG2RAD)
#             self._fault_plane_centers[1][2] = depth
#         elif self.coordinate_specification == 'bottom center':
#             self._fault_plane_corners[0][2] = depth
#             self._fault_plane_corners[1][2] = depth     \
#                                  - dimensions[1] * numpy.sin(self.dip * DEG2RAD)
#             self._fault_plane_centers[1][2] = depth      \
#                            - 0.5 * dimensions[1] * numpy.sin(self.dip * DEG2RAD)
#         self._fault_plane_corners[2][2] = self._fault_plane_corners[0][2]
#         self._fault_plane_corners[3][2] = self._fault_plane_corners[1][2]
#         self._fault_plane_centers[0][2] = self._fault_plane_corners[0][2]
#         self._fault_plane_centers[2][2] = self._fault_plane_corners[1][2]

#         # Convert dimensions to lat-long
#         dimensions[0] *= 1.0 / lat2meter
#         dimensions[1] *= 1.0 / lat2meter

#         # Calculate xy-plane projected width
#         xy_width = dimensions[1] * numpy.cos(self.dip * DEG2RAD)
        
#         # Locate fault plane in 3D space
#         # Note that the coodinate specification is in reference to the fault 
#         #
#         #    Top edge    Bottom edge
#         #      a ----------- b
#         #      |             |
#         #      |             |
#         #      |             |
#         #      |             |
#         #      1      2      3   L
#         #      |             |
#         #      |             |
#         #      |             |
#         #      |             |
#         #      d ----------- c
#         #            W
#         #
#         xy_corners = [None, None, None, None]
#         xy_centers = [None, None, None]
#         if self.coordinate_specification == 'top center':
#             # Non-rotated locations of corners
#             xy_corners[0] = (self.coordinates[0],
#                              self.coordinates[1] + 0.5 * dimensions[0])
#             xy_corners[1] = (self.coordinates[0] + xy_width,
#                              self.coordinates[1] + 0.5 * dimensions[0])
#             xy_corners[2] = (self.coordinates[0] + xy_width,
#                              self.coordinates[1] - 0.5 * dimensions[0])
#             xy_corners[3] = (self.coordinates[0],
#                              self.coordinates[1] - 0.5 * dimensions[0])

#             # Non-rotated lcoations of center-line coordinates
#             xy_centers[0] = self.coordinates
#             xy_centers[1] = (self.coordinates[0] + 0.5 * xy_width,
#                              self.coordinates[1])
#             xy_centers[2] = (self.coordinates[0] + xy_width,
#                              self.coordinates[1])

#         elif self.coordinate_specification == 'centroid':
#             # Non-rotated locations of corners
#             xy_corners[0] = (self.coordinates[0] - 0.5 * xy_width,
#                              self.coordinates[1] + 0.5 * dimensions[0])
#             xy_corners[1] = (self.coordinates[0] + 0.5 * xy_width,
#                              self.coordinates[1] + 0.5 * dimensions[0])
#             xy_corners[2] = (self.coordinates[0] + 0.5 * xy_width,
#                              self.coordinates[1] - 0.5 * dimensions[0])
#             xy_corners[3] = (self.coordinates[0] - 0.5 * xy_width,
#                              self.coordinates[1] - 0.5 * dimensions[0])

#             # Non-rotated lcoations of center-line coordinates
#             xy_centers[0] = (self.coordinates[0] - 0.5 * xy_width,
#                              self.coordinates[1])
#             xy_centers[1] = self.coordinates
#             xy_centers[2] = (self.coordinates[0] + 0.5 * xy_width,
#                              self.coordinates[1])

#         elif self.coordinate_specification == 'bottom center':
#             # Non-rotated locations of corners
#             xy_corners[0] = (self.coordinates[0] - xy_width,
#                              self.coordinates[1] + 0.5 * dimensions[0])
#             xy_corners[1] = (self.coordinates[0],
#                              self.coordinates[1] + 0.5 * dimensions[0])
#             xy_corners[2] = (self.coordinates[0],
#                              self.coordinates[1] - 0.5 * dimensions[0])
#             xy_corners[3] = (self.coordinates[0] - xy_width,
#                              self.coordinates[1] - 0.5 * dimensions[0])

#             # Non-rotated lcoations of center-line coordinates
#             xy_centers[0] = (self.coordinates[0] - xy_width,
#                              self.coordinates[1])
#             xy_centers[1] = (self.coordinates[0] - 0.5 * xy_width,
#                              self.coordinates[1])
#             xy_centers[2] = self.coordinates

#         else:
#             raise ValueError("Unknown coordinate specification '%s'." % self.coordinate_specification)

#         # Rotate fault plane corners and store
#         # top_center = [xy_corners[0][0], xy_corners[0][1] - 0.5 * dimensions[0]]
#         for (n, corner) in enumerate(xy_corners):
#             self._fault_plane_corners[n][0:2] = \
#                        self.transform(self.coordinates, corner, self.strike * DEG2RAD)
#         for (n, center) in enumerate(xy_centers):
#             self._fault_plane_centers[n][0:2] = \
#                        self.transform(self.coordinates, center, self.strike * DEG2RAD)


#     def create_coordinate_arrays(self, resolution=60, buffer_size=0.5):
#         r"""Create coordinate arrays containing subfault.

#         Input
#         -----
#          - *resolution* (int) - Number of grid points per degree.  Defaults to
#            1" resolution.
#          - *buffer_size* (float) - Buffer distance around edge of fault in 
#            degrees, defaults to 0.5 degrees.

#         """

#         rect = self.containing_rect()
#         rect[0] -= buffer_size
#         rect[1] += buffer_size
#         rect[2] -= buffer_size
#         rect[3] += buffer_size

#         self.delta = float(1.0 / resolution)
#         N = [int((rect[1] - rect[0]) * resolution),
#              int((rect[3] - rect[2]) * resolution)]

#         self._x = numpy.linspace(rect[0],rect[1],N[0])
#         self._y = numpy.linspace(rect[2],rect[3],N[1])


#     def create_deformation_array(self, domain=None):
#         r"""Create deformation array dZ.

#         Use Okada model to calculate deformation from subfault parameters 
#         contained in this object.

#         Currently only calculates the vertical displacement.

#         """
#         dimensions, depth, slip = self.convert2meters(["dimensions", "depth", "slip"])

#         # Construct dictionary that okadamap is looking for
#         okada_params = {}
#         okada_params["depth"] = depth
#         okada_params["length"] = dimensions[0]
#         okada_params["width"] = dimensions[1]
#         okada_params["slip"] = slip
#         okada_params["strike"] = self.strike
#         okada_params["dip"] = self.dip
#         okada_params["rake"] = self.rake
#         if self.coordinate_specification == 'bottom center':
#             # okada_map does not support the bottom center specification
#             okada_params["longitude"] = self.fault_plane_centers[1][0]
#             okada_params["latitude"] = self.fault_plane_centers[1][1]
#             okada_params["latlong_location"] = 'centroid'
#         else:
#             okada_params["longitude"] = self.coordinates[0]
#             okada_params["latitude"] = self.coordinates[1]
#             okada_params["latlong_location"] = self.coordinate_specification
        
#         self._dZ = okadamap(okada_params, self.x, self.y)

#         # Calculate fault on different domain if requested
#         if domain is not None:
#             return okadamap(okada_params, domain[0], domain[1])
#         else:
#             return self._dZ


#     def plot_rake(self, axes=None, color='r', markerstyle="o", linestyle='-'):
#         r"""Plot fault rectangle.

#         Input
#         -----
#          - *axes* (`matplotlib.axes.Axes`) - 

#         Output
#         ------
#          - (`matplotlib.axes.Axes`) - 

#         """

#         import matplotlib.pyplot as plt
        
#         # Create axes object if needed
#         if axes is None:
#             fig = plt.figure()
#             axes = fig.add_subplot(1,1,1)

#         centroid = self.fault_plane_centers[1][:2]
#         top_edge = self.fault_plane_centers[0][:2]
#         r = numpy.sqrt((top_edge[0] - self.fault_plane_corners[0][0])**2 +
#                        (top_edge[1] - self.fault_plane_corners[0][1])**2 )
#         theta = (self.strike + self.rake) * DEG2RAD
#         xy_rake = (r * numpy.cos(-theta + 1.5 * numpy.pi) + centroid[0], 
#                    r * numpy.sin(-theta + 1.5 * numpy.pi) + centroid[1])

#         axes.annotate("",
#             xy=xy_rake, xycoords='data',
#             xytext=self.fault_plane_centers[1][:2], textcoords='data',
#             arrowprops=dict(arrowstyle="->", connectionstyle="arc3") )

#         return axes


#     def plot_subfault_depth(self, axes=None):
#         r"""Plot the depth of each subfault vs. x in one plot and vs. y in a second plot.
        
#         """

#         raise NotImplemented("Subfault depth plots not implemented.")


# class UCSBFault(Fault):

#     r"""Fault subclass for reading in subfault format models from UCSB

#     Read in subfault format models produced by Chen Ji's group at UCSB,
#     downloadable from:  

#         http://www.geol.ucsb.edu/faculty/ji/big_earthquakes/home.html

#     """

#     # @property
#     # def slip(self):
#     #     r"""Array of subfault slips

#     #     :TODO:
#     #      - This function needs to be checked to see if it works still
#     #     """
#     #     if self._slip is None and       \
#     #        self.num_cells[0] is not None and self.num_cells[1] is not None:
#     #         self._slip = numpy.ndarray(self.num_cells)
#     #         for (n,subfault) in enumerate(self.subfaults):
#     #             slip = subfault.convert2meters(['slip'])
#     #             # We assume here that the units should be in cm
#     #             j = (n+1) % self.num_cells[1]
#     #             i = int(j / (n+1))
#     #             self._slip[i,j] = slip / 100.0 
#     #     return self._slip
#     # @slip.setter
#     # def slip(self, value):
#     #     self._slip = value
#     # @slip.deleter
#     # def slip(self):
#     #     del self._slip


#     def __init__(self, path=None):
#         r"""UCSBFault initialization routine.
        
#         See :class:`UCSBFault` for more info.

#         """

#         self.num_cells = [None, None]

#         super(UCSBFault, self).__init__(path=path)


#     def read(self, path):
#         r"""Read in subfault specification at *path*.

#         Creates a list of subfaults from the subfault specification file at
#         *path*.

#         """

#         # Read header of file
#         regexp_dx = re.compile(r"Dx=[ ]*(?P<dx>[^k]*)")
#         regexp_dy = re.compile(r"Dy=[ ]*(?P<dy>[^k]*)")
#         regexp_nx = re.compile(r"nx[^=]*=[ ]*(?P<nx>[^D]*)")
#         regexp_ny = re.compile(r"ny[^=]*=[ ]*(?P<ny>[^D]*)")
#         found_subfault_discretization = False
#         found_subfault_boundary = False
#         header_lines = 0
#         with open(path, 'r') as subfault_file:
#             # Find fault secgment discretization
#             for (n,line) in enumerate(subfault_file):
#                 result_dx = regexp_dx.search(line)
#                 result_dy = regexp_dy.search(line)
#                 result_nx = regexp_nx.search(line)
#                 result_ny = regexp_ny.search(line)

#                 if result_dx and result_dy:
#                     dx = float(result_dx.group('dx'))
#                     dy = float(result_dy.group('dy'))
#                     self.num_cells[0] = int(result_nx.group('nx'))
#                     self.num_cells[1] = int(result_ny.group('ny'))
#                     found_subfault_discretization = True
#                     break
#             header_lines += n

#             # Parse boundary
#             in_boundary_block = False
#             boundary_data = []
#             for (n,line) in enumerate(subfault_file):
#                 if line[0].strip() == "#":
#                     if in_boundary_block and len(boundary_data) == 5:
#                         found_subfault_boundary = True
#                         break
#                 else:
#                     in_boundary_block = True
#                     boundary_data.append([float(value) for value in line.split()])

#         # Assume that there is a column label right underneath the boundary
#         # specification
#         header_lines += n + 2

#         # Check to make sure last boundary point matches, then throw away
#         if boundary_data[0] != boundary_data[4]:
#             raise ValueError("Boundary specified incomplete: ",
#                              "%s" % boundary_data)

#         # Locate fault plane in 3D space - see SubFault `calculate_geometry` for
#         # a schematic of where these points are
#         self._fault_plane_corners = [None, # a 
#                                      None, # b
#                                      None, # c
#                                      None] # d
#         self._fault_plane_centers = [[0.0, 0.0, 0.0], # 1
#                                      [0.0, 0.0, 0.0], # 2 
#                                      [0.0, 0.0, 0.0]] # 3
#         # :TODO: Is the order of this a good assumption?
#         self._fault_plane_corners[0] = boundary_data[0]
#         self._fault_plane_corners[3] = boundary_data[1]
#         self._fault_plane_corners[2] = boundary_data[2]
#         self._fault_plane_corners[1] = boundary_data[3]
        
#         # Calculate center by averaging position of appropriate corners
#         for (n, corner) in enumerate(self._fault_plane_corners):
#             for i in xrange(3):
#                 self._fault_plane_centers[1][i] += corner[i] / 4
#             if n == 0 or n == 4:
#                 for i in xrange(3):
#                     self._fault_plane_centers[0][i] += corner[i] / 2
#             else:
#                 for i in xrange(3):
#                     self._fault_plane_centers[2][i] += corner[i] / 2


#         if not (found_subfault_boundary and found_subfault_discretization):
#             raise ValueError("Could not find base fault characteristics in ",
#                              "subfault specification file at %s." % path)

#         # Calculate center of fault

#         column_map = {"coordinates":(1,0), "depth":2, "slip":3, "rake":4, 
#                       "strike":5, "dip":6, 'mu':10}

#         super(UCSBFault, self).read(path, column_map, skiprows=header_lines)

#         # Set general data
#         for subfault in self.subfaults:
#             subfault.dimensions = [dx, dy]
#             subfault.units.update({"slip":"cm", "depth":"km", 'mu':"dyne/cm^2",
#                                    "dimensions":"km"})


#     def plot_slip(self, axes=None):
#         r"""Plot each subfault with color based on slip magnitude."""
#         raise NotImplemented("")



# class CSVFault(Fault):

#     r"""Fault subclass for reading in CSV formatted files

#     Assumes that the first row gives the column headings
#     """

#     def read(self, path):
#         r"""Read in subfault specification at *path*.

#         Creates a list of subfaults from the subfault specification file at
#         *path*.

#         """

#         # Read header of file
#         with open(path, 'r') as subfault_file:
#             header_line = subfault_file.readline().split(",")
#             column_map = {'coordinates':[None, None], 'dimensions':[None, None]}
#             for (n,column_heading) in enumerate(header_line):
#                 if "(" in column_heading:
#                     # Strip out units if present
#                     unit_start = column_heading.find("(")
#                     unit_end = column_heading.find(")")
#                     column_key = column_heading[:unit_start].lower()
#                     self.units[column_key] = column_heading[unit_start:unit_end]
#                 else:
#                     column_key = column_heading.lower()
#                 column_key = column_key.strip()
#                 if column_key in ("depth", "strike", "dip", "rake", "slip"):
#                     column_map[column_key] = n
#                 elif column_key == "longitude":
#                     column_map['coordinates'][0] = n
#                 elif column_key == "latitude":
#                     column_map['coordinates'][1] = n
#                 elif column_key == "length":
#                     column_map['dimensions'][0] = n
#                 elif column_key == "width":
#                     column_map['dimensions'][1] = n

#         super(CSVFault, self).read(path, column_map=column_map, skiprows=1,
#                                          delimiter=",")

# # path, column_map, coordinate_specification="centroid",
# #                          rupture_type="static", t=[0.0, 0.5, 1.0], skiprows=0


# Alternate version of Okada tools
# alpha = poisson

# def calculate_deformations(fault_data, x, y):
#     r""""""

#     delta = fault_data.dip_angle * numpy.pi / 180.0
#     L = fault_data.length
#     w = fault_data.width
#     U = fault_data.dislocation

#     # Construct grid
#     X, Y = numpy.meshgrid(x,y)
#     u = numpy.empty((3,y.shape[0],x.shape[0]))

#     # Calculate relevant quantities
#     d = fault_data.depth
#     p = Y * cos(delta) + d * sin(delta)
#     q = Y * sin(delta) - d * cos(delta)

#     # u_B_1 =  dip_f_b(X + L, p, q, delta, i=1)        \
#     #        - dip_f_b(X + L, p-w, q, delta, i=1)      \
#     #        - dip_f_b(X - L, p, q, delta, i=1)        \
#     #        + dip_f_b(X - L, p-w, q, delta, i=1)
#     # u_B_2 =  dip_f_b(X + L, p, q, delta, i=2)        \
#     #        - dip_f_b(X + L, p-w, q, delta, i=2)      \
#     #        - dip_f_b(X - L, p, q, delta, i=2)        \
#     #        + dip_f_b(X - L, p-w, q, delta, i=2)
#     # u_B_3 =  dip_f_b(X + L, p, q, delta, i=3)        \
#     #        - dip_f_b(X + L, p-w, q, delta, i=3)      \
#     #        - dip_f_b(X - L, p, q, delta, i=3)        \
#     #        + dip_f_b(X - L, p-w, q, delta, i=3)
#     # u[0,:,:] = U / (2.0 * numpy.pi) * (u_B_1)
#     # u[1,:,:] = U / (2.0 * numpy.pi)                 \
#     #                     * (u_B_2 * cos(delta) + u_B_3 * sin(delta))
#     # u[2,:,:] = U / (2.0 * numpy.pi)                 \
#     #                     * (u_B_2 * sin(delta) - u_B_3 * cos(delta))

#     u_B_1 =   dip_f_b(    X,     p, q, delta, i=1)     \
#             - dip_f_b(    X, p - w, q, delta, i=1)     \
#             - dip_f_b(X - L,     p, q, delta, i=1)     \
#             + dip_f_b(X - L, p - w, q, delta, i=1)

#     u_B_2 =   dip_f_b(    X,     p, q, delta, i=2)     \
#             - dip_f_b(    X, p - w, q, delta, i=2)     \
#             - dip_f_b(X - L,     p, q, delta, i=2)     \
#             + dip_f_b(X - L, p - w, q, delta, i=2)

#     u_B_3 =   dip_f_b(    X,     p, q, delta, i=3)     \
#             - dip_f_b(    X, p - w, q, delta, i=3)     \
#             - dip_f_b(X - L,     p, q, delta, i=3)     \
#             + dip_f_b(X - L, p - w, q, delta, i=3)

#     u[0,:,:] = U / (2.0 * numpy.pi) * u_B_1
#     u[1,:,:] = U / (2.0 * numpy.pi) * (u_B_2 * cos(delta) - u_B_3 * sin(delta))
#     u[2,:,:] = U / (2.0 * numpy.pi) * (u_B_2 * sin(delta) + u_B_3 * cos(delta))

#     return u


# def dip_f_b(xi, eta, q, delta, i=2):
#     r"""Calculate the dip component of f^B_i

#     """

#     import numpy

#     d_bar = eta * sin(delta) - q * cos(delta)
#     R = numpy.sqrt(xi**2 + eta**2 + q**2)

#     if i == 1:
#         y_bar = eta * cos(delta) + q * sin(delta)
        
#         I_3 = 1.0 / cos(delta) * y_bar / (R + d_bar)                     \
#             - 1.0 / cos(delta)**2                                        \
#             * (numpy.log(R + eta) - sin(delta) * numpy.log(R + d_bar))

#         return -q / R + (1.0 + alpha) / alpha * I_3 * sin(delta) * cos(delta)
    
#     elif i == 2:
#         X_11 = 1.0 / (R * (R + xi))
#         theta = numpy.arctan(xi * eta / (q * R))
        
#         return -eta * q * X_11 - theta - (1.0 - alpha) / alpha * xi / (R + d_bar) * sin(delta) * cos(delta)
    
#     elif i == 3:
#         X_11 = 1.0 / (R * (R + xi))
#         X = numpy.sqrt(xi**2 + q**2)
#         I_4 = numpy.tan(delta) * xi / (R + d_bar) - 2.0 / cos(delta)**2 * numpy.arctan((eta * (X + q * cos(delta)) + X * (R + X) * sin(delta)) / (xi * (R + X) * cos(delta)))

#         return q**2 * X_11 + (1.0 - alpha) / alpha * I_4 * sin(delta) * cos(delta)
#     else:
#         raise ValueError("Invalid component %s requested.  Valid range = [1,3]" % i)

# ==============================================================================
