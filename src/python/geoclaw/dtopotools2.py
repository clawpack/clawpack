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
    
:Functions:


:TODO:
 - List functions and classes in this docstring?
 - This file contains both the older dtopo functionality and a new class, should
   merge as much functionality into the class as possible ensuring nothing is
   left behind.
 - Refactor Okada functionality
"""

import os, sys, copy, re
import numpy
import clawpack.geoclaw.topotools as topotools

# ==============================================================================
# Constants
# ==============================================================================

from data import Rearth
DEG2RAD = numpy.pi / 180.0
RAD2DEG = 180.0 / numpy.pi
lat2meter = Rearth * DEG2RAD
poisson = 0.25  # Poisson ratio for Okada 

# ==============================================================================
#  General utility functions
# ==============================================================================

def plot_dz_contours(x,y,dz,dz_interval=0.5,verbose=False):
    r"""For plotting seafloor deformation dz"""
    import matplotlib.pyplot as plt
    dzmax = max(dz.max(), -dz.min()) + dz_interval
    clines1 = numpy.arange(dz_interval, dzmax, dz_interval)
    clines = list(-numpy.flipud(clines1)) + list(clines1)

    if len(clines) > 0:
        if verbose:
            print "Plotting contour lines at: ",clines
        plt.contour(x,y,dz,clines,colors='k')
    else:   
        print "No contours to plot"

def plot_dz_colors(x,y,dz,cmax_dz=None,dz_interval=None,verbose=False):
    r"""
    Plot sea floor deformation dz as colormap with contours
    """

    from clawpack.visclaw import colormaps
    import matplotlib.pyplot as plt

    dzmax = abs(dz).max()
    if cmax_dz is None:
        if dzmax < 1.e-12:
            cmax_dz = 0.1
        else:
            cmax_dz = dzmax
    cmap = colormaps.blue_white_red
    plt.pcolor(x, y, dz, cmap=cmap)
    plt.clim(-cmax_dz,cmax_dz)
    cb2 = plt.colorbar(shrink=1.0)
    
    if dz_interval is None:
        dz_interval = cmax_dz/10.
    clines1 = numpy.arange(dz_interval, dzmax + dz_interval, dz_interval)
    clines = list(-numpy.flipud(clines1)) + list(clines1)
    if len(clines) > 0:
        if verbose:
            print "Plotting contour lines at: ",clines
        plt.contour(x,y,dz,clines,colors='k',linestyles='solid')
    elif verbose:
        print "No contours to plot"

    y_ave = 0.5*(y.min() + y.max())
    plt.gca().set_aspect(1./numpy.cos(y_ave*numpy.pi/180.))

    plt.ticklabel_format(format='plain',useOffset=False)
    plt.xticks(rotation=80)
    plt.title('Seafloor deformation')


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
    

def rise_fraction(t, t0, t_rise, t_rise_ending=None):
    """
    A continuously differentiable piecewise quadratic function of t that is 
       0 for t <= t0, 
       1 for t >= t0 + t_rise + t_rise_ending 
    with maximum slope at t0 + t_rise.
    For specifying dynamic fault ruptures:  Subfault files often contain these
    parameters for each subfault for an earthquake event.
    """

    t = numpy.array(t)
    if t_rise_ending is None: 
        t_rise_ending = t_rise

    t1 = t0+t_rise
    t2 = t1+t_rise_ending

    rf = numpy.where(t<=t0, 0., 1.)
    if t2==t0:
        return rf

    t20 = float(t2-t0)
    t10 = float(t1-t0)
    t21 = float(t2-t1)

    c1 = t21 / (t20*t10*t21) 
    c2 = t10 / (t20*t10*t21) 

    rf = numpy.where((t>t0) & (t<=t1), c1*(t-t0)**2, rf)
    rf = numpy.where((t>t1) & (t<=t2), 1. - c2*(t-t2)**2, rf)

    return rf


# ==============================================================================
#  Classes: 
# ==============================================================================

class DTopography(object):

    def __init__(self, path=None, dtopo_type=None):
        self.dz_list = []
        self.times = []
        self.x = None
        self.y = None
        self.X = None
        self.Y = None
        self.delta = None
        if path:
            self.read(path, dtopo_type)

    def read(self, path, dtopo_type):
        r"""
        Read in a dtopo file and use to set attributes of this object.

        input
        -----
         - *path* (path) - Path to the output file to written to.
         - *dtopo_type* (int) - Type of topography file to write out.  Default
           is 1.
        """
        if dtopo_type == 2 or dtopo_type == 3:
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

        elif dtopo_type==1:
            d = numpy.loadtxt(path)
            print "Loaded file %s with %s lines" %(path,d.shape[0])
            t = list(set(d[:,0]))
            t.sort()
            print "times found: ",t
            ntimes = len(t)
            tlast = t[-1]
            lastlines = d[d[:,0]==tlast]
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
                dz = numpy.reshape(d[i1:i2,3],(my,mx))
                dz = numpy.flipud(dz)
                dz_list.append(dz)
            self.X = X
            self.Y = Y
            self.x = X[0,:]
            self.y = Y[:,0]
            self.times = t
            self.dz_list = dz_list

        else:
            raise NotImplementedError("*** Not implemented for dtopo_type: %s" % dtopo_type)


    def write(self, path, dtopo_type=None):
        r"""Write out subfault resulting dtopo to file at *path*.

        input
        -----
         - *path* (path) - Path to the output file to written to.
         - *dtopo_type* (int) - Type of topography file to write out.  Default
           is 1.

        """

        if dtopo_type is None:
            # Try to look at suffix for type
            extension = os.path.splitext(path)[1][1:]
            if extension[:2] == "tt":
                dtopo_type = int(extension[2])
            elif extension == 'xyz':
                dtopo_type = 0
            elif extension == 'txyz':
                dtopo_type = 1
            else:
                # Default to 3
                dtopo_type = 3

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
                raise ValueError("Only topography types 1, 2, and 3 are supported.")


    def dz(self,t):
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


    def plot_dz_colors(self,t,cmax_dz=None,dz_interval=None):
        """
        Interpolate dz_list to specified time t and then call module function
        plot_dz_colors.
        """
        plot_dz_colors(self.X,self.Y,self.dz(t),cmax_dz=cmax_dz, \
                       dz_interval=dz_interval)


    def plot_dz_contours(self,t,dz_interval=0.5):
        """
        Interpolate dz_list to specified time t and then call module function
        plot_dz_contours.
        """
        plot_dz_contours(self.X,self.Y,self.dz(t),dz_interval=dz_interval)

    def animate_dz_colors(self,times=None,cmax_dz=None,dz_interval=None):
        """
        Animate seafloor motion for time-dependent ruptures.
        Interpolate dz_list to specified times and then call module function
        plot_dz_colors.
        """
        from time import sleep
        import matplotlib.pyplot as plt

        if times is None:
            times = self.times

        dz_list = []
        max_dz = 0.
        for t in times:
            dz_t = self.dz(t)
            dz_list.append(dz_t)
            max_dz = max(max_dz, abs(dz_t).max())
        if cmax_dz is None:
            cmax_dz = max_dz
        print "max abs(dz(t)) = ",max_dz
        plt.figure()
        for k,t in enumerate(times):
            plt.clf()
            plot_dz_colors(self.X,self.Y,dz_list[k],cmax_dz=cmax_dz, \
                       dz_interval=dz_interval)
            plt.title("Seafloor deformation at t = %12.6f" % t)
            plt.draw()
            sleep(0.5)


class Fault(object):

    def __init__(self, path=None, subfaults=None, units={}):

        # Parameters for subfault specification
        self.rupture_type = 'static' # 'static' or 'dynamic'
        #self.times = numpy.array([0., 1.])   # or just [0.] ??
        self.dtopo = None

        # Default units of each parameter type
        self.units = {}
        self.units.update(units)
        
        if path is not None:
            # Read in file at path assuming it is a subfault specification
            self.read(path)
        elif subfaults is not None:
            if not isinstance(subfaults, list):
                raise ValueError("Input parameter subfaults must be a list.")
            self.subfaults = subfaults


    def read(self, path, column_map, coordinate_specification="centroid",
                         rupture_type="static", skiprows=0, delimiter=None,
                         units={}, defaults=None):
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
            of the subfault.  See the documentation for *SubFault.okada*.
          - *rupture_type* (str) either "static" or "dynamic"
          - *skiprows* (int) number of header lines to skip before data
          - *delimiter* (str) e.g. ',' for csv files
          - *units* (dict) indicating units for length, width, slip, depth,
                           and for rigidity mu.
          - *defaults* (dict) default values for all subfaults, for values not
                       included in subfault file on each line.

        """

        # Read in rest of data
        # (Use genfromtxt to deal with files containing strings, e.g. unit
        # source name, in some column)
        data = numpy.genfromtxt(path, skiprows=skiprows, delimiter=delimiter)
        if len(data.shape) == 1:
            data = numpy.array([data])

        self.subfaults = []
        for n in xrange(data.shape[0]):

            new_subfault = SubFault()
            new_subfault.coordinate_specification = coordinate_specification
            new_subfault.units = units
            #new_subfault.rupture_type = rupture_type

            for (var, column) in column_map.iteritems():
                if isinstance(column, tuple) or isinstance(column, list):
                    setattr(new_subfault, var, [None for k in column])
                    for (k, index) in enumerate(column):
                        getattr(new_subfault, var)[k] = data[n, index]
                else:
                    setattr(new_subfault, var, data[n, column])

            self.subfaults.append(new_subfault)
        if defaults is not None:
            for subfault in self.subfaults:
                pass
            raise NotImplementedError("*** need to add support for defaults")



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
        r""" Calculate the moment magnitude for a fault composed of subfaults."""
        return Mw(self.Mo())

    
    def create_deformation_array(self, x, y, times=[0.,1.]):
        r"""Create deformation array dZ.

        Use subfaults' `okada` routine and add all 
        deformations together.

        """

        dtopo = DTopography()
        dtopo.x = x
        dtopo.y = y
        X,Y = numpy.meshgrid(x,y)
        dtopo.X = X
        dtopo.Y = Y
        dtopo.times = times

        print "Making Okada dz for each of %s subfaults" \
                % len(self.subfaults)

        for k,subfault in enumerate(self.subfaults):
                sys.stdout.write("%s.." % k)
                sys.stdout.flush()
                subfault.okada(x,y)  # sets subfault.dtopo with times=[0]
                                     # and len(subfault.dtopo.dz_list) == 1
        sys.stdout.write("\nDone\n")

        if self.rupture_type == 'static':
            if len(times) > 2:
                raise ValueError("For static deformation, need len(times) <= 2")
            dz = numpy.zeros(X.shape)
            for subfault in self.subfaults:
                dz += subfault.dtopo.dz_list[0]

            if len(times)==1:
                dtopo.dz_list = [dz]   # only final deformation stored
            elif len(times)==2:
                dz0 = numpy.zeros(X.shape)
                dtopo.dz_list = [dz0, dz]
            self.dtopo = dtopo
            return dtopo

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
            self.dtopo = dtopo
            return dtopo

        else:   
            raise Exception("Unrecognized rupture_type: %s" % self.rupture_type)



    def plot_subfaults(self, plot_centerline=False, slip_color=False, \
            cmap_slip=None, cmin_slip=None, cmax_slip=None, \
            plot_rake=False, xylim=None, plot_box=True):

        r"""Plot subfault rectangle and slip for all subfaults"""
        raise NotImplementedError("To appear.")


    def plot_subfaults_depth(self):
        """
        Plot the depth of each subfault vs. x in one plot and vs. y in a second plot.
        """
        raise NotImplementedError("To appear.")


class SubFault(object):

    def __init__(self):
        self.strike = None
        self.length = None
        self.width = None
        self.depth = None
        self.slip = None
        self.rake = None
        self.dip = None
        self.latitude = None
        self.longitude = None
        self.coordinate_specification = "top center"
        self.mu = 4e11  # default value for rigidity = shear modulus
        self.units = {'mu':"dyne/cm^2"}


    def convert2meters(self, parameters): 
        r"""Convert relevant lengths to correct units.

        Returns converted (length, width, depth, slip) 

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
        r""" Calculate the seismic moment for a single subfault, in units N-m.
"""

        # Convert units of rigidity mu to Pascals 
        # (1 Pa = 1 Newton / m^2 = 10 dyne / cm^2)
        if self.units["mu"] == "Pa":
            mu = self.mu
        if self.units["mu"] == "GPa":
            # e.g. mu = 40 GPa
            mu = 1e9 * self.mu
        if self.units["mu"] == "dyne/cm^2":
            # e.g. mu = 4e11 dyne/cm^2
            mu = 0.1 * self.mu
        elif self.units["mu"] == 'dyne/m^2':
            # Does anyone use this unit?  
            mu = 1e-5 * self.mu 
        else:
            raise ValueError("Unknown unit for rigidity %s." % self.units['mu'])

        length, width, slip = self.convert2meters(["length","width","slip"])
        total_slip = length * width * slip
        Mo = mu * total_slip
        return Mo


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
        
        Some routines adapted from fortran routines written by
        Xiaoming Wang.
        
        okadamap function riginally written in Python by Dave George for
        Clawpack 4.6 okada.py routine.
        Rewritten and made more flexible by Randy LeVeque:

        Note: *self.coordinate_specification* (str) specifies the location on each
            subfault that corresponds to the (longitude,latitude) and depth 
            of the subfault.  See the documentation for *SubFault.okada*.
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
        length, width, depth, slip = self.convert2meters(["length","width", \
                                    "depth","slip"])

        hh =  depth
        L  =  length
        w  =  width
        d  =  slip
        th =  self.strike
        dl =  self.dip
        rd =  self.rake
        x0 =  self.longitude
        y0 =  self.latitude
        location =  self.coordinate_specification

    
        ang_dip = DEG2RAD*dl
        ang_slip = DEG2RAD*rd
        ang_strike = DEG2RAD*th
        halfL = 0.5*L
    
    
        if location == "bottom center":
            del_x = 0.
            del_y = 0.
            depth_bottom = hh
            x_bottom = x0
            y_bottom = y0

        elif location == "top center":
    
            # Convert focal depth used for Okada's model
            # from top of fault plane to bottom:
    
            depth_top = hh
            hh = hh + w*numpy.sin(ang_dip)
            depth_bottom = hh
    
            # Convert fault origin from top of fault plane to bottom:
            del_x = w*numpy.cos(ang_dip)*numpy.cos(ang_strike) / \
                    (lat2meter*numpy.cos(y0*DEG2RAD))
            del_y = -w*numpy.cos(ang_dip)*numpy.sin(ang_strike) / lat2meter
            
            x_top = x0
            y_top = y0 
            x_bottom = x0+del_x
            y_bottom = y0+del_y
            x_centroid = x0+0.5*del_x
            y_centroid = y0+0.5*del_y
    
    
        elif location == "centroid":
    
            # Convert focal depth used for Okada's model
            # from middle of fault plane to bottom:
            depth_top = hh - 0.5*w*numpy.sin(ang_dip)
            hh = hh + 0.5*w*numpy.sin(ang_dip)
            depth_bottom = hh
    
            # Convert fault origin from middle of fault plane to bottom:
            del_x = 0.5*w*numpy.cos(ang_dip)*numpy.cos(ang_strike) / \
                    (lat2meter*numpy.cos(y0*DEG2RAD))
            del_y = -0.5*w*numpy.cos(ang_dip)*numpy.sin(ang_strike) / lat2meter
    
            x_centroid = x0
            y_centroid = y0
            x_top = x0-del_x
            y_top = y0-del_y
            x_bottom = x0+del_x
            y_bottom = y0+del_y
    
        if location == "noaa sift":
            depth_top = hh
            hh = hh + w*numpy.sin(ang_dip)
            depth_bottom = hh

            del_x = 0.
            del_y = 0.
            depth_bottom = hh
            x_bottom = x0
            y_bottom = y0

        else:
            raise ValueError("Unrecognized coordinate_specification" \
                    % coordinate_specification)
    
        # adjust x0,y0 to bottom center of fault plane:
        x0 = x0 + del_x
        y0 = y0 + del_y
    
        # distance along strike from center of an edge to corner:
        dx2 = 0.5*L*numpy.sin(ang_strike) / \
                (lat2meter*numpy.cos(y_bottom*DEG2RAD))
        dy2 = 0.5*L*numpy.cos(ang_strike) / lat2meter
    
        X,Y = numpy.meshgrid(x,y)   # use convention of upper case for 2d
    
        # Convert distance from (X,Y) to (x_bottom,y_bottom) from degrees to
        # meters:
        xx = lat2meter*numpy.cos(DEG2RAD*Y)*(X-x_bottom)   
        yy = lat2meter*(Y-y_bottom)
    
    
        # Convert to distance along strike (x1) and dip (x2):
        x1 = xx*numpy.sin(ang_strike) + yy*numpy.cos(ang_strike) 
        x2 = xx*numpy.cos(ang_strike) - yy*numpy.sin(ang_strike) 
    
        # In Okada's paper, x2 is distance up the fault plane, not down dip:
        x2 = -x2
    
        p = x2*numpy.cos(ang_dip) + hh*numpy.sin(ang_dip)
        q = x2*numpy.sin(ang_dip) - hh*numpy.cos(ang_dip)
    
        f1=self._strike_slip (x1+halfL,p,  ang_dip,q)
        f2=self._strike_slip (x1+halfL,p-w,ang_dip,q)
        f3=self._strike_slip (x1-halfL,p,  ang_dip,q)
        f4=self._strike_slip (x1-halfL,p-w,ang_dip,q)
    
        g1=self._dip_slip (x1+halfL,p,  ang_dip,q)
        g2=self._dip_slip (x1+halfL,p-w,ang_dip,q)
        g3=self._dip_slip (x1-halfL,p,  ang_dip,q)
        g4=self._dip_slip (x1-halfL,p-w,ang_dip,q)
    
        # Displacement in direction of strike and dip:
        ds = d*numpy.cos(ang_slip)
        dd = d*numpy.sin(ang_slip)
    
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
    


class UCSBFault(Fault):

    r"""Fault subclass for reading in subfault format models from UCSB

    Read in subfault format models produced by Chen Ji's group at UCSB,
    downloadable from:  

        http://www.geol.ucsb.edu/faculty/ji/big_earthquakes/home.html

    """

    def __init__(self, path=None, rupture_type='static'):
        r"""UCSBFault initialization routine.
        
        See :class:`UCSBFault` for more info.
        Subfault format contains info for dynamic rupture, so can specify 
        rupture_type = 'static' or 'dynamic'

        """

        self.num_cells = [None, None]

        super(UCSBFault, self).__init__(path=path)
        self.rupture_type = rupture_type


    def read(self, path):
        r"""Read in subfault specification at *path*.

        Creates a list of subfaults from the subfault specification file at
        *path*.

        """

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

        super(UCSBFault, self).read(path, column_map, skiprows=header_lines,
                                coordinate_specification="centroid")

        # Set general data
        for subfault in self.subfaults:
            subfault.length = dx
            subfault.width = dy
            subfault.units.update({"slip":"cm", "depth":"km", 'mu':"dyne/cm^2",
                                   "length":"km", "width":"km"})


class CSVFault(Fault):

    r"""Fault subclass for reading in CSV formatted files

    Assumes that the first row gives the column headings
    """

    def read(self, path, units={}, coordinate_specification="top center",
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
                    self.units[column_name] = column_heading[unit_start:unit_end]
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
                                delimiter=",", units=units,
                                coordinate_specification=coordinate_specification,
                                rupture_type=rupture_type)
