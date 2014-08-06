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
#  General utility functions
# ==============================================================================
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
                             verbose=False, fig_kwargs={}):
    r"""
    Plot sea floor deformation dz as colormap with contours
    """

    from clawpack.visclaw import colormaps
    import matplotlib.pyplot as plt

    if axes is None:
        fig = plt.figure(**fig_kwargs)
        axes = fig.add_subplot(1, 1, 1)

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
    axes.set_xticklabels([label.set_rotation(80) 
                                           for label in axes.get_xticklabels()])
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


    def read(self, path=None, dtopo_type=None):
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


    def animate_dz_colors(self, times=None, cmax_dz=None, dz_interval=None,
                                style='loop', axes=None, fig_kwargs={}):
        """
        Animate seafloor motion for time-dependent ruptures.
        Interpolate dz_list to specified times and then call module function
        plot_dz_colors.
        
        *style* == 'loop'  ==> plot in a simple loops, works in IPython shell
        *style* == 'html'  ==> use JSAnimation and embed in html file
        *style* == 'notebook'  ==> use JSAnimation and display in notebook
        """

        from time import sleep
        import matplotlib.pyplot as plt

        # Create axes if needed
        if axes is None:
            fig = plt.figure(**fig_kwargs)
            axes = fig.add_subplot(111)

        
        if style in ['html', 'notebook']:
            from clawpack.visclaw.JSAnimation import IPython_display
            import clawpack.visclaw.JSAnimation.JSAnimation_frametools as J
            plotdir = '_plots'
            J.make_plotdir(plotdir, clobber=True)

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
            axes = plot_dz_colors(self.X,self.Y,dz_list[k],axes=axes,
                        cmax_dz=cmax_dz, dz_interval=dz_interval, 
                        fig_kwargs=fig_kwargs)
            plt.title("Seafloor deformation at t = %12.6f" % t)
            plt.draw()
            if style == 'loop':
                sleep(0.5)
            else:
                J.save_frame(k, verbose=True)

        if style in ['html', 'notebook']:
            anim = J.make_anim(plotdir)
        if style == 'html':
            file_name="dtopo_movie.html"
            J.make_html(anim, file_name=file_name, title="Moving dtopo")
        if style == 'notebook':
            return anim

                

# ==============================================================================
#  Generic Fault Class
# ==============================================================================
class Fault(object):

    def __init__(self, subfaults=None, units={}):

        # Parameters for subfault specification
        self.rupture_type = 'static' # 'static' or 'dynamic'
        #self.times = numpy.array([0., 1.])   # or just [0.] ??
        self.dtopo = None

        # Default units of each parameter type
        self.units = {}
        self.units.update(units)
        
        if subfaults is not None:
            if not isinstance(subfaults, list):
                raise ValueError("Input parameter subfaults must be a list.")
            self.subfaults = subfaults


    def read(self, path, column_map, coordinate_specification="centroid",
                                     rupture_type="static", skiprows=0, 
                                     delimiter=None, units={}, defaults=None):
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
        r"""Calculate the moment magnitude for a fault composed of subfaults."""
        return Mw(self.Mo())

    
    def create_dtopography(self, x, y, times=[0., 1.], verbose=False):
        r"""Compute change in topography and construct a dtopography object.

        Use subfaults' `okada` routine and add all 
        deformations together.

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
            raise Exception("Unrecognized rupture_type: %s" % self.rupture_type)

        # Store for user
        self.dtopo = dtopo

        return dtopo


    
    def plot_subfaults(self, axes=None, plot_centerline=False, slip_color=False,
                             cmap_slip=None, cmin_slip=None, cmax_slip=None,
                             plot_rake=False, xylim=None, plot_box=True):
        """
        Plot each subfault projected onto the surface.
        Describe parameters...
        """
    
        import matplotlib
        import matplotlib.pyplot as plt
    
        if axes is None:
            fig = plt.figure()
            axes = fig.add_subplot(1, 1, 1)
        # for testing purposes, make random slips:
        test_random = False
    
        max_slip = 0.
        min_slip = 0.
        for subfault in self.subfaults:
            if test_random:
                subfault.slip = 10.*numpy.rand()  # for testing
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
            if test_random:
                print "*** test_random == True so slip and rake have been randomized"
            
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
                if test_random:
                    subfault.rake = 90. + 30.*(rand()-0.5)  # for testing
                tau = (subfault.rake - 90) * numpy.pi/180.
                axes.plot([x_centroid],[y_centroid],'go',markersize=5,label="Centroid")
                dxr = x_top - x_centroid
                dyr = y_top - y_centroid
                x_rake = x_centroid + numpy.cos(tau)*dxr - numpy.sin(tau)*dyr
                y_rake = y_centroid + numpy.sin(tau)*dxr + numpy.cos(tau)*dyr
                axes.plot([x_rake,x_centroid],[y_rake,y_centroid],'g-',linewidth=1)
            if slip_color:
                slip = subfault.slip
                s = min(1, max(0, (slip-cmin_slip)/(cmax_slip-cmin_slip)))
                c = cmap_slip(s*.99)  # since 1 does not map properly with jet
                axes.fill(x_corners,y_corners,color=c,edgecolor='none')
            if plot_box:
                axes.plot(x_corners, y_corners, 'k-')
    
        slipax = axes
            
        y_ave = y_ave / len(self.subfaults)
        slipax.set_aspect(1./numpy.cos(y_ave*numpy.pi/180.))
        axes.ticklabel_format(format='plain',useOffset=False)

        ## RJL: setting labels like this gives None's as labels:
        #axes.set_xticklabels([label.set_rotation(80) 
        #                                   for label in axes.get_xticklabels()])
        if xylim is not None:
            axes.set_xlim(xylim[:2])
            axes.set_ylim(xylim[2:])
        axes.set_title('Fault planes')
        if slip_color:
            cax,kw = matplotlib.colorbar.make_axes(slipax)
            norm = matplotlib.colors.Normalize(vmin=cmin_slip,vmax=cmax_slip)
            cb1 = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap_slip, norm=norm)
        plt.sca(slipax) # reset the current axis to the main figure
    


    def plot_subfaults_depth(self):
        """
        Plot the depth of each subfault vs. x in one plot and vs. y in a second plot.
        """
    
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(nrows=2, ncols=1)
    
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



# ==============================================================================
#  Sub-Fault Class
# ==============================================================================
class SubFault(object):
    r"""Basic sub-fault specification.

    """

    def __init__(self, units={}):
        r"""SubFault initialization routine.
        
        See :class:`SubFault` for more info.

        """
        
        super(SubFault, self).__init__()

        self.strike = None
        r"""Strike direction of subfault in degrees."""
        self.length = None
        r"""Length of subfault in specified units."""
        self.width = None
        r"""Width of subfault in specified units."""
        self.depth = None
        r"""Depth of subfault based on *coordinate_specification*."""
        self.slip = None
        r"""Slip on subfault in strike direction in specified units."""
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
        self.units = {'mu':"dyne/cm^2", 'length':'km', 'width':'km', 
                      'depth':'km'}
        r"""Dictionary of units for the relevant parameters."""

        self.units.update(units)


        self._geometry = None


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
        r"""Calculate the seismic moment for a single subfault

        Returns in units of in units N-m. 
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

        length, width, slip = self.convert2meters(["length", "width", "slip"])
        total_slip = length * width * slip
        Mo = mu * total_slip
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
        length, width, depth, slip = self.convert2meters(["length","width", \
                                    "depth","slip"])

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
        ang_strike = DEG2RAD*th
        halfL = 0.5*L
    
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
        if self.units['depth'] == 'km':
            depth_top = depth_top / 1000.
            depth_bottom = depth_bottom / 1000.

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
        length, width, depth, slip = self.convert2meters(["length","width", \
                                    "depth","slip"])

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

        super(UCSBFault, self).read(path, column_map, skiprows=header_lines,
                                coordinate_specification="centroid")

        # Set general data
        for subfault in self.subfaults:
            subfault.length = dx
            subfault.width = dy
            subfault.units.update({"slip":"cm", "depth":"km", 'mu':"dyne/cm^2",
                                   "length":"km", "width":"km"})



# ==============================================================================
#  CSV sub-class of Fault
# ==============================================================================
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
        units = {'length':'km', 'width':'km', 'depth':'km', 'slip':'m',
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
                subfault.units = units
                # subfault.mu = ??  ## currently using SubFault default
                self.sift_subfaults[name] = subfault

