# Need to fix matplotlib imports

#!/usr/bin/env python
# encoding: utf-8

r"""GeoClaw Moving Topography Tools Module

Module provides several functions for dealing with moving topography (usually
due to earthquakes) including reading sub-fault specifications, writing out 
dtopo files, and calculating Okada based deformations.

:Functions:

:MovingTopography Class:

:TODO:
 - This file contains both the older dtopo functionality and a new class, should
   merge as much functionality into the class as possible ensuring nothing is
   left behind.
 - Refactor Okada functionality
"""

import os
import sys
import copy
import re

import numpy
import numpy as np   # Need to unify: numpy and np both used below
from numpy import sin,cos,sqrt

import clawpack.geoclaw.topotools as topotools

# Constants
from data import Rearth
DEG2RAD = numpy.pi / 180.0
RAD2DEG = 180.0 / numpy.pi
lat2meter = Rearth * DEG2RAD

# ==============================================================================
#  Older Stuff
# ==============================================================================
class DTopo(object):

    def __init__(self, dtopo_params={}):
        self.dtopo_params = dtopo_params
        self.subfaults = []
        self.dz_list = []
        self.times = []
        self.x = None
        self.y = None


def read_subfault_model(fname, columns, units=None, \
                    defaults = {'latlong_location': 'centroid'}, \
                    skiprows=0, delimiter=None):
    """
    Read a subfault model and return a list of dictionaries specifying the
    Okada parameters for each subfault.
    The dictionary may also contain entries 'rupture_time', 'rise_time' and
    perhaps 'rise_time_ending' if these are specified in the file.

    The file *fname* is read in using loadtxt.  The first *skiprow* rows may be
    comments, after that there should be a row for each subfault.  The
    contents of the columns is specified by the input parameter *columns*, 
    each of whose elements is one of the following:
        'latitude', 'longitude', 'length', 'width', 'depth', 
        'strike', 'dip', 'rake', 'slip', 
        'rupture_time', 'rise_time', 'rise_time_ending', 'ignore'
    Columns labelled 'ignore' will be ignored, others will be used to set
    the corresponding elements of the parameter dictionary for each subfault.

    If some Okada parameters are missing from the columns but have the same
    value for each subfault (e.g. 'latlong_location', 'length' or 'width'), 
    set these in the *defaults* dictionary, e.g. 
    defaults = {'latlong_location': 'centroid', 'length': 100e3, 'width': 50e3}
    
    A dictionary *units* can also be provided.  By default the units for
    param = 'length', 'depth', 'width', and 'slip' are assumed to be meters, 
    but you can set units[param] = 'm' or 'km' or 'cm'.

    *delimiter* is the delimiter separating columns.  By default it is any
    whitespace but it could be ',' for a csv file, for example.

    """

    valid_labels = """latitude longitude strike dip rake slip length width
                      depth rupture_time rise_time rise_time_ending ignore""".split()

    if units is None:
        units = {}
        units['slip'] = 'm'
        units['length'] = 'm'
        units['width'] = 'm'
        units['depth'] = 'm'

    usecols = []
    for j,label in enumerate(columns):
        if label not in valid_labels:
            raise Exception("Unrecognized label in columns list: %s" % label)
        if label != 'ignore':
            usecols.append(j)

    try:
        data = numpy.genfromtxt(fname, skiprows=skiprows, delimiter=delimiter,\
                                usecols=usecols)
    except:
        import os
        raise Exception("Unable to load file %s" % fname)
        
    try:
        ncols = data.shape[1]
    except:
        # if only one row in data file, convert to 2d array
        data = numpy.array([data])  
        ncols = data.shape[1]
        
    if 0:
        # no longer true since usecols is used
        if data.shape[1] != len(columns):
            raise Exception("len(columns) = %s but %s columns read from file" \
                     % (len(columns), ncols))

    nfaults = data.shape[0]
    print "Read %s faultslip datasets from %s" % (nfaults,fname)

    subfaults = []
    total_slip = 0.
    for k in range(nfaults):
        subfault_params = copy.copy(defaults)
        for j in range(ncols):
            jj = usecols[j]
            #print "+++ j, jj, param: ",j,jj,columns[jj]
            if columns[jj] != 'ignore':
                subfault_params[columns[jj]] = data[k,j]
            else:
                raise Exception("This shouldn't happen!")
        #import pdb; pdb.set_trace()
        for param in ['slip','length','width','depth']:
            # convert to meters if necessary:
            if units.get(param, 'm') == 'km':
                subfault_params[param] = subfault_params[param] * 1e3 
                #units[param] = 'm'
            elif units.get(param, 'm') == 'cm':
                subfault_params[param] = subfault_params[param] * 1e-2
                #units[param] = 'm'
                
        subfaults.append(subfault_params)
        subfault_slip = subfault_params['slip']*subfault_params['length'] \
                        * subfault_params['width']
        total_slip += subfault_slip
        #print "Subfault slip*length*width = ",subfault_slip
        #import pdb; pdb.set_trace()

    print "Total slip*length*width = ",total_slip
    if 0:
        for mu in [3e11, 4e11, 6e11]:
            Mo = 0.1*mu*total_slip  # 0.1 factor to convert to Nm
            Mw = 2./3. * (log10(Mo) - 9.1)
            print "With mu = %6.1e, moment magnitude is Mw = %5.2f" % (mu,Mw)
    return subfaults
    

def read_subfault_model_csv(path, delimiter=","):
    r"""Read a subfault model from a csv file

    Assumes that the first row gives the column headings, which should agree 
    with names in *valid_labels* defined in *read_subfault_model*.

    :Input:

    :Output:

    """

    sub_fault = SubFault()

    with open(path, 'r') as data_file:
        # Parse column and units
        columns = [label.strip() for label in data_file.readline().split(delimiter)]
        units = [unit.strip() for unit in data_file.readline().split(delimiter)]

    return read_subfault_model(path, columns, units, delimiter=delimiter, skiprows=2)

    return sub_fault



    fid = open(fname)
    c1 = fid.readline().split(',')
    columns = [c.strip() for c in c1]
    u1 = fid.readline().split(',')
    units_list = [u.strip() for u in u1]
    units = {}
    for i,c in enumerate(columns):
        units[c] = units_list[i]
    fid.close()
    subfaults = read_subfault_model(fname, columns, units, delimiter=',', \
                 skiprows=2)
    return subfaults


def read_subfault_model_ucsb(fname):
    """
    Read in subfault format models produced by Chen Ji's group at UCSB,
    downloadable from:
     http://www.geol.ucsb.edu/faculty/ji/big_earthquakes/home.html
    """
    lines = open(fname,'r').readlines()
    regexp_dx = re.compile(r"Dx=[ ]*(?P<dx>[^k]*)")
    regexp_dy = re.compile(r"Dy=[ ]*(?P<dy>[^k]*)")
    regexp_nx = re.compile(r"nx[^=]*=[ ]*(?P<nx>[^D]*)")
    regexp_ny = re.compile(r"ny[^=]*=[ ]*(?P<ny>[^D]*)")
    
    i_dxdy = -1
    for i,line in enumerate(lines):
        result_dx = regexp_dx.search(line)
        result_dy = regexp_dy.search(line)
        result_nx = regexp_nx.search(line)
        result_ny = regexp_ny.search(line)
        if result_dx and result_dy:
            i_dxdy = i
            dx = float(result_dx.group('dx'))
            dy = float(result_dy.group('dy'))
            nx = int(result_nx.group('nx'))
            ny = int(result_ny.group('ny'))
            print "Found dx = %s, dy = %s, nx = %s, ny = %s in line %s" \
                  % (dx, dy, nx, ny, i_dxdy)
            
    if i_dxdy == -1:
        print "*** Did not find a line containing both Dx and Dy"
        raise Exception()
    
    
    defaults = {'length': dx, 'width' : dy, 'latlong_location':'centroid'}
    units = {'slip': 'cm', 'depth': 'km', 'length': 'km', 'width': 'km'}

    columns = """latitude longitude depth slip rake strike dip rupture_time
                 rise_time rise_time_ending ignore""".split()

    subfaults = read_subfault_model(fname, columns=columns, units=units, \
                 defaults = defaults, skiprows= 9+i_dxdy)

    return subfaults

def set_fault_xy(faultparams):
    longitude = faultparams['longitude'] 
    latitude = faultparams['latitude'] 
    dip = faultparams['dip'] 
    strike = faultparams['strike'] 
    slip = faultparams['slip'] 
    rake = faultparams['rake'] 
    length = faultparams['length'] 
    width = faultparams['width'] 
    depth = faultparams['depth'] 
    location = faultparams['latlong_location'] 
    ang_strike = strike*DEG2RAD
    ang_dip = dip*DEG2RAD
    ang_rake = rake*DEG2RAD
    x0 = longitude
    y0 = latitude

    if location == "top center":

        depth_top = depth
        depth_centroid = depth + 0.5*width*sin(ang_dip)
        depth_bottom = depth + width*sin(ang_dip)

        # Convert fault origin from top of fault plane to bottom:
        del_x = width*cos(ang_dip)*cos(ang_strike) / (lat2meter*cos(y0*DEG2RAD))
        del_y = -width*cos(ang_dip)*sin(ang_strike) / lat2meter
        
        x_top = x0
        y_top = y0 
        x_bottom = x0+del_x
        y_bottom = y0+del_y
        x_centroid = x0+0.5*del_x
        y_centroid = y0+0.5*del_y

    elif location == "centroid":

        depth_top = depth - 0.5*width*sin(ang_dip)
        depth_centroid = depth 
        depth_bottom = depth + 0.5*width*sin(ang_dip)

        del_x = 0.5*width*cos(ang_dip)*cos(ang_strike) / (lat2meter*cos(y0*DEG2RAD))
        del_y = -0.5*width*cos(ang_dip)*sin(ang_strike) / lat2meter

        x_centroid = x0
        y_centroid = y0
        x_top = x0-del_x
        y_top = y0-del_y
        x_bottom = x0+del_x
        y_bottom = y0+del_y

    else:
        raise ValueError("Unrecognized latlong_location %s " % location)

    # distance along strike from center of an edge to corner:
    dx2 = 0.5*length*sin(ang_strike) / (lat2meter*cos(y_bottom*DEG2RAD))
    dy2 = 0.5*length*cos(ang_strike) / lat2meter
    x_corners = [x_bottom-dx2,x_top-dx2,x_top+dx2,x_bottom+dx2,x_bottom-dx2]
    y_corners = [y_bottom-dy2,y_top-dy2,y_top+dy2,y_bottom+dy2,y_bottom-dy2]

    paramlist = """x_top y_top x_bottom y_bottom x_centroid y_centroid
        depth_top depth_bottom depth_centroid x_corners y_corners""".split()

    for param in paramlist:
        cmd = "faultparams['%s'] = %s" % (param,eval(param))
        exec(cmd)


def make_okada_final_dz(faultparamss, dtopo_params):
    """
    Create dz at final time (sum of dz from all subfaults).
    """

    mx = dtopo_params['mx']
    my = dtopo_params['my']
    
    x=numpy.linspace(dtopo_params['xlower'],dtopo_params['xupper'],mx)
    y=numpy.linspace(dtopo_params['ylower'],dtopo_params['yupper'],my)
    dz = numpy.zeros((my,mx))
    
    print "Making Okada dz for each of %s subfaults" \
            % len(subfaults)

    for k,subfault in enumerate(subfaults):
            sys.stdout.write("%s.." % k)
            sys.stdout.flush()
            dz = dz + okadamap(subfault, x, y)
            
    sys.stdout.write("\nDone\n")
    return x,y,dz


def rise_fraction(t, t0, t_rise, t_rise_ending=None):
    """
    A continuously differentiable piecewise quadratic function of t that is 
       0 for t <= t0, 
       1 for t >= t0 + t_rise + t_rise_ending 
    with maximum slope at t0 + t_rise.
    """

    if t_rise_ending is None: 
        t_rise_ending = t_rise

    t1 = t0+t_rise
    t2 = t1+t_rise_ending

    rf = where(t<=t0, 0., 1.)
    if t2==t0:
        return rf

    t20 = float(t2-t0)
    t10 = float(t1-t0)
    t21 = float(t2-t1)

    c1 = t21 / (t20*t10*t21) 
    c2 = t10 / (t20*t10*t21) 

    rf = where((t>t0) & (t<=t1), c1*(t-t0)**2, rf)
    rf = where((t>t1) & (t<=t2), 1. - c2*(t-t2)**2, rf)

    return rf

def test_rise_fraction(t0=0,t_rise=5,t_rise_ending=20):
    t2 = t0 + 2.*t_rise 
    t = numpy.linspace(t0-10., t2+10., 1001)
    rf = rise_fraction(t,t0,t_rise, t_rise_ending)
    figure(0)
    clf()
    plot(t, rf, 'b')
    ylim(-.1, 1.1)

def write_dtopo_header(dtopo_params):
    """
    Write header for dtopotype 3 files.
    """

    dx = float(dtopo_params['xupper'] - dtopo_params['xlower']) \
            / (dtopo_params['mx'] - 1)
    dy = float(dtopo_params['yupper'] - dtopo_params['ylower']) \
            / (dtopo_params['my'] - 1)
    dt = float(dtopo_params['tfinal'] - dtopo_params['t0']) \
            / (dtopo_params['ntimes'] - 1)

    fid = open(dtopo_params['fname'],'w')
    fid.write("%7i       mx \n" % dtopo_params['mx'])
    fid.write("%7i       my \n" % dtopo_params['my'])
    fid.write("%7i       mt \n" % dtopo_params['ntimes'])
    fid.write("%20.14e   xlower\n" % dtopo_params['xlower'])
    fid.write("%20.14e   ylower\n" % dtopo_params['ylower'])
    fid.write("%20.14e   t0\n" % dtopo_params['t0'])
    fid.write("%20.14e   dx\n" % dx)
    fid.write("%20.14e   dy\n" % dy)
    fid.write("%20.14e   dt\n" % dt)
    return fid

def make_dtopo_from_subfaults(subfaults, dtopo_params):
    
    """
    Create a dtopo file named *fname* based on applying the Okada
    model to each of the subfaults specified in the list *subfaults*.
    This should be a list of dictionaries that specify the necessary 
    Okada parameters and optionally 'rupture_time', 'rise_time',
    'rise_time_ending'.

    The dictionary *dtopo_params* specifies the parameters for the
    dtopo file created.  It must have entries for
       'fname' : file name for dtopo file
    and
       'mx', 'my', 'xlower', 'xupper', 'ylower', 'yupper' 
    specifying the area over which the dz values will be output on an mx by my
    grid.  Note that dz values are at vertices, so to output on a 1 minute
    grid, for example, set
       mx = 60*(xupper - xlower) + 1 

    Also 
        dtopo_params['t0'] initial time for output
        dtopo_params['tfinal'] final time for output
        dtopo_params['ntimes'] >= 2, number of output times
    The seafloor deformation dz will be output at ntimes equally spaced
    times between t0 and tfinal.  

    If dtopo_params['faulttype'] == 'static'
        then dz at t0 will be 0 and dz at tfinal will be the total dz
        determined by all the subfaults, ignoring the rupture_time and
        rise_time parameters, if any.  Normally in this case it is best to
        specify ntimes = 2 and t0=0, tfinal=1 second, for example.
        If larger ntimes is specified then linear interpolation is used to
        determine dz at intermediate times.
    If dtopo_params['faulttype'] == 'kinematic' or 'dynamic'
        then it is assumed that the rupture_time and rise_time are set and
        will be used to interpolate to the specified output times.

    dtopo_params['dtopotype'] = 3 is assumed currently but this could be
        extended to allow output with other dtopotypes.
    """

    fname = dtopo_params['fname']
    faulttype = dtopo_params.get('faulttype', 'static')
    dtopotype = dtopo_params.get('dtopotype', 3)
    t0 = dtopo_params.get('t0', 0.)
    tfinal = dtopo_params.get('tfinal', 1.)
    ntimes = dtopo_params.get('ntimes', 2)
    mx = dtopo_params['mx']
    my = dtopo_params['my']
    xlower = dtopo_params['xlower']
    ylower = dtopo_params['ylower']
    xupper = dtopo_params['xupper']
    yupper = dtopo_params['yupper']
    x=numpy.linspace(xlower,xupper,mx)
    y=numpy.linspace(ylower,yupper,my)
    times = numpy.linspace(t0,tfinal,ntimes)


    plot_rupture = False

    dZ = numpy.zeros((my,mx))
    dz_list = []

    if dtopotype in [3]:
        fid = write_dtopo_header(dtopo_params)
    else:
        raise Exception("Unsupported dtopotype: %s" % dtopotype)

    if faulttype == 'static':

        print "Making Okada dz for each of %s subfaults" \
                % len(subfaults)
        for k,subfault in enumerate(subfaults):
            dZk = okadamap(subfault, x, y)
            sys.stdout.write("%s.." % k)
            sys.stdout.flush()
            dZ = dZ + dZk
        sys.stdout.write("\nDone\n")

        for t in times:
            alpha=(t-t0)/(tfinal-t0)
            dz_list.append(alpha*dZ)
            for j in range(my-1, -1, -1):
                format = mx*'%12.6e  '
                fid.write(format  % tuple(alpha*dZ[j,:]))
                fid.write('\n')


    elif faulttype in ['kinematic','dynamic']:

        t_prev = -1.e99
        for t in times:
            if plot_rupture:
                figure(5)
                clf()
            for k,subfault in enumerate(subfaults):
                t0 = subfault.get('rupture_time',0)
                t1 = subfault.get('rise_time',0.5)
                t2 = subfault.get('rise_time_ending',0)
                rf = rise_fraction(t0,t1,t2)
                dfrac = rf(t) - rf(t_prev)
                if dfrac > 0.:
                    if subfault.get('dZ',None) is None:
                        subfault.dZ = okadamap(subfault, x, y)
                        print '+++ Applying Okada to subfault %s at t = %s' \
                            % (k,t)
                    dZ = dZ + dfrac * subfault.dZ

                if plot_rupture:
                    xc = subfault.longitude
                    yc = subfault.latitude
                    rise = rf(t)
                    if rise==0: plot(xc,yc,'wo')
                    elif rise==1: plot(xc,yc,'ko')
                    else: plot(xc,yc,'ro')
            if plot_rupture:
                clines = list(arange(-10.5,0.,0.5)) + list(arange(0.5,11,0.5))
                contour(x,y,dZ,clines)
                title('time t = %8.3f' % t)
                draw()

            dz_list.append(dZ)

            for j in range(my-1, -1, -1):
                fid.write(mx*'%012.6e  '  % dZ[j,:])
                fid.write('\n')
    else:   
        raise Exception("Unrecognized faulttype: %s" % faulttype)

    print "Created ",fname
    fid.close()

    dtopo = DTopo()
    dtopo.dtopo_params = dtopo_params
    dtopo.x = x
    dtopo.y = y
    dtopo.times = times
    dtopo.dz_list = dz_list
    dtopo.subfaults = subfaults

    return dtopo


def write_dz(fname,X,Y,dZ,t0=0.,tend=1.,ntimes=2):
    """
    Create dtopo file with instantaneous slip distributed at ntimes between
    t0 and tend.
    """
    fid = open(fname, 'w')
    for t in numpy.linspace(t0,tend,ntimes):
        alpha=(t-t0)/(tend-t0)
        
        for jj in range(len(Y)):
            j=-1-jj
            for i in range(len(X)) :
                fid.write('%012.6e %012.6e %012.6e %012.6e \n' \
                       % (t,X[i],Y[j],alpha*dZ[j,i]))

    fid.close()
    print "Created ",fname

    
def write_dz_witht(fname,X,Y,dZ,times):
    fid = open(fname, 'w')
    for it in range(len(times)):
        for jj in range(len(Y)):
            j=-1-jj
            for i in range(len(X)) :
                fid.write('%012.6e %012.6e %012.6e %012.6e \n' \
                     % (times[it],X[i],Y[j],dZ[it,j,i]))

    fid.close()
    print "Created ",fname


def read_dtopo_old(fname, deftype=None, only_last=True):
    """
    Read in a deformation file and return arrays X,Y,dZ.
    If deftype == static then data is at a single time and columns of file are x,y,dz.
    If deftype == dynamic then multiple times are in the file and columns are t,x,y,dz.
       In this case dZ for the last time is returned.   **Need to improve!**
    If deftype is None, then attempt to determine from file extension.
    """
    if deftype is None:
        # try to determine from file extension:
        fext = os.path.splitext(fname)[1]
        if fext in ['.tt1', '.txyz', '.txydz']:
            deftype = 'dynamic'
        elif fext=='.xyz':
            deftype = 'static'
        else:
            print "*** Error determining deftype from file extension"
            return

    if deftype == 'dynamic':
        d = np.loadtxt(fname)
        print "Loaded file %s with %s lines" %(fname,d.shape[0])
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
        if only_last:
            print "Using only last time with mx=%s and my=%s" % (mx,my)
            if mx*my != len(lastlines):
                raise Exception("*** Error in determining mx and my!\nlen(lastlines)=%s" \
                                  % len(lastlines))      
            X = np.reshape(lastlines[:,1],(my,mx))
            Y = np.reshape(lastlines[:,2],(my,mx))
            dZ = np.reshape(lastlines[:,3],(my,mx))
        else:
            X = np.reshape(lastlines[:,1],(my,mx))
            Y = np.reshape(lastlines[:,2],(my,mx))
            dZ = []
            print "Returning dZ as a list of mx*my arrays"
            for n in range(ntimes):
                i1 = n*mx*my
                i2 = (n+1)*mx*my
                dZ.append(np.reshape(d[i1:i2,3],(my,mx)))
    elif deftype == 'static':  
        d = np.loadtxt(fname)  
        xvals = list(set(d[:,1]))
        xvals.sort()
        mx = len(xvals)
        my = len(d) / mx
        print "Read dtopo: mx=%s and my=%s" % (mx,my)
        if mx*my != len(d):
            raise Exception("*** Error in determining mx and my!\nlen(d)=%s" \
                              % len(d))
        X = np.reshape(d[:,0],(my,mx))
        Y = np.reshape(d[:,1],(my,mx))
        dZ = np.reshape(d[:,2],(my,mx))
    else:
        print "*** Unrecognized deftype: ",deftype
        raise

    return X,Y,dZ

def read_dtopo(fname, dtopotype):
    if dtopotype==3:
        fid = open(fname)
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

        dZvals = np.loadtxt(fname, skiprows=9)
        dz_list = []
        for k,t in enumerate(times):
            dZk = np.reshape(dZvals[k*my:(k+1)*my, :], (my,mx))
            dZk = np.flipud(dZk)
            dz_list.append(dZk)
            
        dtopo = DTopo()
        dtopo.mx = mx
        dtopo.my = my
        dtopo.x = x
        dtopo.y = y
        dtopo.times = times
        dtopo.dz_list = dz_list
    else:
        raise Exception("*** Unrecognized dtopotype: %s" % dtopotype)

    return dtopo
        

def plot_subfaults(subfaults, plot_centerline=False, slip_color=False, \
            cmap_slip=None, cmin_slip=None, cmax_slip=None, \
            plot_rake=False, xylim=None, plot_box=True):

    """
    Plot each subfault projected onto the surface.
    Describe parameters...
    """

    #figure(44,(6,12)) # For CSZe01
    #clf()

    # for testing purposes, make random slips:
    test_random = False


    max_slip = 0.
    min_slip = 0.
    for subfault in subfaults:
        if test_random:
            subfault['slip'] = 10.*rand()  # for testing
            #subfault['slip'] = 8.  # uniform
        slip = subfault['slip']
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
    for subfault in subfaults:

        set_fault_xy(subfault)

        # unpack parameters:
        paramlist = """x_top y_top x_bottom y_bottom x_centroid y_centroid
            depth_top depth_bottom depth_centroid x_corners y_corners""".split()

        for param in paramlist:
            cmd = "%s = subfault['%s']" % (param,param)
            exec(cmd)

        y_ave += y_centroid


        # Plot projection of planes to x-y surface:
        if plot_centerline:
            plot([x_top],[y_top],'bo',label="Top center")
            plot([x_centroid],[y_centroid],'ro',label="Centroid")
            plot([x_top,x_centroid],[y_top,y_centroid],'r-')
        if plot_rake:
            if test_random:
                subfault['rake'] = 90. + 30.*(rand()-0.5)  # for testing
            tau = (subfault['rake'] - 90) * pi/180.
            plot([x_centroid],[y_centroid],'go',markersize=5,label="Centroid")
            dxr = x_top - x_centroid
            dyr = y_top - y_centroid
            x_rake = x_centroid + cos(tau)*dxr - sin(tau)*dyr
            y_rake = y_centroid + sin(tau)*dxr + cos(tau)*dyr
            plot([x_rake,x_centroid],[y_rake,y_centroid],'g-',linewidth=1)
        if slip_color:
            slip = subfault['slip']
            #c = cmap_slip(0.5*(cmax_slip + slip)/cmax_slip)
            #c = cmap_slip(slip/cmax_slip)
            s = min(1, max(0, (slip-cmin_slip)/(cmax_slip-cmin_slip)))
            c = cmap_slip(s*.99)  # since 1 does not map properly with jet
            fill(x_corners,y_corners,color=c,edgecolor='none')
        if plot_box:
            plot(x_corners, y_corners, 'k-')

    slipax = gca()
        
    y_ave = y_ave / len(subfaults)
    slipax.set_aspect(1./cos(y_ave*pi/180.))
    ticklabel_format(format='plain',useOffset=False)
    xticks(rotation=80)
    if xylim is not None:
        axis(xylim)
    title('Fault planes')
    if slip_color:
        cax,kw = mpl.colorbar.make_axes(slipax)
        norm = mpl.colors.Normalize(vmin=cmin_slip,vmax=cmax_slip)
        cb1 = mpl.colorbar.ColorbarBase(cax, cmap=cmap_slip, norm=norm)
        #import pdb; pdb.set_trace()
    sca(slipax) # reset the current axis to the main figure


def plot_subfaults_depth(subfaults):
    """
    Plot the depth of each subfault vs. x in one plot and vs. y in a second plot.
    """

    for subfault in subfaults:

        set_fault_xy(subfault)

        # unpack parameters:
        paramlist = """x_top y_top x_bottom y_bottom x_centroid y_centroid
            depth_top depth_bottom depth_centroid x_corners y_corners""".split()

        for param in paramlist:
            cmd = "%s = subfault['%s']" % (param,param)
            exec(cmd)

        # Plot planes in x-z and y-z to see depths:
        subplot(211)
        plot([x_top,x_bottom],[-depth_top,-depth_bottom])
        subplot(212)
        plot([y_top,y_bottom],[-depth_top,-depth_bottom])

    subplot(211)
    title('depth vs. x')
    subplot(212)
    title('depth vs. y')

def plot_dz_contours(x,y,dz,dz_interval=0.5):
    dzmax = max(dz.max(), -dz.min()) + dz_interval
    clines1 = np.arange(dz_interval, dzmax, dz_interval)
    clines = list(-np.flipud(clines1)) + list(clines1)

    print "Plotting contour lines at: ",clines
    contour(x,y,dz,clines,colors='k')

def plot_dz_colors(x,y,dz,cmax_dz=None,dz_interval=None):
    """
    Plot sea floor deformation dtopo as colormap with contours
    """

    from clawpack.visclaw import colormaps

    dzmax = abs(dz).max()
    if cmax_dz is None:
        cmax_dz = dzmax
    cmap = colormaps.blue_white_red
    pcolor(x, y, dz, cmap=cmap)
    clim(-cmax_dz,cmax_dz)
    cb2 = colorbar(shrink=1.0)
    
    if dz_interval is None:
        dz_interval = cmax_dz/10.
    clines1 = np.arange(dz_interval, dzmax + dz_interval, dz_interval)
    clines = list(-np.flipud(clines1)) + list(clines1)
    print "Plotting contour lines at: ",clines
    contour(x,y,dz,clines,colors='k',linestyles='solid')
    y_ave = 0.5*(y.min() + y.max())
    gca().set_aspect(1./cos(y_ave*pi/180.))

    ticklabel_format(format='plain',useOffset=False)
    xticks(rotation=80)
    title('Seafloor deformation')


def strike_direction(x1,y1,x2,y2):
    """
    Calculate strike direction between two points.
    Actually calculates "initial bearing" from (x1,y1) in direction
    towards (x2,y2), following
        http://www.movable-type.co.uk/scripts/latlong.html
    """

    from numpy import pi,sin,cos,arctan2
    x1 = x1*pi/180.
    y1 = y1*pi/180.
    x2 = x2*pi/180.
    y2 = y2*pi/180.
    dx = x2-x1
    theta = arctan2(sin(dx)*cos(y2), \
             cos(y1)*sin(y2) - sin(y1)*cos(y2)*cos(dx))
    s = theta*180./pi
    if s<0:
        s = 360+s
    return s


# ==============================================================================
#  Okada Functionality
#
#  Okada model is a mapping from several fault parameters
#  to a surface deformation.
#  See Okada 1985, or Okada 1992, Bull. Seism. Soc. Am.
#  
#  Some routines adapted from fortran routines written by
#  Xiaoming Wang.
#  
#  okadamap function riginally written in Python by Dave George in okada.py.
#  Rewritten and made more flexible by Randy LeVeque:
#  Location can be specified as "top center" or "centroid".
#  
#  The main function is okadamap(okadaparams,X,Y).
# ==============================================================================
# Constants
poisson = 0.25

def okadamap(okadaparams,X,Y):

    """
    create displacement matrix dZ for a surface displacement
    over gridded region defined by X,Y, vectors of length nx,ny
    given okadaparams
    """

    # rad = pi/180.       # conversion factor from degrees to radians
    # rr = 6.378e6       # radius of earth -- original code
    #rr = Rearth         # should use this instead!   
    # lat2meter = rr*rad  # conversion factor from degrees latitude to meters

    hh =  okadaparams["depth"]
    L  =  okadaparams["length"]
    w  =  okadaparams["width"]
    d  =  okadaparams["slip"]
    th =  okadaparams["strike"]
    dl =  okadaparams["dip"]
    rd =  okadaparams["rake"]
    y0 =  okadaparams["latitude"]
    x0 =  okadaparams["longitude"]
    location =  okadaparams.get("latlong_location", "top center")

    ang_dip = DEG2RAD*dl
    ang_slip = DEG2RAD*rd
    ang_strike = DEG2RAD*th
    halfL = 0.5*L

    plot_plane = False
    print_xy = False

    if plot_plane:
        figure(202)
        #clf()

    if print_xy:
        print "x0,y0: ",x0,y0


    if location == "top center":

        # Convert focal depth used for Okada's model
        # from top of fault plane to bottom:

        depth_top = hh
        hh = hh + w*sin(ang_dip)
        depth_bottom = hh

        # Convert fault origin from top of fault plane to bottom:
        del_x = w*cos(ang_dip)*cos(ang_strike) / (lat2meter*cos(y0*DEG2RAD))
        del_y = -w*cos(ang_dip)*sin(ang_strike) / lat2meter
        
        x_top = x0
        y_top = y0 
        x_bottom = x0+del_x
        y_bottom = y0+del_y
        x_centroid = x0+0.5*del_x
        y_centroid = y0+0.5*del_y


    elif location == "centroid":

        # Convert focal depth used for Okada's model
        # from middle of fault plane to bottom:
        depth_top = hh - 0.5*w*sin(ang_dip)
        hh = hh + 0.5*w*sin(ang_dip)
        depth_bottom = hh

        # Convert fault origin from middle of fault plane to bottom:
        del_x = 0.5*w*cos(ang_dip)*cos(ang_strike) / (lat2meter*cos(y0*DEG2RAD))
        del_y = -0.5*w*cos(ang_dip)*sin(ang_strike) / lat2meter

        x_centroid = x0
        y_centroid = y0
        x_top = x0-del_x
        y_top = y0-del_y
        x_bottom = x0+del_x
        y_bottom = y0+del_y

    else:
        raise ValueError("Unrecognized latlong_location" % location)

    # adjust x0,y0 to bottom center of fault plane:
    x0 = x0 + del_x
    y0 = y0 + del_y

    # distance along strike from center of an edge to corner:
    dx2 = 0.5*L*sin(ang_strike) / (lat2meter*cos(y_bottom*DEG2RAD))
    dy2 = 0.5*L*cos(ang_strike) / lat2meter

    if print_xy:
        print "del_x, del_y: ",del_x,del_y
        print "original x0,y0: ",x0,y0
        print "bottom: ",x_bottom, y_bottom
        print "centroid: ",x_centroid, y_centroid
        print "top: ",x_top, y_top
        print "dx2,dy2: ",dx2,dy2
    if plot_plane:
        figure(203)
        subplot(211)
        plot([x_top,x_bottom],[-depth_top,-depth_bottom])
        title('depth vs. x')
        subplot(212)
        plot([y_top,y_bottom],[-depth_top,-depth_bottom])
        title('depth vs. y')
        #ylim([-100,0])
        figure(202)
        plot([x_top],[y_top],'bo',label="Top center")
        plot([x_centroid],[y_centroid],'ro',label="Centroid")
        plot([x_top,x_centroid],[y_top,y_centroid],'r-')
        plot([x_bottom-dx2,x_top-dx2,x_top+dx2,x_bottom+dx2,x_bottom-dx2],\
             [y_bottom-dy2,y_top-dy2,y_top+dy2,y_bottom+dy2,y_bottom-dy2],'b-')
        axis('scaled')
        axis([X[0],X[-1],Y[0],Y[-1]])
        title("Blue: top center, Red: centroid of subfault")


    x,y = numpy.meshgrid(X,Y)

    # Convert distance from (x,y) to (x_bottom,y_bottom) from degrees to meters:
    xx = lat2meter*cos(DEG2RAD*y)*(x-x_bottom)   
    yy = lat2meter*(y-y_bottom)


    # Convert to distance along strike (x1) and dip (x2):
    x1 = xx*sin(ang_strike) + yy*cos(ang_strike) 
    x2 = xx*cos(ang_strike) - yy*sin(ang_strike) 

    # In Okada's paper, x2 is distance up the fault plane, not down dip:
    x2 = -x2

    if 0:
        figure(203)
        clf()
        plot([xx[0,0],xx[0,-1],xx[-1,-1],xx[-1,0],xx[0,0]], \
             [yy[0,0],yy[0,-1],yy[-1,-1],yy[-1,0],yy[0,0]], 'k-')
        
        plot([x1[0,0],x1[0,-1],x1[-1,-1],x1[-1,0],x1[0,0]], \
             [x2[0,0],x2[0,-1],x2[-1,-1],x2[-1,0],x2[0,0]], 'b-')
        
    p = x2*cos(ang_dip) + hh*sin(ang_dip)
    q = x2*sin(ang_dip) - hh*cos(ang_dip)

    f1=strike_slip (x1+halfL,p,  ang_dip,q)
    f2=strike_slip (x1+halfL,p-w,ang_dip,q)
    f3=strike_slip (x1-halfL,p,  ang_dip,q)
    f4=strike_slip (x1-halfL,p-w,ang_dip,q)

    g1=dip_slip (x1+halfL,p,  ang_dip,q)
    g2=dip_slip (x1+halfL,p-w,ang_dip,q)
    g3=dip_slip (x1-halfL,p,  ang_dip,q)
    g4=dip_slip (x1-halfL,p-w,ang_dip,q)

    # Displacement in direction of strike and dip:
    ds = d*cos(ang_slip)
    dd = d*sin(ang_slip)

    us = (f1-f2-f3+f4)*ds
    ud = (g1-g2-g3+g4)*dd

    dZ = (us+ud)

    if 0:
        contour(x,y,dZ,numpy.linspace(-8,8,17),colors='k')

    return dZ


#===========================================================================
def strike_slip (y1,y2,ang_dip,q):
    """
    !.....Used for Okada's model
    !.. ..Methods from Yoshimitsu Okada (1985)
    !-----------------------------------------------------------------------
    """
    sn = sin(ang_dip)
    cs = cos(ang_dip)
    d_bar = y2*sn - q*cs
    r = sqrt(y1**2 + y2**2 + q**2)
    xx = sqrt(y1**2 + q**2)
    a4 = 2.0*poisson/cs*(numpy.log(r+d_bar) - sn*numpy.log(r+y2))
    f = -(d_bar*q/r/(r+y2) + q*sn/(r+y2) + a4*sn)/(2.0*3.14159)

    return f


#============================================================================
def dip_slip (y1,y2,ang_dip,q):
    """
    !.....Based on Okada's paper (1985)
    !.....Added by Xiaoming Wang
    !-----------------------------------------------------------------------
    """
    sn = sin(ang_dip)
    cs = cos(ang_dip)

    d_bar = y2*sn - q*cs;
    r = sqrt(y1**2 + y2**2 + q**2)
    xx = sqrt(y1**2 + q**2)
    a5 = 4.*poisson/cs*numpy.arctan((y2*(xx+q*cs)+xx*(r+xx)*sn)/y1/(r+xx)/cs)
    f = -(d_bar*q/r/(r+y1) + sn*numpy.arctan(y1*y2/q/r) - a5*sn*cs)/(2.0*3.14159)

    return f


#============================================================================
def filtermask (dZ,faultparams):
    """
    borrowed from code written by Xiaoming Wang and Tom Logan at ARSC

    !.....Filter the deformation using a circular mask centered
    !.....at the epicenter using a calculated radius
    !.....Removes small numerical artifacts away from the epicenter
    """
    filterindices=[]

    osixty = 0.016666666667
    #rad = 0.01745329252
    #rr = 6.378e6       # original code
    #rr = Rearth         # should use this instead!   

    xo = faultparams['xlower']
    yo = faultparams['ylower']
    nx = faultparams['mx']
    ny = faultparams['my']
    spacing = faultparams['dx']

    x0 = faultparams['longitude']
    y0 = faultparams['latitude']
    l =  faultparams['length']
    w =  faultparams['width']
    dl = faultparams['dip']


    ang_dip = DEG2RAD*dl # convert degree to radian

    #!-- fault origin in pixels -----------
    ypix = (y0-yo)/spacing
    xpix = (x0-xo)/spacing

    #!-- conversion from meters to pixels ---
    tmpd=spacing*DEG2RAD
    xdist = tmpd*Rearth

    #!-- size of the fault in pixels --------
    npix_x = l/xdist
    npix_y = w/xdist

    #!-- set the range (radius) of the filter circle --------
    #!----- for small dip angles, use the length and width --
    #!----- for larger dip angles, use only the length ------

    if dl<30.0:
        drange = 1.5 * cos(ang_dip)*sqrt(npix_x*npix_x+npix_y*npix_y)
    else:
        drange = 1.2 * npix_x

    print("Filtering deformation using a circle of radius %s" % (drange))

    #!-- Create the filtering mask ----------
    for i in xrange(nx):
        for j in xrange(ny) :
            dist = sqrt((i+1-xpix)**2+(j+1-ypix)**2)
            if dist > drange :
                filterindices.append((j,i))

    #!-- apply the filter to the actual deformation ------
    dZ = filterdata(dZ,filterindices,radius=2)

    return dZ


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
#  Sub-Fault Class 
# ==============================================================================
class SubFault(object):

    r"""Class representing a single subfault.

    :TODO:
     - Support something other than lat-long
     - Provide plots (and other plot types)
     - Provide detailed documentation


    Subfault Parameters
    -------------------
     - *dimensions* (list) - Dimensions of the fault plane.
     - *coordinates* (list) - Longitude-latitude of some point on the fault 
       plane.  Point is specified by *coordinate_specification*.
     - *coordinate_specification* (string) - Specifies location relative to the
       fault of the *coordinates* values.  Valid options include "top center",
       "bottom center", and "centroid".
     - *depth* (float) - Depth of the specified point below the sea-floor.
     - *strike* (float) - Orientation of the top edge, measured in degrees 
       clockwise from North.  Between 0 and 360. The fault plane dips downward 
       to the right when moving along the top edge in the strike direction.
     - *dip* (float) - Angle at which the plane dips downward from the top edge,
       a positive angle between 0 and 90 degrees.
     - *rake* (float) - Angle in the fault plane in which the slip occurs,
       measured in degrees counterclockwise from the strike direction. Between 
       -180 and 180.
     - *slip* (float) Positive distance the hanging block moves relative to the 
       foot block in the direction specified by the rake. The "hanging block" is
       the one above the dipping fault plane (or to the right if you move in the 
       strike direction).

    Attributes
    ----------
     - *units* (dict) - Dictionary containing unit specifications for the 
       *coordinates*, *dimenstions*, *slip*, and *depth* subfault parameters.  
       Defaults to "lat-long", "m", "m", "m" respectively.

    Properties
    ----------
     :Note: All properties are in meters and do not match the units dictionary.

    """

    @property
    def x(self):
        r"""Coordinate array (x) for subfault."""
        if self._x is None:
            self.create_coordinate_arrays()
        return self._x
    @x.setter
    def x(self, value):
        self._x = value
    @x.deleter
    def x(self):
        del self._x
    @property
    def X(self):
        r"""2d x-coordinate array."""
        if self._X is None:
            self.create_2d_coordinate_arrays()
        return self._X
    @X.deleter
    def X(self):
        del self._X

    @property
    def y(self):
        r"""Coordinate array (y) for subfault."""
        if self._y is None:
            self.create_coordinate_arrays()
        return self._y
    @y.setter
    def y(self, value):
        self._y = value
    @y.deleter
    def y(self):
        del self._y
    @property
    def Y(self):
        r"""2d y-coordinate array."""
        if self._Y is None:
            self.create_2d_coordinate_arrays()
        return self._Y
    @Y.deleter
    def Y(self):
        del self._Y

    @property
    def dZ(self):
        r"""Deformation dZ of subfault."""
        if self._dZ is None:
            self.create_deformation_array()
        return self._dZ
    @dZ.setter
    def dZ(self, value):
        self._dZ = value
    @dZ.deleter
    def dZ(self):
        del self._dZ

    # Calculated geometry
    @property
    def fault_plane_corners(self):
        r"""Coordinates of the corners of the fault plane"""
        if self._fault_plane_corners is None:
            self.calculate_geometry()
        return self._fault_plane_corners
    @property
    def fault_plane_centers(self):
        r"""Coordinates along the center-line of the fault plane"""
        if self._fault_plane_centers is None:
            self.calculate_geometry()
        return self._fault_plane_centers

    def __init__(self, path=None, units={}):
        r"""SubFault initialization routine.
        
        See :class:`SubFault` for more info.

        """

        super(SubFault, self).__init__()

        # Object defered storage
        self._x = None
        self._X = None
        self._y = None
        self._Y = None
        self._dZ = None
        self._fault_plane_corners = None
        self._fault_plane_centers = None

        # Parameters for subfault specification
        self.coordinates = [] # longitude, latitude
        self.coordinate_specification = 'centroid' # 'centroid', 'top center', 'epicenter'
        self.dimensions = [0.0, 0.0] # [length, width]
        self.rake = None
        self.strike = None
        self.dip = None
        self.depth = None
        self.slip = None
        self.rupture_type = 'static' # 'static', 'dynamic', 'kinetic'
        self.t = numpy.array([0.0, 5.0, 10.0])

        # Default units of each parameter type
        self.units = {'coordinates':'lat-long', 'dimensions':'m', 'slip':'m',
                      'depth':"km"}
        self.units.update(units)

        # Read in file at path if provided
        if path is not None:
            self.read(path)


    def __str__(self):
        output = "Subfault Characteristics:\n"
        return output


    def convert2meters(self): 
        r"""Convert relevant lengths to correct units.

        Returns converted (dimensions, depth, slip) 

        """
        dimensions = [0.0, 0.0]
        depth = 0.0
        slip = 0.0

        conversion_dict = {"km":1e3, "cm":1e-2, "nm":1852.0, "m":1.0}

        dimensions = None
        depth = None
        slip = None
        if self.dimensions is not None:
            dimensions = [self.dimensions[0] * conversion_dict[self.units["dimensions"]],
                          self.dimensions[1] * conversion_dict[self.units["dimensions"]]]
        if self.depth is not None:
            depth = self.depth * conversion_dict[self.units["depth"]]
        if self.slip is not None:
            slip = self.slip * conversion_dict[self.units["slip"]]

        return dimensions, depth, slip


    def Mw(self, mu=5e11):
        r"""Calculate the effective moment magnitude of subfault."""

        dimensions, depth, slip = self.convert2meters()
        total_slip = dimensions[0] * dimensions[1] * slip
        Mo = 0.1 * mu * total_slip
        return 2.0 / 3.0 * (numpy.log10(Mo) - 9.0)


    def calculate_slip(self, Mw, mu=5e11):
        r"""Set slip based on a moment magnitude *Mw*."""

        Mo = 10.**(Mw * 3.0 / 2.0 + 9.1)

        dimensions, depth, slip = self.convert2meters()
        subfault_area = dimensions[0] * dimensions[1]
        
        self.slip = Mo / (0.1 * mu * subfault_area)

        # Convert back to requested units
        if self.units['slip'] == 'cm':
            self.slip *= 1e2
        elif self.units["slip"] == 'km':
            self.slip *= 1e-3


    def containing_rect(self):
        r"""Find containing rectangle of subfault in x-y plane.

        Returns tuple of x-limits and y-limits.

        """

        extent = [numpy.infty, -numpy.infty, numpy.infty, -numpy.infty]
        for corner in self.fault_plane_corners:
            extent[0] = min(corner[0], extent[0])
            extent[1] = max(corner[0], extent[1])
            extent[2] = min(corner[1], extent[2])
            extent[3] = max(corner[1], extent[3])

        return extent


    def transform(self, origin, point, theta):
        r"""Transform to rotated coordinate system (via strike).

        Input
        -----
         - *origin* (tuple) - Point being rotated about, (x,y).
         - *point* (tuple) - Point to be rotated, (x,y).
         - *theta* (float) - Angle of rotation in radians.

        Output
        ------
         - (tuple) The result of the transforming *point* about the point
           *origin*.

        """
        return (  (point[0] - origin[0]) * numpy.cos(-theta) 
                - (point[1] - origin[1]) * numpy.sin(-theta) 
                + origin[0],
                  (point[0] - origin[0]) * numpy.sin(-theta) 
                + (point[1] - origin[1]) * numpy.cos(-theta) 
                + origin[1]) 


    def calculate_geometry(self):
        r"""Calculate the fault geometry.

        Routine calculates the class attributes *fault_plane_corners* and 
        *fault_plane_centers* which are the corners of the fault plane and 
        points along the centerline respecitvely in 3D space.

        """

        # Simple conversion factor of latitude to meters
        lat2meter = topotools.dist_latlong2meters(0.0, 1.0)[1]

        # Setup coordinate arrays
        self._fault_plane_corners = [[None, None, None], # a 
                                     [None, None, None], # b
                                     [None, None, None], # c
                                     [None, None, None]] # d
        self._fault_plane_centers = [[None, None, None], # 1
                                     [None, None, None], # 2 
                                     [None, None, None]] # 3

        # Convert values to meters
        dimensions, depth, slip = self.convert2meters()

        # Depths (in meters)
        if self.coordinate_specification == 'top center':
            self._fault_plane_corners[0][2] = depth
            self._fault_plane_corners[1][2] = depth     \
                                 + dimensions[1] * sin(self.dip * DEG2RAD)
            self._fault_plane_centers[1][2] = depth      \
                           + 0.5 * dimensions[1] * sin(self.dip * DEG2RAD)
        elif self.coordinate_specification == 'centroid':
            self._fault_plane_corners[0][2] = depth     \
                                 - dimensions[1] * sin(self.dip * DEG2RAD)
            self._fault_plane_corners[1][2] = depth     \
                                 + dimensions[1] * sin(self.dip * DEG2RAD)
            self._fault_plane_centers[1][2] = depth
        elif self.coordinate_specification == 'bottom center':
            self._fault_plane_corners[0][2] = depth
            self._fault_plane_corners[1][2] = depth     \
                                 - dimensions[1] * sin(self.dip * DEG2RAD)
            self._fault_plane_centers[1][2] = depth      \
                           - 0.5 * dimensions[1] * sin(self.dip * DEG2RAD)
        self._fault_plane_corners[2][2] = self._fault_plane_corners[0][2]
        self._fault_plane_corners[3][2] = self._fault_plane_corners[1][2]
        self._fault_plane_centers[0][2] = self._fault_plane_corners[0][2]
        self._fault_plane_centers[2][2] = self._fault_plane_corners[1][2]

        # Convert dimensions to lat-long
        dimensions[0] *= 1.0 / lat2meter
        dimensions[1] *= 1.0 / lat2meter

        # Calculate xy-plane projected width
        xy_width = dimensions[1] * cos(self.dip * DEG2RAD)
        
        # Locate fault plane in 3D space
        # Note that the coodinate specification is in reference to the fault 
        #
        #    Top edge    Bottom edge
        #      a ----------- b
        #      |             |
        #      |             |
        #      |             |
        #      |             |
        #      1      2      3   L
        #      |             |
        #      |             |
        #      |             |
        #      |             |
        #      d ----------- c
        #            W
        #
        xy_corners = [None, None, None, None]
        xy_centers = [None, None, None]
        if self.coordinate_specification == 'top center':
            # Non-rotated locations of corners
            xy_corners[0] = (self.coordinates[0],
                             self.coordinates[1] + 0.5 * dimensions[0])
            xy_corners[1] = (self.coordinates[0] + xy_width,
                             self.coordinates[1] + 0.5 * dimensions[0])
            xy_corners[2] = (self.coordinates[0] + xy_width,
                             self.coordinates[1] - 0.5 * dimensions[0])
            xy_corners[3] = (self.coordinates[0],
                             self.coordinates[1] - 0.5 * dimensions[0])

            # Non-rotated lcoations of center-line coordinates
            xy_centers[0] = self.coordinates
            xy_centers[1] = (self.coordinates[0] + 0.5 * xy_width,
                             self.coordinates[1])
            xy_centers[2] = (self.coordinates[0] + xy_width,
                             self.coordinates[1])

        elif self.coordinate_specification == 'centroid':
            # Non-rotated locations of corners
            xy_corners[0] = (self.coordinates[0] - 0.5 * xy_width,
                             self.coordinates[1] + 0.5 * dimensions[0])
            xy_corners[1] = (self.coordinates[0] + 0.5 * xy_width,
                             self.coordinates[1] + 0.5 * dimensions[0])
            xy_corners[2] = (self.coordinates[0] + 0.5 * xy_width,
                             self.coordinates[1] - 0.5 * dimensions[0])
            xy_corners[3] = (self.coordinates[0] - 0.5 * xy_width,
                             self.coordinates[1] - 0.5 * dimensions[0])

            # Non-rotated lcoations of center-line coordinates
            xy_centers[0] = (self.coordinates[0] - 0.5 * xy_width,
                             self.coordinates[1])
            xy_centers[1] = self.coordinates
            xy_centers[2] = (self.coordinates[0] + 0.5 * xy_width,
                             self.coordinates[1])

        elif self.coordinate_specification == 'bottom center':
            # Non-rotated locations of corners
            xy_corners[0] = (self.coordinates[0] - xy_width,
                             self.coordinates[1] + 0.5 * dimensions[0])
            xy_corners[1] = (self.coordinates[0],
                             self.coordinates[1] + 0.5 * dimensions[0])
            xy_corners[2] = (self.coordinates[0],
                             self.coordinates[1] - 0.5 * dimensions[0])
            xy_corners[3] = (self.coordinates[0] - xy_width,
                             self.coordinates[1] - 0.5 * dimensions[0])

            # Non-rotated lcoations of center-line coordinates
            xy_centers[0] = (self.coordinates[0] - xy_width,
                             self.coordinates[1])
            xy_centers[1] = (self.coordinates[0] - 0.5 * xy_width,
                             self.coordinates[1])
            xy_centers[2] = self.coordinates

        else:
            raise ValueError("Unknown coordinate specification '%s'." % self.coordinate_specification)

        # Rotate fault plane corners and store
        # top_center = [xy_corners[0][0], xy_corners[0][1] - 0.5 * dimensions[0]]
        for (n, corner) in enumerate(xy_corners):
            self._fault_plane_corners[n][0:2] = \
                       self.transform(self.coordinates, corner, self.strike * DEG2RAD)
        for (n, center) in enumerate(xy_centers):
            self._fault_plane_centers[n][0:2] = \
                       self.transform(self.coordinates, center, self.strike * DEG2RAD)


    def create_coordinate_arrays(self, resolution=60, buffer_size=0.5):
        r"""Create coordinate arrays containing subfault.

        Input
        -----
         - *resolution* (int) - Number of grid points per degree.  Defaults to
           1" resolution.
         - *buffer_size* (float) - Buffer distance around edge of fault in 
           degrees, defaults to 0.5 degrees.

        """

        rect = self.containing_rect()
        rect[0] -= buffer_size
        rect[1] += buffer_size
        rect[2] -= buffer_size
        rect[3] += buffer_size

        self.delta = float(1.0 / resolution)
        N = [int((rect[1] - rect[0]) * resolution),
             int((rect[3] - rect[2]) * resolution)]

        self._x = numpy.linspace(rect[0],rect[1],N[0])
        self._y = numpy.linspace(rect[2],rect[3],N[1])


    def create_2d_coordinate_arrays(self):
        r"""Create 2d-coodrinate arrays."""
        self._X, self._Y = numpy.meshgrid(self.x, self.y)


    def create_deformation_array(self):
        r"""Create deformation array dZ.

        Use Okada model to calculate deformation from subfault parameters 
        contained in this object.

        Currently only calculates the vertical displacement.

        """
        import pdb; pdb.set_trace()
        dimensions, depth, slip = self.convert2meters()

        # Construct dictionary that okadamap is looking for
        okada_params = {}
        okada_params["depth"] = depth
        okada_params["length"] = dimensions[0]
        okada_params["width"] = dimensions[1]
        okada_params["slip"] = slip
        okada_params["strike"] = self.strike
        okada_params["dip"] = self.dip
        okada_params["rake"] = self.rake
        okada_params["longitude"] = self.coordinates[0]
        okada_params["latitude"] = self.coordinates[1]
        okada_params["latlong_location"] = self.coordinate_specification
        self._dZ = okadamap(okada_params, self.x, self.y)


    def write(self, path, topo_type=None):
        r"""Write out subfault characterization file to *path*.

        input
        -----
         - *path* (path) - Path to the output file to written to.
         - *topo_type* (int) - Type of topography file to write out.  Default is 1.

        """

        if topo_type is None:
            # Try to look at suffix for type
            extension = os.path.splitext(path)[1][1:]
            if extension[:2] == "tt":
                topo_type = int(extension[2])
            elif extension == 'xyz':
                topo_type = 1
            else:
                # Default to 3
                topo_type = 3

        if self.rupture_type != "static":
            raise NotImplemented("Only the 'static' rupture type is supported.")

        # Construct each interpolating function and evaluate at new grid
        with open(path, 'w') as data_file:

            if topo_type == 1:
                # Topography file with 4 columns, t, x, y, dz written from the upper
                # left corner of the region
                Y_flipped = numpy.flipud(self.Y)
                for (n, time) in enumerate(self.t):
                    alpha = (time - self.t[0]) / self.t[-1]
                    dZ_flipped = numpy.flipud(alpha * self.dZ[:,:])
                    for j in xrange(self.Y.shape[0]):
                        for i in xrange(self.X.shape[1]):
                            data_file.write("%s %s %s %s\n" % (self.t[n], self.X[j,i], Y_flipped[j,i], dZ_flipped[j,i]))
        
            elif topo_type == 2 or topo_type == 3:
                # Write out header
                data_file.write("%7i       mx \n" % self.x.shape[0])
                data_file.write("%7i       my \n" % self.y.shape[0])
                data_file.write("%7i       mt \n" % self.t.shape[0])
                data_file.write("%20.14e   xlower\n" % self.x[0])
                data_file.write("%20.14e   ylower\n" % self.y[0])
                data_file.write("%20.14e   t0\n" % self.t[0])
                data_file.write("%20.14e   dx\n" % self.delta)
                data_file.write("%20.14e   dy\n" % self.delta)
                data_file.write("%20.14e   dt\n" % float(self.t[1] - self.t[0]))

                if topo_type == 2:
                    raise ValueError("Topography type 2 is not yet supported.")
                elif topo_type == 3:
                    for (n, time) in enumerate(self.t):
                        alpha = (time - self.t[0]) / (self.t[-1])
                        for j in range(self.Y.shape[0]-1, -1, -1):
                            data_file.write(self.X.shape[1] * '%012.6e  ' 
                                                  % tuple(alpha * self.dZ[j,:]))
                            data_file.write("\n")

            else:
                raise ValueError("Only topography types 1, 2, and 3 are supported.")


    def read(self, path, file_type="usgs"):
        r"""Read in subfault specification from file at *path*

        :File Types:
         - USGS Sub-Fault Specification *USGS*.
         - Comma seperated values *CSV*.
         - Topography file *TOPOx* where *x* is the type of topography file.

        :Input:
         - *path* (path) - Path to subfault specification.
         - *file_type* (string) - String describing type of subfault 
           specification.  Default is "USGS".

        """

        raise NotImplemented("Subfault file reading not implemented yet.")


    def plot(self, axes=None, region_extent=None, contours=None, 
                   coastlines=None, limits=None, cmap=None):
        r"""Plot subfault deformation.


        Input
        -----
         - *axes* (`matplotlib.axes.Axes`) -
         - *region_extent* (list) - 
         - *contours* (list) -
         - *coastlines* (path) -
         - *limits* (list) -
         - *cmap* (`matplotlib.colors.Colormap`) -

        Output
        ------
         - *axes* (`matplotlib.axes.Axes`) - Axes used for plot, either
           this is created in this method if the input argument *axes* is None
           or the same object is passed back.

        """

        import matplotlib.pyplot as plt
        import matplotlib.colors as colors
        import clawpack.visclaw.colormaps as colormaps
        
        # Create axes object if needed
        if axes is None:
            fig = plt.figure()
            axes = fig.add_subplot(1,1,1)

        # Calculate plot extent and limits if not provided
        if region_extent is None:
            region_extent = ( numpy.min(self.x), numpy.max(self.x),
                              numpy.min(self.y), numpy.max(self.y) )
        if limits is None:
            depth_extent = [numpy.min(self.dZ),numpy.max(self.dZ)]
        else:
            depth_extent = limits

        # Setup axes labels, ticks and aspect
        axes.ticklabel_format(format="plain", useOffset=False)
        mean_lat = 0.5 * (region_extent[3] - region_extent[2])
        axes.set_aspect(1.0 / cos(numpy.pi / 180.0 * mean_lat))
        axes.set_title("Subfault Deformation")
        axes.set_xlabel("Longitude")
        axes.set_ylabel("Latitude")

        # Colormap and color norm
        if cmap is None:
            if depth_extent[0] >= 0.0:
                cmap = colormaps.make_colormap({0.0:'w', 1.0:'r'})
                extend = 'top'
            elif depth_extent[1] <= 0.0:
                cmap = colormaps.make_colormap({0.0:'b', 1.0:'w'})
                extend = 'bottom'
            else:
                cmap = colormaps.make_colormap({0.0:'b', 0.5:'w', 1.0:'r'})
                extend = 'both'
        # Equalize color extents
        if depth_extent[0] >= 0.0:
            depth_extent[0] = 0.0
        elif depth_extent[1] <= 0.0:
            depth_extent[1] = 0.0
        else:
            depth_extent[1] = max(-depth_extent[0], depth_extent[1])
            depth_extent[0] = -depth_extent[1]
        color_norm = colors.Normalize(depth_extent[0], depth_extent[1], clip=True)

        # Plot data
        if contours is not None:
            plot = axes.contourf(self.x, self.y, self.dZ, contours, cmap=cmap,
                                 extend=extend)
        else:
            plot = axes.imshow(numpy.flipud(self.dZ), vmin=depth_extent[0], 
                                                      vmax=depth_extent[1],
                                                      extent=region_extent,
                                                      cmap=cmap,
                                                      norm=color_norm)

        cbar = plt.colorbar(plot, ax=axes)
        cbar.set_label("Deformation (m)")

        # Plot coastlines
        if coastlines is not None:
            coastline_data = Topography(coastlines)
            axes.contour(coastline_data.X, coastline_data.Y, 
                         coastline_data.Z, levels=[0.0],colors='r')

        axes.set_xlim(region_extent[0:2])
        axes.set_ylim(region_extent[2:])

        return axes


    def plot_fault_rect(self, axes=None, color='r', markerstyle="o", 
                                         linestyle='-'):
        r"""Plot fault rectangle.

        Input
        -----
         - *axes* (`matplotlib.axes.Axes`) - 

        Output
        ------
         - (`matplotlib.axes.Axes`) - 

        """
        
        # Create axes object if needed
        if axes is None:
            fig = plt.figure()
            axes = fig.add_subplot(1,1,1)

        # Plot corners
        style = color + markerstyle
        for (n,corner) in enumerate(self.fault_plane_corners):
            axes.plot(corner[0], corner[1], style)
            axes.text(corner[0], corner[1], str(n+1))

        # Plot edges
        style = color + linestyle
        edges = []
        for edge in xrange(len(self.fault_plane_corners) - 1):
            edges.append([self.fault_plane_corners[edge][:2], 
                          self.fault_plane_corners[edge+1][:2]])
        edges.append([self.fault_plane_corners[-1][:2], 
                      self.fault_plane_corners[0][:2]])
        for edge in edges:
            axes.plot([edge[0][0], edge[1][0]], [edge[0][1], edge[1][1]], style)

        return axes


    def plot_rake(self, axes=None, color='r', markerstyle="o", linestyle='-'):
        r"""Plot fault rectangle.

        Input
        -----
         - *axes* (`matplotlib.axes.Axes`) - 

        Output
        ------
         - (`matplotlib.axes.Axes`) - 

        """
        
        # Create axes object if needed
        if axes is None:
            fig = plt.figure()
            axes = fig.add_subplot(1,1,1)

        centroid = self.fault_plane_centers[1][:2]
        top_edge = self.fault_plane_centers[0][:2]
        r = numpy.sqrt((top_edge[0] - self.fault_plane_corners[0][0])**2 +
                       (top_edge[1] - self.fault_plane_corners[0][1])**2 )
        theta = (self.strike + self.rake) * DEG2RAD
        xy_rake = (r * cos(-theta + 1.5 * numpy.pi) + centroid[0], 
                   r * sin(-theta + 1.5 * numpy.pi) + centroid[1])

        axes.annotate("",
            xy=xy_rake, xycoords='data',
            xytext=self.fault_plane_centers[1][:2], textcoords='data',
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3") )

        return axes


    def plot_subfault_depth(self, axes=None):
        r"""Plot the depth of each subfault vs. x in one plot and vs. y in a second plot.
        
        """

        raise NotImplemented("Subfault depth plots not implemented.")


    def plot_contours(self, axes=None, dz_interval=0.5):
        r"""Plot sea-floor deformation contours

        """

        dzmax = numpy.max(self.dZ.max(), -self.dZ.min()) + dz_interval
        clines1 = numpy.arange(dz_interval, dzmax, dz_interval)
        clines = list(-np.flipud(clines1)) + list(clines1)

        self.plot(axes=axes, contours=clines)

    def plot_seafloor(self, axes=None, cmax_dz=None, dz_interval=None):
        r"""Plot sea floor deformation dtopo as colormap
        
        """

        self.plot(axes=axes)



class Fault(object):

    r"""Basic object for representing a fault composed of many subfaults.

    """

    def __init__(self, path=None, subfaults=None):
        r"""Fault initialization routine.
        
        See :class:`Fault` for more info.

        """

        super(Fault, self).__init__()

        self.subfaults = subfaults
        
        if path is not None:
            # Read in file at path assuming it is a subfault specification
            self.read(path)
        elif subfaults is not None:
            if not isinstance(subfaults, list):
                raise ValueError("Input parameter subfaults must be a list.")
        else:
            raise ValueError("Must provide either a path to a subfault ",
                             "specification or a list of subfault objects.")


    def read(self, path, column_map, coordinate_specification="centroid",
                         rupture_type="static", t=[0.0, 0.5, 1.0]):
        r"""Read in subfault specification at *path*.

        Creates a list of subfaults from the subfault specification file at
        *path*.

        column_map = {"coordinates":(1,0), "depth":2, "slip":3, "rake":4, 
                          "strike":5, "dip":6}

        """

        # Read in rest of data
        data = numpy.loadtxt(path, skiprows=skiprows)

        self.subfaults = []
        for n in xrange(data.shape[0]):

            new_subfault = SubFault()
            new_subfault.coordinate_specification = coordinate_specification
            new_subfault.rupture_type = rupture_type
            new_subfault.t = numpy.array(t) 
            
            for (var, column) in column_map.iteritems():
                if isinstance(column, tuple) or isinstance(column, list):
                    for (k, index) in enumerate(column):
                        getattr(new_subfault, var)[k] = data[n, index]
                else:
                    setattr(new_subfault, var, data[n, column])

            self.subfaults.append(new_subfault)


    def write(self, path, topo_type=3):
        r"""Write out a dtopo file created by the combination of all subfaults

        """

        raise NotImplemented("Writing of a compound fault not implemneted yet.")

        if topo_type == 1:
            pass
        elif topo_type in [2,3]:
            pass
        else:
            raise ValueError("Unknown topography output type %s." % topo_type)


class UCSBFault(Fault):

    r"""Fault subclass for reading in subfault format models from UCSB

    Read in subfault format models produced by Chen Ji's group at UCSB,
    downloadable from:  

        http://www.geol.ucsb.edu/faculty/ji/big_earthquakes/home.html

    """

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
        found_data = False
        with open(path, 'r') as subfault_file:
            for line in subfault_file:
                result_dx = regexp_dx.search(line)
                result_dy = regexp_dy.search(line)
                result_nx = regexp_nx.search(line)
                result_ny = regexp_ny.search(line)

                if result_dx and result_dy:
                    i_dxdy = i
                    dx = float(result_dx.group('dx'))
                    dy = float(result_dy.group('dy'))
                    nx = int(result_nx.group('nx'))
                    ny = int(result_ny.group('ny'))
                    found_data = True
                    break
        if not found_data:
            raise ValueError("Could not find base fault characteristics in ",
                             "subfault specification file at %s." % path)

        column_map = {"coordinates":(1,0), "depth":2, "slip":3, "rake":4, 
                      "strike":5, "dip":6}

        super(UCSBFault, self).read(path, column_map=column_map, units={})

        # Set general data
        for fault in self.subfaults:
            fault.dimensions = [dx, dy]



class CSVFault(Fault):

    r"""Fault subclass for reading in CSV formatted files

    Assumes that the first row gives the column headings
    """

    def read(self, path):
        r"""Read in subfault specification at *path*.

        Creates a list of subfaults from the subfault specification file at
        *path*.

        """
        valid_column_labels = ()

        # Read header of file
        with open(path, 'r') as subfault_file:
            header_line = subfault_file.readline().split()
            column_map = {'coordinates':[None, None]}
            for (n,column_heading) in enumerate(header_line):
                if column_heading in valid_column_labels:
                    column_map[column_heading] = n
                elif column_heading == "longitude":
                    column_map['coordinates'][0] = n
                elif column_heading == "latitude":
                    column_map['coordinates'][1] = n

        super(CSVFault, self).read(path, column_map=column_map, units={})

