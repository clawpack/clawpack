import numpy as np
import os

from pylab import *
from okada2 import okadamap
import copy

import matplotlib 

import re,sys

rad = np.pi/180.       # conversion factor from degrees to radians
rr = 6.378e6        # radius of earth
lat2meter = rr*rad  # conversion factor from degrees latitude to meters

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
        #data = loadtxt(fname, skiprows=skiprows, delimiter=delimiter)
        data = genfromtxt(fname, skiprows=skiprows, delimiter=delimiter,usecols=usecols)
    except:
        raise Exception("Unable to load file %s" % fname)
        
    try:
        ncols = data.shape[1]
    except:
        # if only one row in data file, convert to 2d array
        data = array([data])  
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
    

def read_subfault_model_csv(fname):
    """
    Read a subfault model from a csv file whose first row gives the column headings, 
    which should agree with names in valid_labels defined in read_subfault_model.
    """

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
    ang_strike = strike*rad
    ang_dip = dip*rad
    ang_rake = rake*rad
    x0 = longitude
    y0 = latitude

    if location == "top center":

        depth_top = depth
        depth_centroid = depth + 0.5*width*sin(ang_dip)
        depth_bottom = depth + width*sin(ang_dip)

        # Convert fault origin from top of fault plane to bottom:
        del_x = width*cos(ang_dip)*cos(ang_strike) / (lat2meter*cos(y0*rad))
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

        del_x = 0.5*width*cos(ang_dip)*cos(ang_strike) / (lat2meter*cos(y0*rad))
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
    dx2 = 0.5*length*sin(ang_strike) / (lat2meter*cos(y_bottom*rad))
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
    
    x=linspace(dtopo_params['xlower'],dtopo_params['xupper'],mx)
    y=linspace(dtopo_params['ylower'],dtopo_params['yupper'],my)
    dz = zeros((my,mx))
    
    print "Making Okada dz for each of %s subfaults" \
            % len(subfaults)

    for k,subfault in enumerate(subfaults):
            sys.stdout.write("%s.." % k)
            sys.stdout.flush()
            dz = dz + okadamap(subfault, x, y)
            #import pdb; pdb.set_trace()
            
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
    t = linspace(t0-10., t2+10., 1001)
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
    x=linspace(xlower,xupper,mx)
    y=linspace(ylower,yupper,my)
    times = linspace(t0,tfinal,ntimes)


    plot_rupture = False

    dZ = zeros((my,mx))
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
    for t in linspace(t0,tend,ntimes):
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


def test_make_dz_subfault_1():
    fname_subfaults = 'testdata/subfault-1.txt'
    fname_dtopo = 'testdata/subfault-1.tt1'
    dtopo_params = {}
    dtopo_params['mx'] = 101
    dtopo_params['my'] = 101
    dtopo_params['xlower'] = -1
    dtopo_params['xupper'] = 3
    dtopo_params['ylower'] = -1
    dtopo_params['yupper'] = 3
    subfaults = read_subfault_model_ucsb(fname_subfaults)
    x,y,dz = make_okada_dz(subfaults, dtopo_params)
    write_dz(fname_dtopo, x,y,dz)
    return x,y,dz,subfaults

def test_make_dz_subdivided_1():
    fname_subfaults = 'subdividedfault.csv'
    fname_dtopo = 'subdividedfault.tt1'
    dtopo_params = {}
    dtopo_params['mx'] = 101
    dtopo_params['my'] = 101
    dtopo_params['xlower'] = -2
    dtopo_params['xupper'] = 2
    dtopo_params['ylower'] = 48
    dtopo_params['yupper'] = 52
    subfaults = read_subfault_model_csv(fname_subfaults)
    x,y,dz = make_okada_dz(subfaults, dtopo_params)
    write_dz(fname_dtopo, x,y,dz)
    return x,y,dz,subfaults

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
        x=linspace(xlower,xupper,mx)
        y=linspace(ylower,yupper,my)
        times = linspace(t0, t0+(mt-1)*dt, mt)

        dZvals = loadtxt(fname, skiprows=9)
        dz_list = []
        for k,t in enumerate(times):
            dZk = reshape(dZvals[k*my:(k+1)*my, :], (my,mx))
            dZk = flipud(dZk)
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
            plot([x_centroid],[y_centroid],'go',label="Centroid")
            dxr = x_top - x_centroid
            dyr = y_top - y_centroid
            x_rake = x_centroid + cos(tau)*dxr - sin(tau)*dyr
            y_rake = y_centroid + sin(tau)*dxr + cos(tau)*dyr
            plot([x_rake,x_centroid],[y_rake,y_centroid],'g-',linewidth=2)
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

