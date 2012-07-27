
from pyclaw.data import Data
import numpy as np 
import pylab

def read_header(fname):
    """
    Read a generic header of the form
      value  key
    """

    header = {}
    header_lines = ''
    f = open(fname,'r')

    for i in range(20):
        line = f.readline()
        tokens = line.split()
        #print "+++ tokens: ",tokens
        if ((len(tokens) != 2) or (type(tokens[1]) != str)):
            #print "+++ Done"
            break

        header_lines = header_lines + line
        try:
            try:
                value = int(tokens[0])
            except:
                value = float(tokens[0])
        except:
            raise ValueError("*** Cannot parse line ",line)

        key = tokens[1]
        header[key] = value

    header['header_lines'] = header_lines
    print header_lines
    f.close()
    return header


def read_dtopo(fname, topotype=1):
    header = read_header(fname)
    if (len(header['header_lines']) == 0) and (abs(topotype)>1):
        raise ValueError("*** topotype = %s but no header found" % topotype)
    if (len(header['header_lines']) > 0) and (abs(topotype)==1):
        raise ValueError("*** topotype = %s but header found" % topotype)
    header_length = header['header_lines'].count('\n')
    
    if topotype==1:
        d = np.loadtxt(fname, skiprows=header_length)
        ncol = d.shape[1]
        if ncol == 3:
            tcol = None
            xcol = d[:,0]
            ycol = d[:,1]
            dzcol = d[:,2]
        else:
            tcol = d[:,0]
            xcol = d[:,1]
            ycol = d[:,2]
            dzcol = d[:,3]
        xdiff = np.diff(xcol)
        inddiff = pylab.find(xdiff<0)
        mx = inddiff[0]+1
        mt = 1
        if tcol is not None:
            # check if there is more than one time in this file:
            tdiff = np.diff(tcol)
            inddiff = pylab.find(tdiff>0)
            if len(inddiff) > 0:
                mt = len(inddiff)+1
        my = len(xcol)/(mx*mt)
        if my != int(len(xcol)/float(mx*mt)):
            print "*** found mx = %s, mt = %s" % (mx,mt)
            raise ValueError("*** mx*mt does not divide length of data")
        print "Found mt = %s times with mx = %s, my = %s" % (mt,mx,my)
        if tcol is None:
            t = []
        else:
            t = tcol[0::mx*my]

        X = np.reshape(xcol[:mx*my],(my,mx))
        Y = np.reshape(ycol[:mx*my],(my,mx))
        dZ = np.reshape(dzcol,(mt,my,mx))

    elif topotype==3:
        try:
            my = header['nrows']
            mx = header['ncols']
            mt = header['mtimes']
            xll = header['xll']
            yll = header['yll']
            t0 = header['t0']
            dx = header['dx']
            dy = header['dy']
            dt = header['dt']
        except:
            raise ValueError("*** Missing information from header")
        if (dt==0) or (mt==1):
            t = np.array([t0])
        else:
            t = np.arange(t0,t0+mt*dt,dt)
        x = np.linspace(xll,xll + mx*dx, mx)
        y = np.linspace(yll,yll + my*dy, my)
        X,Y = np.meshgrid(x,y)
        X = X.T
        Y = Y.T

        fdata = np.loadtxt(fname, skiprows=header_length)
        dzcol = np.reshape(fdata,(mx*my*mt, 1))
        dZ = np.empty((mt,mx,my))
        for i in range(mt):
            dzi = np.reshape(dzcol[i*mx*my:(i+1)*mx*my], (my,mx))
            dzi = np.flipud(dzi)
            dZ[i,:,:] = dzi.T

    else:
        raise ValueError("*** Cannot read topotype = " % topotype)

    dtopo = Data()
    dtopo.t = t
    dtopo.X = X
    dtopo.Y = Y
    dtopo.dZ = dZ
    return dtopo
