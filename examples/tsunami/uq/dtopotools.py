import numpy as np
import os

def read_dtopo(fname, deftype=None, only_last=True):
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

def plot_dynamic(X,Y,dZlist):
    import time
    import pylab
    pylab.figure(300)
    for dz in dZlist:
        pylab.clf()
        pylab.contour(X,Y,dz,np.linspace(-19.5,19.5,40),colors='b')
        pylab.draw()
        time.sleep(1)
        

