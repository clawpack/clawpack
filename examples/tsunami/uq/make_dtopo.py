
from pylab import *
import os
from pyclaw.geotools import topotools
import dtopotools


from pyclaw.geotools import okada
dtopo_fname = 'fault1.tt1'
dtopo_cfg = 'fault1.cfg'
if os.path.exists(dtopo_fname):
    print "*** Not regenerating dtopo file (already exists): %s" % dtopo_fname
else:
    print "Using Okada model to create %s " % dtopo_fname
    okada.builddynamicdeffile(dtopo_cfg, dtopo_cfg, dtopo_fname)

dtopo_fname = 'fault2.tt1'
dtopo_cfg = 'fault2.cfg'
if os.path.exists(dtopo_fname):
    print "*** Not regenerating dtopo file (already exists): %s" % dtopo_fname
else:
    print "Using Okada model to create %s " % dtopo_fname
    okada.builddynamicdeffile(dtopo_cfg, dtopo_cfg, dtopo_fname)


def make(d1,d2):
    """
    Create fault.tt1 file using fault1 with displacement d1 and
    fault2 with displacement d2.
    """
    inputfile = 'fault1.tt1'
    X,Y,dZ1 = dtopotools.read_dtopo(inputfile, deftype='dynamic',
               only_last=True)
    inputfile = 'fault2.tt1'
    X,Y,dZ2 = dtopotools.read_dtopo(inputfile, deftype='dynamic',
               only_last=True)

    dZ = d1*dZ1 + d2*dZ2
    outputfile = 'fault.tt1'
    fid = open(outputfile,'w')
    for t in [0,1]:
        for j in xrange(dZ.shape[0]):
            for i in xrange(dZ.shape[1]):
                fid.write('%012.6e %012.6e %012.6e %012.6e \n' %
                       (t,X[j,i],Y[j,i],t*dZ[j,i]))

    fid.close()

    print "Created %s using displacements d1 = %s, d2 = %s" \
            % (outputfile, d1, d2)
    return X,Y,dZ

def plot(X,Y,dZ):
    figure(200)
    clf()
    clines = linspace(-10,10,41)
    contour(X,Y,dZ,clines,colors='k')
    title('Contours of dtop at 0.5m increments')


if __name__=="__main__":
    try:
        d1 = float(sys.argv[1])
        d2 = float(sys.argv[2])
    except:
        print "*** Need to provide two float arguments d1 and d2"
        sys.exit()
    X,Y,dZ = make(d1,d2)
    plot(X,Y,dZ)
    savefig('fault.png')
    print "dtopo plot saved to fault.png"



