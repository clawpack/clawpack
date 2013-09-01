from numpy import *
import os


def make_xyB_from_fgmax(outdir='_output'):
    # Create file for fixed grid with x,y,B.
    
    if not os.path.isdir(outdir):
        raise Exception("Missing directory: %s" % dirname)
        
    fname = outdir + '/fgmax.data'
    
    try:
        fid = open(fname)
    except:
        raise Exception("cannot open %s" % fname)

    for i in range(7):
        line = fid.readline()
    line = fid.readline().split()
    fid.close()
    mx = int(line[1])
    my = int(line[2])


    fname = outdir + '/fort.FG1.valuemax' 
    print "Reading %s ..." % fname
    try:
        d = loadtxt(fname)
    except:
        raise Exception("*** Cannot read file: %s" % fname)

    h = reshape(d[:,3],(mx,my),order='F')

    x = reshape(d[:,0],(mx,my),order='F')
    y = reshape(d[:,1],(mx,my),order='F')
    level = reshape(d[:,2].astype('int'),(mx,my),order='F')
    
    fname = outdir + '/fort.FG1.aux1' 
    print "Reading %s ..." % fname
    daux = loadtxt(fname)
    topo = []
    for i in range(2,9):
        topoi = reshape(daux[:,i],(mx,my),order='F')
        topoi = ma.masked_where(topoi < -1e50, topoi)
        topo.append(topoi)

    B = ma.masked_where(level==0, topo[0])  # level==0 ==> never updated
    levelmax = level.max()
    for i in range(levelmax):
        B = where(level==i+1, topo[i], B)

    fname = outdir + '/fgmax_xyB.txt'
    ofile = open(fname, 'w')
    for j in range(my):
        for i in range(mx):
            ofile.write(3*'%22.14e ' %(x[i,j],y[i,j],B[i,j]) + '\n')
    ofile.close()
    print "Created ",fname
    
        
    
if __name__=="__main__":
    import sys
    make_xyB_from_fgmax(*sys.argv[1:])
