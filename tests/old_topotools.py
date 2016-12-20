r"""Basic functions from the old version of topotools for testing purposes."""

from __future__ import absolute_import
from __future__ import print_function
import numpy
import re

import nose
from six.moves import range

def topo1writer(outfile, topo, xlower, xupper, ylower, yupper, nxpoints,
                             nypoints):
    fout=open(outfile, 'w')
    dx = float(xupper-xlower)/(nxpoints-1)
    dy = float(yupper-ylower)/(nypoints-1)

    x = numpy.linspace(xlower,xupper,nxpoints)
    y = numpy.linspace(ylower,yupper,nypoints)
    X,Y = numpy.meshgrid(x,y)
    Z = topo(X,Y).T
    

    for jj in range(0,nypoints):
        y = yupper - jj*dy
        for i in range(0,nxpoints):
            x =  xlower + i*dx
            j = nypoints - 1 - jj
            z = Z[i,j]
            fout.write("%22.15e  %22.15e  %22.15e\n" % (x,y,z)) 

    fout.close

def topoheaderread(inputfile, closefile=True):
    topoheader={'ncols':0,'nrows':0,'xll':0.0,'yll':0.0,'cellsize':0.0,'nodata_value':0}
    keylist=list(topoheader.keys())

    keymap = {'ncols':'ncols','nrows':'nrows','xll':'xll','yll':'yll','cellsize':'cellsize','nodata_value':'nodata_value', \
    'xllcenter':'xll','yllcenter':'yll','xllcorner':'xll','yllcorner':'yll'}

    fid=open(inputfile,'r')
    keyleft=len(keylist)
    while keyleft> 0 :
        line=fid.readline().split()
        if line:
            if line[0].lower() in list(keymap.keys()):
                topoheader[keymap[line[0].lower()]]= convertd2e(line[1])
                keyleft=keyleft-1
            if line[1].lower() in list(keymap.keys()):
                topoheader[keymap[line[1].lower()]]= convertd2e(line[0])
                keyleft=keyleft-1

    #check if passes convert strings values to numeric
    for key in keylist :
        if not key in topoheader:
            print(('ERROR: topoheader not fully specified in %s' % (inputfile)))
            exit
        else:
            if '.' in topoheader[key] or 'nan' in topoheader[key].lower() or 'e' in topoheader[key].lower():
                topoheader[key]=float(topoheader[key])
            else:
                topoheader[key]=int(topoheader[key])


    if closefile:
        fid.close()
        return topoheader
    else:
        return (fid,topoheader)


def datafile2array(datafile=" ",sep=None, dtype="float",skiplines=0, \
skipfirstcols=0, skiplastcols=0):
    fid=open(datafile)
    data=fid.readlines()
    fid.close()

    dataarray=[]
    for row in range(skiplines,len(data)):
        data[row]=convertd2e(data[row])
        data[row]=data[row].split(sep)
        if data[row]!=[]:
            if dtype!=" ":
                for col in range(skipfirstcols,len(data[row])-skiplastcols) :
                    if dtype=="float":
                        data[row][col]=float(data[row][col])
                    elif dtype=="int":
                        data[row][col]=int(data[row][col])
            if dataarray!=[]: 
                if len(data[row])-skipfirstcols-skiplastcols==len(dataarray[0]):
                    dataarray.append(data[row][skipfirstcols:len(data[row])-skiplastcols])
            else:
                dataarray.append(data[row][skipfirstcols:len(data[row])-skiplastcols])
    
    dataarray=numpy.array(dataarray)
    return dataarray


def topofile2griddata(inputfile, topotype=2):
    
    try:
        import pylab
    except ImportError:
        raise nose.SkipTest("Skipping test since matplotlib was not found.")

    if topotype>1:
        (fin,topoheader)=topoheaderread(inputfile,closefile=False)
        zdata=fin.readlines()
        fin.close()
        for row in range(len(zdata)):
            zdata[row]=convertd2e(zdata[row])
            zdata[row]=zdata[row].split()
            for col in range(len(zdata[row])) :
                zdata[row][col]=float(zdata[row][col])

        Z=numpy.array(zdata)
        Z=numpy.reshape(Z,(topoheader['nrows'],topoheader['ncols']))

        xlower=topoheader['xll']
        xupper=xlower+ topoheader['cellsize']*(topoheader['ncols']-1)

        ylower = topoheader['yll']
        yupper = ylower+ topoheader['cellsize']*(topoheader['nrows']-1)

        x=numpy.linspace(xlower,xupper,topoheader['ncols'])
        y=numpy.linspace(ylower,yupper,topoheader['nrows'])
        [X,Y]=numpy.meshgrid(x,y)
        Y=numpy.flipud(Y)
    else:
        a=datafile2array(inputfile)
        xdiff=numpy.diff(a[:,0])
        #inddiff=pylab.find(xdiff<0)  
        inddiff = numpy.nonzero(xdiff<0)[0]  # rewrite above line without find
        xlength=inddiff[0]+1
        ylength=len(a[:,0])/xlength
        x=a[:,0]
        y=a[:,1]
        z=a[:,2]

        X=numpy.reshape(x,(ylength,xlength))
        Y=numpy.reshape(y,(ylength,xlength))
        Z=numpy.reshape(z,(ylength,xlength))

    return X,Y,Z

def convertd2e (numberstring=" "):
    Dd=re.compile("[Dd]")
    newstring=Dd.sub("e",numberstring)
    return newstring
