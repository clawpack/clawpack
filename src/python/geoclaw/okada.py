#!/usr/bin/python
"""
okada
========================================================
This version of the Okada routines is deprecated but is
being kept around to test geoclaw/examples/tsunami/chile2010
Going forward, dtopotools.okadamap should be used instead.
========================================================

    python module for getting Okada maps.

    Okada model is a mapping from several fault parameters
    to a surface deformation.
    See Okada 1985, or Okada 1992, Bull. Seism. Soc. Am.

    some routines adapted from fortran routines written by
    Xiaoming Wang.

"""
from __future__ import absolute_import
from __future__ import print_function
import pdb
import numpy
#import scipy
from numpy import *
#from scipy import *
from matplotlib import *
import matplotlib.mlab as mlab
import matplotlib.pyplot as pyplot
import os
import string
from .datatools import *
from six.moves import range

#=================================================================================
def builddeffile (okadaparamfile,faultparamfile,outfile):

    faultparams=getokadaparams(okadaparamfile)
    faultparams.update(getfaultparams(faultparamfile))

    fid=open(outfile,'w')

    X=linspace(faultparams['xlower'],faultparams['xupper'],faultparams['mx'])
    Y=linspace(faultparams['ylower'],faultparams['yupper'],faultparams['my'])

    dZ=okadamap(faultparams,X,Y)
    ind=fixdata.findbadindices(dZ)
    if ind:
        dZ=fixdata.fillbaddata(dZ,ind)

    dZ = filtermask(dZ,faultparams)
    #pdb.set_trace()
    for jj in range(faultparams['my']):
        j=-1-jj
        for i in range(faultparams['mx']) :
            fid.write('%012.6e %012.6e %012.6e \n' % (X[i],Y[j],dZ[j,i]))

    fid.close()
    return

#=================================================================================
def builddynamicdeffile (okadaparamfile,faultparamfile,outfile,t0=0.0, tend=1.0, nt = 2):

    faultparams=getokadaparams(okadaparamfile)
    faultparams.update(getfaultparams(faultparamfile))

    fid=open(outfile,'w')

    X=linspace(faultparams['xlower'],faultparams['xupper'],faultparams['mx'])
    Y=linspace(faultparams['ylower'],faultparams['yupper'],faultparams['my'])

    T=linspace(t0,tend,nt)

    dZ=okadamap(faultparams,X,Y)
    ind=fixdata.findbadindices(dZ)
    if ind:
        dZ=fixdata.fillbaddata(dZ,ind)

    dZ = filtermask(dZ,faultparams)
    #pdb.set_trace()
    for it in T:
        alpha=(it-t0)/(tend-t0)
        for jj in range(faultparams['my']):
            j=-1-jj
            for i in range(faultparams['mx']) :
                fid.write('%012.6e %012.6e %012.6e %012.6e \n' % (it,X[i],Y[j],alpha*dZ[j,i]))

    fid.close()
    return

#=================================================================================
def getokadaparams (infile):

    """
    obtain parameters necessary for okada map from config file: infile

        file format:
        parameter names and values should appear on the same single line seperated by a space
    """

    keylist=["Focal_Depth","Fault_Length","Fault_Width","Dislocation","Strike_Direction", \
             "Dip_Angle","Slip_Angle","Epicenter_Latitude","Epicenter_Longitude"]

    okadaparams={}
    fid=open(infile,'r')
    keyleft=len(keylist)
    while keyleft> 0 :
        line=string.split(fid.readline())
        if line:
            if line[0] in keylist:
                okadaparams[line[0]]=float(line[1])
                keyleft=keyleft-1
            if line[1] in keylist:
                okadaparams[line[1]]=float(line[0])
                keyleft=keyleft-1

    for key in keylist :
        if not key in okadaparams:
            print(('ERROR: parameters for okada fault not fully specified in %s' % (infile)))
            exit

    fid.close()
    return okadaparams
    #end getokadaparams=============================================================

#===================================================================================
def getfaultparams (infile):

    """
    obtain params from a file that specify a fault grid from infile
    params are xlower,ylower,dx,dy,mx,my, OR
    xlower,ylower,xupper,yupper,mx,my

    file format:
        parameter names and values should appear on the same single line seperated by a space
    """

    keylist=["xlower","ylower","xupper","yupper","dx","dy","mx","my"]

    faultgridparams={}
    fid=open(infile,'r')
    keyleft=len(keylist)-2
    while keyleft> 0 :
        line=string.split(fid.readline())
        if line:
            if line[0] in keylist:
                faultgridparams[line[0]]=float(line[1])
                keyleft=keyleft-1
            if line[1] in keylist:
                faultgridparams[line[1]]=float(line[0])
                keyleft=keyleft-1

    faultgridparams['mx'] = int(faultgridparams['mx'])
    faultgridparams['my'] = int(faultgridparams['my'])

    if ('dx' in faultgridparams)& ('dy' in faultgridparams):
        faultgridparams['xupper'] = faultgridparams['xlower'] + faultgridparams['dx']*(faultgridparams['mx']-1)
        faultgridparams['yupper'] = faultgridparams['ylower'] + faultgridparams['dy']*(faultgridparams['my']-1)
    elif ('xupper' in faultgridparams)&('yupper' in faultgridparams):
        faultgridparams['dx'] = (faultgridparams['xupper']-faultgridparams['xlower'])/(faultgridparams['mx']-1)
        faultgridparams['dy'] = (faultgridparams['yupper']-faultgridparams['ylower'])/(faultgridparams['my']-1)
    else:
        print(('ERROR: parameters for fault grid not fully specified in %s' % (infile)))
        exit

    for key in keylist :
        if not key in faultgridparams:
            print(('ERROR: parameters for fault grid not fully specified in %s' % (infile)))
            exit

    fid.close()
    return faultgridparams
    #end getfaultparams===========================================================


#==================================================================================
def  okadamap(okadaparams,X,Y):

    """
    create displacement matrix dZ for a surface displacement
    over gridded region defined by X,Y, vectors of length nx,ny
    given okadaparams
    """

    xo=X[0]
    yo=Y[0]
    nx=len(X)
    ny=len(Y)

    zero = 0.0
    osixty = 0.016666666667
    rad = 0.01745329252
    rr = 6.378e6

    hh =  okadaparams["Focal_Depth"]
    l  =  okadaparams["Fault_Length"]
    w  =  okadaparams["Fault_Width"]
    d  =  okadaparams["Dislocation"]
    th =  okadaparams["Strike_Direction"]
    dl =  okadaparams["Dip_Angle"]
    rd =  okadaparams["Slip_Angle"]
#    yo =  okadaparams["Domain_Latitude"]
#    xo =  okadaparams["Domain_Longitude"]
    y0 =  okadaparams["Epicenter_Latitude"]
    x0 =  okadaparams["Epicenter_Longitude"]

    ang_l = rad*dl
    ang_r = rad*rd
    ang_t = rad*th
    halfl = 0.5*l

    #Calculate focal depth used for Okada's model
    hh = hh+w*sin(ang_l)
    #displacement due to different epicenter definition
    del_x = w*cos(ang_l)*cos(ang_t)
    del_y = w*cos(ang_l)*sin(ang_t)
    xl = rr*cos(rad*yo)*(x0-xo)*rad + del_x
    yl = rr*(y0-yo)*rad - del_y
    ds = d*cos(ang_r)
    dd = d*sin(ang_r)
    sn = sin(ang_l)
    cs = cos(ang_l)

    dZ=numpy.eye(ny,nx)

    # Vectorized by RJL, 12/10:
    x,y = meshgrid(X,Y)

    #!-----added by Xiaoming Wang, Nov 11 2006----
    xl = rr*cos(rad*y)*(x0-xo)*rad + del_x  # TAL - Fixed sign, 9/07
    #!---------------------------------------------
    yy = rr*(y-yo)*rad
    xx = rr*cos(rad*y)*(x-xo)*rad

    x1 = (xx-xl)*sin(ang_t)+(yy-yl)*cos(ang_t)
    x2 = (xx-xl)*cos(ang_t)-(yy-yl)*sin(ang_t)
    x2 = -x2
    x3 = zero
    p = x2*cs+hh*sn

    f1=strike_slip (x1,x2,x3,x1+halfl,p,ang_l,hh)
    f2=strike_slip (x1,x2,x3,x1+halfl,p-w,ang_l,hh)
    f3=strike_slip (x1,x2,x3,x1-halfl,p,ang_l,hh)
    f4=strike_slip (x1,x2,x3,x1-halfl,p-w,ang_l,hh)
    g1=dip_slip (x1,x2,x3,x1+halfl,p,ang_l,hh)
    g2=dip_slip (x1,x2,x3,x1+halfl,p-w,ang_l,hh)
    g3=dip_slip (x1,x2,x3,x1-halfl,p,ang_l,hh)
    g4=dip_slip (x1,x2,x3,x1-halfl,p-w,ang_l,hh)


    us = (f1-f2-f3+f4)*ds
    ud = (g1-g2-g3+g4)*dd

    dZ = (us+ud)

    return dZ
    #=========================================================================


#==============================================================================
def strike_slip (x1,x2,x3,y1,y2,dp,dd):
    """
    !.....Used for Okada's model
    !.. ..Methods from Yoshimitsu Okada (1985)
    !-----------------------------------------------------------------------
    """
    sn = sin(dp)
    cs = cos(dp)
    p = x2*cs + dd*sn
    q = x2*sn - dd*cs
    d_bar = y2*sn - q*cs
    r = sqrt(y1**2 + y2**2 + q**2)
    xx = sqrt(y1**2 + q**2)
    a4 = 0.5*1/cs*(log(r+d_bar) - sn*log(r+y2))
    f = -(d_bar*q/r/(r+y2) + q*sn/(r+y2) + a4*sn)/(2.0*3.14159)

    return f
    #==========================================================================


#==============================================================================
def dip_slip (x1,x2,x3,y1,y2,dp,dd):
    """
    !.....Based on Okada's paper (1985)
    !.....Added by Xiaoming Wang
    !-----------------------------------------------------------------------
    """
    sn = sin(dp)
    cs = cos(dp)

    p = x2*cs + dd*sn
    q = x2*sn - dd*cs
    d_bar = y2*sn - q*cs;
    r = sqrt(y1**2 + y2**2 + q**2)
    xx = sqrt(y1**2 + q**2)
    a5 = 0.5*2/cs*arctan((y2*(xx+q*cs)+xx*(r+xx)*sn)/y1/(r+xx)/cs)
    f = -(d_bar*q/r/(r+y1) + sn*arctan(y1*y2/q/r) - a5*sn*cs)/(2.0*3.14159)

    return f
    #===========================================================================


#================================================================================
def filtermask (dZ,faultparams):
    """
    borrowed from code written by Xiaoming Wang and Tom Logan at ARSC

    !.....Filter the deformation using a circular mask centered
    !.....at the epicenter using a calculated radius
    !.....Removes small numerical artifacts away from the epicenter
    """
    filterindices=[]

    osixty = 0.016666666667
    rad = 0.01745329252
    rr = 6.378e6

    xo = faultparams['xlower']
    yo = faultparams['ylower']
    nx = faultparams['mx']
    ny = faultparams['my']
    spacing = faultparams['dx']

    x0 = faultparams['Epicenter_Longitude']
    y0 = faultparams['Epicenter_Latitude']
    l =  faultparams['Fault_Length']
    w =  faultparams['Fault_Width']
    dl = faultparams['Dip_Angle']


    ang_l = rad*dl # convert degree to radian

    #!-- epicenter in pixels -----------
    ypix = (y0-yo)/spacing
    xpix = (x0-xo)/spacing

    #!-- conversion from meters to pixels ---
    tmpd=spacing*rad
    xdist = tmpd*rr

    #!-- size of the fault in pixels --------
    npix_x = l/xdist
    npix_y = w/xdist

    #!-- set the range (radius) of the filter circle --------
    #!----- for small dip angles, use the length and width --
    #!----- for larger dip angles, use only the length ------

    if dl<30.0:
        drange = 1.5 * cos(ang_l)*sqrt(npix_x*npix_x+npix_y*npix_y)
    else:
        drange = 1.2 * npix_x

    print(("Filtering deformation using a circle of radius %s" % (drange)))

    #!-- Create the filtering mask ----------
    for i in range(nx):
        for j in range(ny) :
            dist = sqrt((i+1-xpix)**2+(j+1-ypix)**2)
            if dist > drange :
                filterindices.append((j,i))

    #!-- apply the filter to the actual deformation ------
    dZ = fixdata.filterdata(dZ,filterindices,radius=2)

    return dZ
