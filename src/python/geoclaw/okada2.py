"""
Module for computing seafloor deformation using Okada model.

Okada model is a mapping from several fault parameters
to a surface deformation.
See Okada 1985, or Okada 1992, Bull. Seism. Soc. Am.

Some routines adapted from fortran routines written by
Xiaoming Wang.

okadamap function riginally written in Python by Dave George in okada.py.
Rewritten and made more flexible by Randy LeVeque:
Location can be specified as "top center" or "centroid".

The main function is okadamap(okadaparams,X,Y).

"""

from pylab import *
import os
import string
from datatools import fixdata

poisson = 0.25  # Poisson ratio


#=============================================================================

def okadamap(okadaparams,X,Y):

    """
    create displacement matrix dZ for a surface displacement
    over gridded region defined by X,Y, vectors of length nx,ny
    given okadaparams
    """

    rad = pi/180.       # conversion factor from degrees to radians
    rr = 6.378e6        # radius of earth
    lat2meter = rr*rad  # conversion factor from degrees latitude to meters

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

    ang_dip = rad*dl
    ang_slip = rad*rd
    ang_strike = rad*th
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
        del_x = w*cos(ang_dip)*cos(ang_strike) / (lat2meter*cos(y0*rad))
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
        del_x = 0.5*w*cos(ang_dip)*cos(ang_strike) / (lat2meter*cos(y0*rad))
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
    dx2 = 0.5*L*sin(ang_strike) / (lat2meter*cos(y_bottom*rad))
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


    x,y = meshgrid(X,Y)

    # Convert distance from (x,y) to (x_bottom,y_bottom) from degrees to meters:
    xx = lat2meter*cos(rad*y)*(x-x_bottom)   
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
        contour(x,y,dZ,linspace(-8,8,17),colors='k')

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
    a4 = 2.0*poisson/cs*(log(r+d_bar) - sn*log(r+y2))
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
    a5 = 4.*poisson/cs*arctan((y2*(xx+q*cs)+xx*(r+xx)*sn)/y1/(r+xx)/cs)
    f = -(d_bar*q/r/(r+y1) + sn*arctan(y1*y2/q/r) - a5*sn*cs)/(2.0*3.14159)

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
    rad = 0.01745329252
    rr = 6.378e6

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


    ang_dip = rad*dl # convert degree to radian

    #!-- fault origin in pixels -----------
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


#============================================================================

def filterdata (Z,filterinds,radius=1):
    """
    filter data in array z, at indice tuples in list filterinds
    by averaging surrounding data, ball with radius=radius in inf-norm
    acts as a low-band pass filter and removes oscillatory data
    """

    m=shape(Z)[0]
    n=shape(Z)[1]

    for ind in filterinds :
        i=ind[0]
        j=ind[1]
        r=radius

        irange=range(max(0,i-r),min(i+r+1,m))
        jrange=range(max(0,j-r),min(j+r+1,n))
        summands=0
        sum=0.
        for ii in irange:
            for jj in jrange:
                ballind=(ii,jj)
                sum = sum + Z[ballind[0],ballind[1]]
                summands=summands+1
        if summands >0 : 
            Z[ind[0],ind[1]] = sum/summands

    return Z

