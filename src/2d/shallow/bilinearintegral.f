c======================================================================
      function bilinearintegral(xim,xip,yjm,yjp,x1,x2,y1,y2,dxx,dyy,
     &                                              z11,z12,z21,z22)
c======================================================================

      implicit none

*     !i/o
      double precision xim,xip,yjm,yjp,x1,x2,y1,y2,dxx,dyy,
     &                                                z11,z12,z21,z22
      double precision bilinearintegral

*     !local
      double precision xlow,ylow,xhi,yhi,area,sumxi,sumeta,a,b,c,d

c#######################################################################
c
c         bilinearintegral integrates the bilinear with values z##
c         over the rectangular region xim <= x <= xip, and
c                                     yjm <= y <= yjp

c                                          written by David L George
c                                          Vancouver, WA April 2010
c#######################################################################



*     !integrate the portion of the bilinear intersected with the
c     !rectangular cell analytically.

c     !find limits of integral (this should already be true?)
      xlow = max(xim,x1)
      ylow = max(yjm,y1)
      xhi  = min(xip,x2)
      yhi =  min(yjp,y2)

*     !find the area of integration
      area = (yhi-ylow)*(xhi-xlow)
      sumxi = (xhi + xlow - 2.d0*x1)/dxx
      sumeta = (yhi + ylow - 2.d0*y1)/dyy

*     !find coefficients of bilinear a*xi + b*eta + c*xi*eta + d
      a = z21-z11
      b = z12-z11
      c = z22-z21-z12+z11
      d = z11

      bilinearintegral = (0.5d0*(a*sumxi + b*sumeta)
     &                     + 0.25d0*c*sumxi*sumeta + d)*area

      return

      end function


c======================================================================
      function bilinearintegral_s(xim,xip,yjm,yjp,x1,x2,y1,y2,dxx,dyy,
     &                                      z11,z12,z21,z22)
c======================================================================

      use geoclaw_module, only: r2d => rad2deg, d2r => deg2rad
      use geoclaw_module, only: Rearth => earth_radius

      implicit none

*     !i/o
      double precision xim,xip,yjm,yjp,x1,x2,y1,y2,dxx,dyy,
     &                                 z11,z12,z21,z22
      double precision bilinearintegral_s

*     !local
      double precision xlo,ylo,xhi,yhi,delx,dely,a,b,c,d
      double precision xdiffhi,xdifflo,ydiffhi,ydifflo,xdiff2
      double precision adsinint,cbsinint

c#######################################################################
c
c         bilinearintegral integrates the bilinear with values z##
c         over the rectangular region xim <= x <= xip, and
c                                     yjm <= y <= yjp
c         integration is actually done on the surface of a sphere

c                                          written by David L George
c                                          Vancouver, WA April 2010
c#######################################################################



*     !integrate the portion of the bilinear intersected with the
c     !rectangular cell analytically.
c     !find limits of integral (this should already be true?)
      xlo = max(xim,x1)
      ylo = max(yjm,y1)
      xhi  = min(xip,x2)
      yhi =  min(yjp,y2)
      delx = xhi - xlo
      dely = yhi - ylo

*     !find terms for the integration
      xdiffhi = xhi - x1
      xdifflo = xlo - x1
      ydiffhi = yhi - y1
      ydifflo = ylo - y1
      xdiff2 = 0.5d0*(xdiffhi**2 - xdifflo**2)

      cbsinint = (r2d*cos(d2r*yhi) + ydiffhi*sin(d2r*yhi))
     &            -(r2d*cos(d2r*ylo) + ydifflo*sin(d2r*ylo))

      adsinint = r2d*(sin(d2r*yhi) - sin(d2r*ylo))


*     !find coefficients of bilinear a*xi + b*eta + c*xi*eta + d
      a = (z21-z11)/dxx
      b = (z12-z11)/dyy
      c = (z22-z21-z12+z11)/(dxx*dyy)
      d = z11

      bilinearintegral_s = ((a*xdiff2 + d*delx)*adsinint
     &            + r2d*(c*xdiff2 + b*delx)*cbsinint)*(Rearth*d2r)**2

      return

      end function
