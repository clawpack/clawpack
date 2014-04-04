
c===========================================================================
       function topointegral(xim,xip,yjm,yjp,
     &                   xxlow,yylow,dxx,dyy,mxx,myy,zz,intmethod)
c===========================================================================

c###########################################################################
c     topointegral integrates a surface over a rectangular region
c     that is the intersection with a Cartesion grid (of topography data)
c     the surface integrated is defined by a piecewise bilinear through the
c     nodes of the Cartesian grid.

c      The rectangular intersection has coords:
c      xim <= x <= xip, yjm <= y <= yjp
c
c      The Cartesian grid has coords:
c      xxlow <= x <= xxhi, yylow <= y <= yyhi, with grid cell size dxx by dyy
c      and mxx by myy cells.
c
c                                                written by David L. George
c                                                Seattle, WA 7/16/08
c###########################################################################

      use geoclaw_module

      implicit double precision (a-h,o-z)
      double precision zz(1:mxx,1:myy)

c     # initialize:
      theintegral = 0.d0

      xxhi=xxlow+(mxx-1)*dxx
      yyhi=yylow+(myy-1)*dyy

c========TEST FOR SMALL ROUNDING ERROR==========
      if ((xim-xxlow).lt.0.d0.or.(xip-xxhi).gt.0.d0) then
         xim=dmax1(xxlow,xim)
         xip=dmin1(xxhi,xip)
         endif

      if ((yjm-yylow).lt.0.d0.or.(yjp-yyhi).gt.0.d0) then
         yjm=dmax1(yylow,yjm)
         yjp=dmin1(yyhi,yjp)
         endif
c=========================================

      dx=xip-xim
      dy=yjp-yjm

c=============INTEGRATE PIECEWISE BILINEAR OVER RECTANGULAR REGION====
      if (intmethod.eq.1) then !use bilinear method

c         don't waste time looping through the entire grid
c         just find indices that include the rectangular region

         djjstart=(yjm-yylow)/dyy
         jjstart=idint(djjstart)+1

         diistart=(xim-xxlow)/dxx
         iistart=idint(diistart)+1

         diiend=(xip-xxlow)/dxx
         iiend=ceiling(diiend) + 1

         djjend=(yjp-yylow)/dyy
         jjend=ceiling(djjend)+1

         iistart=max(iistart,1)
         jjstart=max(jjstart,1)
         iiend=min(mxx,iiend)
         jjend=min(myy,jjend)


         do jj=jjstart,jjend-1
            y1=yylow + (jj-1.d0)*dyy
            y2=yylow + (jj)*dyy
c           # the array zz is indexed from north to south: jjz is the actual index
c           # of interest in the array zz
            jjz1= myy-jj+1
            jjz2= jjz1-1

            do ii=iistart,iiend-1
               x1=xxlow + (ii-1.d0)*dxx
               x2=xxlow + (ii)*dxx

               z11 = zz(ii,jjz1)
               z12 = zz(ii,jjz2)
               z21 = zz(ii+1,jjz1)
               z22 = zz(ii+1,jjz2)

               if (coordinate_system.eq.1) then !cartesian rectangle
                  theintegral = theintegral + bilinearintegral(
     &                                         xim,xip,yjm,yjp,
     &                                         x1,x2,y1,y2,
     &                                         dxx,dyy,
     &                                         z11,z12,z21,z22)
                elseif (coordinate_system.eq.2) then !integrate on surface of sphere
                  theintegral = theintegral + bilinearintegral_s(
     &                                         xim,xip,yjm,yjp,
     &                                         x1,x2,y1,y2,
     &                                         dxx,dyy,
     &                                         z11,z12,z21,z22)
                else
                  write(*,*)  'TOPOINTEGRAL: coordinate_system error'
                  endif
               enddo
            enddo

      else
         write(*,*) 'TOPOINTEGRAL: only intmethod = 1,2 is supported'
         endif

      topointegral= theintegral
      return
      end function
