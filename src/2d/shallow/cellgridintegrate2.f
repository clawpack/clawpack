c====================================================================
      subroutine cellgridintegrate(topoint,xim,xcell,xip,yjm,ycell,
     &           yjp,xlowtopo,ylowtopo,xhitopo,yhitopo,dxtopo,dytopo,
     &           mxtopo,mytopo,mtopo,i0topo,mtopoorder,
     &           mtopofiles,mtoposize,topo)
c=====================================================================


      implicit double precision (a-h,o-z)

      dimension topo(mtoposize)

      dimension xlowtopo(mtopofiles)
      dimension ylowtopo(mtopofiles)
      dimension xhitopo(mtopofiles)
      dimension yhitopo(mtopofiles)
      dimension dxtopo(mtopofiles)
      dimension dytopo(mtopofiles)

      dimension mxtopo(mtopofiles)
      dimension mytopo(mtopofiles)

      dimension mtopoorder(mtopofiles)
      dimension i0topo(mtopofiles)
      dimension mtopo(mtopofiles)


c     ##############################################################################
c     cellgridintegrate integrates a unique surface, over a rectangular cell
c     defined from data from multiple regular Cartesian grids
c     (using the finest data available in any region)
c
c     The rectangle has coords:
c     xim <= x <= xip, yjm <= y <= yjp, with center (x,y) = (xcell, ycell)
c
c     The intersection (with one particular grid has coords:
c     xintlo <= x <= xinthi, yintlo <= y <= yinthi, with center (x,y) = (xintc, yintc)

c     The _set_ version uses a recursive strategy using the formulas for
c     intersections of sets.

c     # initialize the integral of the surface
      topoint=0.d0

*     !determine the type of integration
      im = 1

c     # first see if the grid cell is entirely in a fine topofile
      do m = 1,mtopofiles
c        !look at topofiles, from fine to coarse
         mfid = mtopoorder(m)
         i0=i0topo(mfid)
c        !check for intersection of grid cell and this topofile
         cellarea = (xip-xim)*(yjp-yjm)
         call intersection(indicator,area,xmlo,xmc,xmhi,
     &       ymlo,ymc,ymhi,xim,xip,yjm,yjp,
     &       xlowtopo(mfid),xhitopo(mfid),ylowtopo(mfid),yhitopo(mfid))
         if (indicator.eq.1) then !cell overlaps grid
            if (area.eq.cellarea) then !cell is entirely in grid
               ! (should we check if they agree to some tolerance??)
c              !integrate surface and get out of here
                topoint = topoint + topointegral(
     &              xmlo,xmc,xmhi,ymlo,ymc,
     &              ymhi,xlowtopo(mfid),ylowtopo(mfid),dxtopo(mfid),
     &              dytopo(mfid),mxtopo(mfid),mytopo(mfid),
     &              topo(i0),im)
               return
            else
               go to 222
            endif
         endif
      enddo

 222  continue

      ! this grid cell intersects only topo grid m and perhaps coarser:
      call rectintegral(xim,xip,yjm,yjp,m,topoint)

      return
      end

c=======================================================================
      subroutine intersection(indicator,area,xintlo,xintc,xinthi,
     &      yintlo,yintc,yinthi,x1lo,x1hi,y1lo,y1hi,x2lo,x2hi,y2lo,y2hi)

c     find the intersection of two rectangles, return the intersection
c     and it's area, and indicator =1
c     if there is no intersection, indicator =0

      implicit none

c     !i/o integer
      integer indicator

c     !i/o doubles
      double precision area,xintlo,xintc,xinthi,yintlo,yintc,yinthi,
     &                 x1lo,x1hi,y1lo,y1hi,x2lo,x2hi,y2lo,y2hi


      xintlo=dmax1(x1lo,x2lo)
      xinthi=dmin1(x1hi,x2hi)
      yintlo=dmax1(y1lo,y2lo)
      yinthi=dmin1(y1hi,y2hi)

      xintc = 0.5d0*(xintlo+xinthi)
      yintc = 0.5d0*(yintlo+yinthi)

      if (xinthi.gt.xintlo.and.yinthi.gt.yintlo) then
         area = (xinthi-xintlo)*(yinthi-yintlo)
         indicator = 1
      else
         area = 0.d0
         indicator = 0
      endif

      return
      end











