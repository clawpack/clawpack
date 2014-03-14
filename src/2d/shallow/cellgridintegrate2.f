c====================================================================
      subroutine cellgridintegrate(topoint,xim,xcell,xip,yjm,ycell,
     &           yjp,xlowtopo,ylowtopo,xhitopo,yhitopo,dxtopo,dytopo,
     &           mxtopo,mytopo,mtopo,i0topo,mtopoorder,
     &           mtopofiles,mtoposize,topo)
c=====================================================================

c *** Note: xcell and ycell are no longer needed -- should be removed.

      use topo_module, only: rectintegral, intersection

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
c     xintlo <= x <= xinthi, yintlo <= y <= yinthi

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
         call intersection(indicator,area,xmlo,xmhi,
     &       ymlo,ymhi,xim,xip,yjm,yjp,
     &       xlowtopo(mfid),xhitopo(mfid),ylowtopo(mfid),yhitopo(mfid))
         if (indicator.eq.1) then !cell overlaps grid
            if (area.eq.cellarea) then !cell is entirely in grid
               ! (should we check if they agree to some tolerance??)
c              !integrate surface and get out of here
                topoint = topoint + topointegral(xmlo,xmhi,ymlo,
     &              ymhi,xlowtopo(mfid),ylowtopo(mfid),dxtopo(mfid),
     &              dytopo(mfid),mxtopo(mfid),mytopo(mfid),
     &              topo(i0),im)
               return
            else
               go to 222
            endif
         endif
      enddo

      write(6,601) xim,xip,yjm,yjp
 601  format('*** Error, grid cell does not overlap any topo grid',/,
     &     '  xim = ',e24.14,'  xip = ',e24.14,      
     &   /,'  yjm = ',e24.14,'  yjp = ',e24.14)
      stop

 222  continue

      ! this grid cell intersects only topo grid m and perhaps coarser:
      call rectintegral(xim,xip,yjm,yjp,m,topoint)

      return
      end

c ## moved subroutine intersection to topo_module.f90
