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

c     ! grid cell integral must come from more than one topofile
      do m = mtopofiles,1,-1
c        ! look for topo files overlapping this grid cell.
c        ! By decending mtopoorder from mtopoorder(mtopofiles) to
c        ! mtopoorder(1) we are decending from the coarsest to the finest topofile.
c        ! by the end of the nested loops, the integral of this topofile
c        ! will be taken over regions not covered by finer topofiles
c        ! nested loops only account for a maximum of 4 intersecting topo files
c        ! if 5 topofiles intersect, then an error is sent.

         mfid = mtopoorder(m)

c        ! note that mfid indicates the topofile number or id for the "m'th" coarsest topofile
         ! within this do loop, we are always integrating mfid
c        ! the nested do loops find intersections with other topofiles
c        ! to determin the areas, but for integration of mfid only!

c        !check for intersection of grid cell and this topofile
         call intersection(indicator,area,xmlo,xmc,xmhi,
     &       ymlo,ymc,ymhi,xim,xip,yjm,yjp,
     &       xlowtopo(mfid),xhitopo(mfid),ylowtopo(mfid),yhitopo(mfid))

         if (indicator.eq.1) then
c           ! cell overlaps the file
            i0=i0topo(mfid)
c           ! integrate surface over intersection of grid and cell
            topoint = topoint + topointegral(
     &              xmlo,xmc,xmhi,ymlo,ymc,
     &              ymhi,xlowtopo(mfid),ylowtopo(mfid),dxtopo(mfid),
     &              dytopo(mfid),mxtopo(mfid),mytopo(mfid),
     &              topo(i0),im)

c           ! loop through grids finer than this one and subtract any
c           ! integrals over intersections
            do mm = m-1,1,-1
               mfidint = mtopoorder(mm)
c              ! check for 2nd intersection
               call intersection(indicator,area,xmmlo,xmmc,xmmhi,
     &                 ymmlo,ymmc,ymmhi,xmlo,xmhi,ymlo,ymhi,
     &                 xlowtopo(mfidint),xhitopo(mfidint),
     &                 ylowtopo(mfidint),yhitopo(mfidint))

               if (indicator.eq.1) then
c                 ! get rid of coarser integral
                  topoint = topoint - topointegral(
     &              xmmlo,xmmc,xmmhi,ymmlo,ymmc,
     &              ymmhi,xlowtopo(mfid),ylowtopo(mfid),dxtopo(mfid),
     &              dytopo(mfid),mxtopo(mfid),mytopo(mfid),
     &              topo(i0),im)
c                  -----------------------------------------------------
c                 ------------------------------------------------------------
c                 ! loop through grids finer than this one and add back any
c                 ! integrals over intersections that will later get subtracted again
                  do mmm = mm-1,1,-1
                     mfidint = mtopoorder(mmm)
c                    ! check for 3rd intersection
                   call intersection(indicator,area,xmmmlo,xmmmc,xmmmhi,
     &                 ymmmlo,ymmmc,ymmmhi,xmmlo,xmmhi,ymmlo,ymmhi,
     &                 xlowtopo(mfidint),xhitopo(mfidint),
     &                 ylowtopo(mfidint),yhitopo(mfidint))

                     if (indicator.eq.1) then
c                       ! add back
                        topoint = topoint + topointegral(
     &                  xmmmlo,xmmmc,xmmmhi,ymmmlo,ymmmc,ymmmhi,
     &                  xlowtopo(mfid),ylowtopo(mfid),dxtopo(mfid),
     &                  dytopo(mfid),mxtopo(mfid),mytopo(mfid),
     &                                             topo(i0),im)

                        do m4 = mmm-1,1,-1
                           mfidint = mtopoorder(m4)
c                          ! check for 4th intersection
                           call intersection(indicator,area,xm4lo,
     &                        xm4c,xm4hi,ym4lo,ym4c,ym4hi,xmmmlo,
     &                        xmmmhi,ymmmlo,ymmmhi,
     &                        xlowtopo(mfidint),xhitopo(mfidint),
     &                        ylowtopo(mfidint),yhitopo(mfidint))

                           if (indicator.eq.1) then
c                          ! add back
                              topoint = topoint - topointegral(
     &                        xm4lo,xm4c,xm4hi,ym4lo,ym4c,ym4hi,
     &                        xlowtopo(mfid),ylowtopo(mfid),
     &                        dxtopo(mfid),dytopo(mfid),mxtopo(mfid),
     &                        mytopo(mfid),topo(i0),im)

                                 do m5 = m4-1,1,-1
                                    mfidint = mtopoorder(m5)
c                                   ! check for 5th intersection
                                    call intersection(indicator,area,
     &                                 xm5lo,xm5c,xm5hi,ym5lo,ym5c,
     &                                 ym5hi,xm4lo,xm4hi,ym4lo,
     &                                 ym4hi,xlowtopo(mfidint),
     &                                 xhitopo(mfidint),
     &                                 ylowtopo(mfidint),
     &                                 yhitopo(mfidint))

                                    if (indicator.eq.1) then
                                       write(*,*) 'CELLGRIDINTEGRATE:'
            write(*,*) 'ERROR: 5 NESTED TOPOGRIDS. MAXIMUM OF 4 ALLOWED'
                                    endif
                                 enddo
                           endif
                        enddo
                     endif
                  enddo

               endif
            enddo
c           ------------------------------------------------------------
         endif

      enddo

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











