c
c -------------------------------------------------------------------
      subroutine flag2refine(mx,my,mbc,meqn,maux,xlower,ylower,dx,dy,
     &                 t,level,tolsp,q,aux,amrflags,DONTFLAG,DOFLAG)
c -------------------------------------------------------------------

c
c ::::::::::::::::::::: flag2refine ::::::::::::::::::::::::::::::::::
c
c User routine to control flagging of points for refinement.
c
c 
c
c The logical function allowflag(x,y,t) is called to 
c check whether further refinement at this level is allowed in this cell
c at this time.  
c
c    q   = grid values including ghost cells (bndry vals at specified
c          time have already been set, so can use ghost cell values too)
c
c  aux   = aux array on this grid patch
c
c amrflags  = array to be flagged with either the value
c             DONTFLAG (no refinement needed)  or
c             DOFLAG   (refinement desired)    
c

c Specific for GeoClaw for flood type problems: ie: dambreaks, river flooding etc.
c  
c This particular version flag2refine_geo_flood.f differs from the routine used
c for tsunami modeling.

c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      use geoclaw_module
      use topo_module

      implicit double precision (a-h, o-z)

      dimension   q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      dimension   aux(1-mbc:mx+mbc,1-mbc:my+mbc,maux)
      dimension   amrflags(1-mbc:mx+mbc,1-mbc:my+mbc)
      logical     allowflag
      external  allowflag
      logical shoreregion,wave,shoreline,land,river,ocean

      include 'regions.i'
      include 'qinit.i'

c     # loop over interior points on this grid:

      do 200 j = 1,my
        y = ylower + mbc*dy + (j-0.5d0)*dy
        y1 = ylower + mbc*dy + (j-1)*dy
        y2 = ylower + mbc*dy + j*dy
        do 100 i = 1,mx
          x = xlower + mbc*dx + (i-0.5d0)*dx
          x1 = xlower + mbc*dx + (i-1)*dx
          x2 = xlower + mbc*dx + i*dx
c         # (i,j) grid cell is [x1,x2] x [y1,y2].

c         # default for each point is not to flag unless some condition
c         # below is satisfied:

          amrflags(i,j) = DONTFLAG
      

          do 30 m=1,mtopofiles
c           # check to see if refinement is forced in any topo file region:
            if (level .lt. minleveltopo(m) .and.
     &          t.ge.tlowtopo(m) .and. t.le.thitopo(m)) then
              xlow = xlowtopo(m)
              xhi = xhitopo(m)
              ylow = ylowtopo(m)
              yhi = yhitopo(m)
                 if (x2.gt.xlow.and.x1.lt.xhi.and.
     &               y2.gt.ylow.and.y1.lt.yhi) then
                    amrflags(i,j) = DOFLAG
                    go to 100 !# flagged, so no need to check anything else
                    endif
              endif
   30       continue

           do 40 m=1,mregions
c           # check to see if refinement is forced in any other region:
            if (level .lt. minlevelregion(m) .and.
     &          t.ge.tlowregion(m) .and. t.le.thiregion(m)) then
              xlow = xlowregion(m)
              xhi = xhiregion(m)
              ylow = ylowregion(m)
              yhi = yhiregion(m)
                 if (x2.gt.xlow.and.x1.lt.xhi.and.
     &               y2.gt.ylow.and.y1.lt.yhi) then
                    amrflags(i,j) = DOFLAG
                    go to 100 !# flagged, so no need to check anything else
                    endif
              endif
   40       continue

         if (mdtopo.gt.0) then
c           # check if we're in the dtopo region and need to refine:
c           # force refinement to level minleveldtopo 
            if (level.lt.minleveldtopo.and.t.le.tfdtopo.and.
     &              x2.gt.xlowdtopo.and.x1.lt.xhidtopo.and.
     &              y2.gt.ylowdtopo.and.y1.lt.yhidtopo) then
                amrflags(i,j)=DOFLAG
                go to 100 !# flagged, so no need to check anything else
                endif
             endif


         if (iqinit.gt.0 .and. t.eq.0.d0) then
c           # check if we're in the region where initial perturbation is
c           # specified and need to force refinement:
            if (level.lt.minlevelqinit.and.
     &              x2.gt.xlowqinit.and.x1.lt.xhiqinit.and.
     &              y2.gt.ylowqinit.and.y1.lt.yhiqinit) then
                amrflags(i,j)=DOFLAG
                go to 100 !# flagged, so no need to check anything else
                endif
             endif

c        -----------------------------------------------------------------

c        # refinement not forced, so check if it is allowed, and if so,
c        # check if there is a reason to flag this point:

         if (allowflag(x,y,t,level)) then

             depth= q(i,j,1)
             momentum = sqrt(q(i,j,2)**2 + q(i,j,3)**2)
             watersurface = q(i,j,1) + aux(i,j,1)
             if (depth.lt.drytolerance) watersurface=0.0 

c            # determine flowgrade
             do iflow=1,mflowgrades
                 if (iflowgradevariable(iflow).eq.1) then
                    flowgradenorm=depth
                    flowgradegrad=depth
                 elseif (iflowgradevariable(iflow).eq.2) then
                    flowgradenorm=momentum
                    flowgradegrad=momentum
                 elseif (iflowgradevariable(iflow).eq.3) then
                    flowgradenorm=dabs(watersurface)
                    flowgradegrad=momentum
                    endif

                 if (iflowgradetype(iflow).eq.1) then
                    flowgrademeasure=flowgradenorm
                 else
                    flowgrademeasure=flowgradegrad
                    endif

                 if (flowgrademeasure.gt.flowgradevalue(iflow)
     &                  .and.level.lt.iflowgrademinlevel(iflow)) then
                    amrflags(i,j)=DOFLAG
                    go to 100
                    endif

                 enddo
             endif

 100     continue  !# end loop on i
 200    continue   !# end loop on j

      return
      end
