c
c -------------------------------------------------------------
c
      subroutine stepgrid(q,fm,fp,gm,gp,mitot,mjtot,mbc,dt,dtnew,dx,dy,
     &                  nvar,xlow,ylow,time,mptr,maux,aux)
c
c
c ::::::::::::::::::: STEPGRID ::::::::::::::::::::::::::::::::::::
c take a time step on a single grid. overwrite solution array q.
c A modified version of the clawpack routine step2 is used.
c
c return fluxes in fm,fp and gm,gp.
c patch has room for ghost cells (mbc of them) around the grid.
c everything is the enlarged size (mitot by mjtot).
c
c mbc       = number of ghost cells  (= lwidth)
c mptr      = grid number  (for debugging)
c xlow,ylow = lower left corner of enlarged grid (including ghost cells).
c dt         = incoming time step
c dx,dy      = mesh widths for this grid
c dtnew      = return suggested new time step for this grid's soln.
c
c
c
c      This version of stepgrid, stepgrid_geo.f allows output on
c      fixed grids specified in setfixedgrids.data
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      use geoclaw_module
      use amr_module
      use fixedgrids_module
      implicit double precision (a-h,o-z)

      external rpn2,rpt2


c      parameter (msize=max1d+4)
c      parameter (mwork=msize*(maxvar*maxvar + 13*maxvar + 3*maxaux +2))

      dimension q(nvar,mitot,mjtot)
      dimension fp(nvar,mitot,mjtot),gp(nvar,mitot,mjtot)
      dimension fm(nvar,mitot,mjtot),gm(nvar,mitot,mjtot)
      dimension aux(maux,mitot,mjtot)
c      dimension work(mwork)

      logical :: debug = .false.
      logical :: dump = .false.
c
      tcfmax = -rinfinity
      level = node(nestlevel,mptr)

      if (dump) then
         write(outunit,*)" at start of stepgrid: dumping grid ",mptr
         do i = 1, mitot
         do j = 1, mjtot
            write(outunit,545) i,j,(q(ivar,i,j),ivar=1,nvar),
     .                         (aux(iaux,i,j),iaux=1,maux)
 545        format(2i4,4e15.7,/,8x,4e15.7)
         end do
         end do
      endif
c
      meqn   = nvar
      mx = mitot - 2*mbc
      my = mjtot - 2*mbc
      maxm = max(mx,my)       !# size for 1d scratch array
      mbig = maxm
      xlowmbc = xlow + mbc*dx
      ylowmbc = ylow + mbc*dy

c     # method(2:7) and mthlim
c     #    are set in the amr2ez file (read by amr)
c
      method(1) = 0
c
c
c     # fluxes initialized in step2
c
C       mwork0 = (maxm+2*mbc)*(12*meqn + mwaves + meqn*mwaves + 2)
C c
C       if (mwork .lt. mwork0) then
C          write(outunit,*) 'CLAW2 ERROR... mwork must be increased to ',
C      &               mwork0
C          write(*      ,*) 'CLAW2 ERROR... mwork must be increased to ',
C      &               mwork0
C          stop
C       endif
c
c     # partition work array into pieces needed for local storage in
c     # step2 routine. Find starting index of each piece:
c
C       i0faddm = 1
C       i0faddp = i0faddm + (maxm+2*mbc)*meqn
C       i0gaddm = i0faddp + (maxm+2*mbc)*meqn
C       i0gaddp = i0gaddm + 2*(maxm+2*mbc)*meqn
C       i0q1d   = i0gaddp + 2*(maxm+2*mbc)*meqn
C       i0dtdx1 = i0q1d + (maxm+2*mbc)*meqn
C       i0dtdy1 = i0dtdx1 + (maxm+2*mbc)
C       i0aux1 = i0dtdy1 + (maxm+2*mbc)
C       i0aux2 = i0aux1 + (maxm+2*mbc)*maux
C       i0aux3 = i0aux2 + (maxm+2*mbc)*maux
C c
C c
C       i0next = i0aux3 + (maxm+2*mbc)*maux    !# next free space
C       mused  = i0next - 1                    !# space already used
C       mwork1 = mwork - mused              !# remaining space (passed to step2)

c

c::::::::::::::::::::::::Fixed Grid Output:::::::::::::::::::::::::::::::::
      tc0=time !# start of computational step
      tcf=tc0+dt !# end of computational step

!$OMP CRITICAL (FixedGrids)
c     # see if any f-grids should be written out
      do ng=1,num_fixed_grids
        if (tc0 > fgrids(ng)%start_time .and. 
     &      fgrids(ng)%last_output_index < fgrids(ng)%num_output) then
c     # fgrid ng may need to be written out
c     # find the first output number that has not been written out and
c     # find the first output number on a fixed grid that is >= tc0
c     # which will not be written out
           if (fgrids(ng)%dt > 0.d0) then
             ioutfgend= 1+max(0,nint((tc0 - fgrids(ng)%start_time)
     &                                 / fgrids(ng)%dt))
           else
             ioutfgend=1
           endif
           ioutfgend = min(ioutfgend,fgrids(ng)%num_output)
           ioutfgstart = fgrids(ng)%last_output_index + 1
c     # write-out fgrid times that are less than tc0, and have not been written yet
c     # these should be the most accurate values at any given point in the fgrid
c     # since tc0> output time
           do ioutfg=ioutfgstart,ioutfgend
             toutfg=fgrids(ng)%start_time+(ioutfg-1)*fgrids(ng)%dt
             if (toutfg < tc0) then
c               # write out the solution for fixed grid ng
c               # test if arrival times should be output
                ioutflag = fgrids(ng)%output_arrival_times*
     &                         (fgrids(ng)%num_output-
     &                          fgrids(ng)%last_output_index)
                call fgrid_out(ng,fgrids(ng),toutfg,ioutfg,ioutflag)

                fgrids(ng)%last_output_time = toutfg
                fgrids(ng)%last_output_index = 
     &                               fgrids(ng)%last_output_index + 1
             endif
           enddo

        endif
      enddo
!$OMP END CRITICAL (FixedGrids)
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

       call b4step2(mbc,mx,my,nvar,q,
     &             xlowmbc,ylowmbc,dx,dy,time,dt,maux,aux)
      
c::::::::::::::::::::::::FIXED GRID DATA before step:::::::::::::::::::::::
c     # fill in values at fixed grid points effected at time tc0
!$OMP CRITICAL (FixedGrids)
      do ng=1,num_fixed_grids

      if ( (fgrids(ng)%x_low < xlowmbc + mx*dx) .and.
     &     (fgrids(ng)%x_hi  > xlowmbc) .and.
     &     (fgrids(ng)%y_low < ylowmbc + my*dy) .and.
     &     (fgrids(ng)%y_hi  > ylowmbc) .and.
     &     (fgrids(ng)%last_output_index < fgrids(ng)%num_output) .and.
     &     (tcf >= fgrids(ng)%start_time) ) then
         
         if (fgrids(ng)%last_output_time + fgrids(ng)%dt >= tc0 .and.
     &       fgrids(ng)%last_output_time + fgrids(ng)%dt <= tcf) then

c        # fixedgrid ng has an output time within [tc0,tcf] interval
c        # and it overlaps this computational grid spatially
         call fgrid_interp(1,fgrids(ng),tc0,q,nvar,mx,my,mbc,dx,dy,
     &                     xlowmbc,ylowmbc,maux,aux,0)
     
c         # routine to spatially interpolate computational solution
c         # at tc0 to the fixed grid spatial points,
c         #saving solution, variables and tc0 at every grid point
         endif

c        # set maxima or minima if this is a new coarse step
c        if (tc0.ge.tcfmax) then

c        # RJL: rewrote to set min/max every time a grid at level 1
c        # is about to be taken.  The previous code failed if there was more than one grid
c        # at level 1.   Note that all grids are up to date at start of step on level 1.
c        # New feature added at end of this routine to check more frequently if
c        # levelcheck > 0.
         if (level .eq. 1) then
         if (fgrids(ng)%output_surface_max 
     &       + fgrids(ng)%output_arrival_times > 0) then
     
         call fgrid_interp(3,fgrids(ng),tc0,q,nvar,mx,my,mbc,dx,dy,
     &                     xlowmbc,ylowmbc,maux,aux,2)
     
         endif
         endif

      endif
      enddo
      tcfmax=max(tcfmax,tcf)

!$OMP END CRITICAL (FixedGrids)

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c     New fixed grid stuff: Update fixed grid info from this patch...

      call fgmax_frompatch(mx,my,nvar,mbc,maux,q,aux,
     &     dx,dy,xlowmbc,ylowmbc,level,time,time+dt)
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c     # take one step on the conservation law:
c
      call step2(mbig,nvar,maux,
     &           mbc,mx,my,
     &              q,aux,dx,dy,dt,cflgrid,
     &              fm,fp,gm,gp,rpn2,rpt2)
c
c            
        mptr_level = node(nestlevel,mptr)

c       write(outunit,811) mptr, mptr_level, cflgrid
c811    format(" Courant # of grid ",i5," level",i3," is ",d12.4)
c

!$OMP  CRITICAL (cflm)

        cflmax = dmax1(cflmax,cflgrid)
        cfl_level = dmax1(cfl_level,cflgrid)

!$OMP END CRITICAL (cflm)

c
!        write(outunit,*)" updating grid ", mptr
c       # update q
        dtdx = dt/dx
        dtdy = dt/dy
         if (mcapa.eq.0) then
          do 50 j=mbc+1,mjtot-mbc
          do 50 i=mbc+1,mitot-mbc
          do 50 m=1,nvar
c
c            # no capa array.  Standard flux differencing:

           q(m,i,j) = q(m,i,j)
     &           - dtdx * (fm(m,i+1,j) - fp(m,i,j))
     &           - dtdy * (gm(m,i,j+1) - gp(m,i,j))
 50       continue
         else
          do 51 j=mbc+1,mjtot-mbc
          do 51 i=mbc+1,mitot-mbc
          do 51 m=1,nvar
c            # with capa array.
           q(m,i,j) = q(m,i,j)
     &           - (dtdx * (fm(m,i+1,j) - fp(m,i,j))
     &           +  dtdy * (gm(m,i,j+1) - gp(m,i,j))) / aux(mcapa,i,j)
 51       continue
!           write(outunit,543) m,i,j,q(m,i,j),fm(m,i+1,j),fp(m,i,j),
!     .        gm(m,i,j+1), gp(m,i,j)
543       format(3i4,5e25.16)

         endif

c 50      continue
c
c     # Copied here from b4step2 since need to do before saving to qc1d:
      forall(i=1:mitot, j=1:mjtot, q(1,i,j) < dry_tolerance)
        q(1,i,j) = max(q(1,i,j),0.d0)
        q(2:meqn,i,j) = 0.d0
      end forall
c
      if (method(5).eq.1) then
c        # with source term:   use Godunov splitting
         call src2(nvar,mbc,mx,my,xlowmbc,ylowmbc,dx,dy,
     &             q,maux,aux,time,dt)
         endif

!$OMP CRITICAL (FixedGrids)
c     ::::::::::::::::::::::::Fixed Grid data afterstep:::::::::::::::::::::::
c     # fill in values at fixed grid points effected at time tcf
      do ng=1,num_fixed_grids
      if ((fgrids(ng)%x_low < xlowmbc + mx * dx) .and.
     &    (fgrids(ng)%x_hi  > xlowmbc) .and.
     &    (fgrids(ng)%y_low < ylowmbc + my * dy) .and.
     &    (fgrids(ng)%y_hi  > ylowmbc) .and.
     &    (fgrids(ng)%last_output_index < fgrids(ng)%num_output) .and.
     &    (tcf >= fgrids(ng)%start_time)) then
      
        if (fgrids(ng)%last_output_time + fgrids(ng)%dt >= tc0 .and.
     &      fgrids(ng)%last_output_time + fgrids(ng)%dt <= tcf) then

c        # fixedgrid ng has an output time within [tc0,tcf] interval
c        # and it overlaps this computational grid spatially
C         i0=i0fg(ng) !# index into the ng grid in the work array

        call fgrid_interp(2,fgrids(ng),tcf,q,nvar,mx,my,mbc,dx,dy,
     &                    xlowmbc,ylowmbc,maux,aux,0)

c            # routine to interpolate solution
c            # at tcf to the fixed grid storage array,
c            #saving solution and tcf at every grid point

        endif

c        # fill in values for eta if they need to be saved for later checking max/mins
c        # check for arrival times
        if (fgrids(ng)%output_surface_max 
     &      + fgrids(ng)%output_arrival_times > 0) then

        call fgrid_interp(3,fgrids(ng),tc0,q,nvar,mx,my,mbc,dx,dy, 
     &                    xlowmbc,ylowmbc,maux,aux,1)
     
        endif
         
c        # RJL: Modified 8/20/11 
c        # If levelcheck > 0 then update max/mins at end of step on this grid.
c        # Note that if there are finer grids then fgridoften will not have been updated
c        # properly yet by those grids.  This modification allows checking max/min more
c        # frequently than the original code (equivalent to levelcheck==0) when you know
c        # what level is most relevant for this fixed grid.  Note also that if there are no
c        # grids at levelcheck overlapping a portion of the fixed grid then the max/min values 
c        # will be updated only at start of next level 1 step.
 
        levelcheck = 0 
        if (level == levelcheck) then
        if (fgrids(ng)%output_arrival_times 
     &      + fgrids(ng)%output_surface_max > 0) then

        call fgrid_interp(3,fgrids(ng),tc0,q,nvar,mx,my,mbc,dx,dy,
     &    xlowmbc,ylowmbc,maux,aux,2)
         endif
         endif


      endif
      enddo
!$OMP END CRITICAL (FixedGrids)
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


c     # output fluxes for debugging purposes:
      if (debug) then
         write(dbugunit,*)" fluxes for grid ",mptr
         do 830 i = mbc+1, mitot-1
            do 830 j = mbc+1, mjtot-1
               write(dbugunit,831) i,j,fm(1,i,j),fp(1,i,j),
     .                                 gm(1,i,j),gp(1,i,j)
               do 830 m = 2, meqn
                  write(dbugunit,832) fm(m,i,j),fp(m,i,j),
     .            gm(m,i,j),gp(m,i,j)
  831          format(2i4,4d16.6)
  832          format(8x,4d16.6)
  830    continue
      endif

c
c
c For variable time stepping, use max speed seen on this grid to
c choose the allowable new time step dtnew.  This will later be
c compared to values seen on other grids.
c
       if (cflgrid .gt. 0.d0) then
           dtnew = dt*cfl/cflgrid
         else
c          # velocities are all zero on this grid so there's no
c          # time step restriction coming from this grid.
            dtnew = rinfinity
          endif

c     # give a warning if Courant number too large...
c
      if (cflgrid .gt. cflv1) then
            write(*,810) cflgrid,cflv1,mptr,mptr_level
            write(outunit,810) cflgrid, cflv1,mptr,mptr_level
  810       format('*** WARNING *** Courant number  =', d12.4,
     &              '  is larger than input cfl_max = ', d12.4,
     &              '  on grid ',i3, ' level ',i3)
            endif
c
      if (dump) then
         write(outunit,*)" at end of stepgrid: dumping grid ",mptr
         do i = mbc+1, mitot-mbc
         do j = mbc+1, mjtot-mbc
            write(outunit,545) i,j,(q(ivar,i,j),ivar=1,nvar)
c            write(*,545) i,j,(q(i,j,ivar),ivar=1,nvar)
         end do
         end do
      endif
c
      return
      end


