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

      implicit double precision (a-h,o-z)
      external rpn2,rpt2

      include  "call.i"
      include  "fixedgrids.i"

      common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom

      parameter (msize=max1d+4)
      parameter (mwork=msize*(maxvar*maxvar + 13*maxvar + 3*maxaux +2))

      dimension q(nvar,mitot,mjtot)
      dimension fp(nvar,mitot,mjtot),gp(nvar,mitot,mjtot)
      dimension fm(nvar,mitot,mjtot),gm(nvar,mitot,mjtot)
      dimension aux(maux,mitot,mjtot)
      dimension work(mwork)

      logical    debug,  dump
      data       debug/.false./,  dump/.false./
c
c     # set tcom = time.  This is in the common block comxyt that could
c     # be included in the Riemann solver, for example, if t is explicitly
c     # needed there.

      tcom = time
      
      level = node(nestlevel,mptr)

      if (dump) then
         write(*,*)" dumping grid ",mptr
         do i = 1, mitot
         do j = 1, mjtot
            write(outunit,545) i,j,(q(ivar,i,j),ivar=1,nvar)
c            write(*,545) i,j,(q(ivar,i,j),ivar=1,nvar)
 545        format(2i4,4e15.7)
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
      mwork0 = (maxm+2*mbc)*(12*meqn + mwaves + meqn*mwaves + 2)
c
      if (mwork .lt. mwork0) then
         write(outunit,*) 'CLAW2 ERROR... mwork must be increased to ',
     &               mwork0
         write(*      ,*) 'CLAW2 ERROR... mwork must be increased to ',
     &               mwork0
         stop
      endif
c
c     # partition work array into pieces needed for local storage in
c     # step2 routine. Find starting index of each piece:
c
      i0faddm = 1
      i0faddp = i0faddm + (maxm+2*mbc)*meqn
      i0gaddm = i0faddp + (maxm+2*mbc)*meqn
      i0gaddp = i0gaddm + 2*(maxm+2*mbc)*meqn
      i0q1d   = i0gaddp + 2*(maxm+2*mbc)*meqn
      i0dtdx1 = i0q1d + (maxm+2*mbc)*meqn
      i0dtdy1 = i0dtdx1 + (maxm+2*mbc)
      i0aux1 = i0dtdy1 + (maxm+2*mbc)
      i0aux2 = i0aux1 + (maxm+2*mbc)*maux
      i0aux3 = i0aux2 + (maxm+2*mbc)*maux
c
c
      i0next = i0aux3 + (maxm+2*mbc)*maux    !# next free space
      mused  = i0next - 1                    !# space already used
      mwork1 = mwork - mused              !# remaining space (passed to step2)

c

c::::::::::::::::::::::::Fixed Grid Output:::::::::::::::::::::::::::::::::
      tc0=time !# start of computational step
      tcf=tc0+dt !# end of computational step

       call b4step2(mx,my,mbc,mx,my,nvar,q,
     &             xlowmbc,ylowmbc,dx,dy,time,dt,maux,aux)

c::::::::::::::::::::::::FIXED GRID DATA before step:::::::::::::::::::::::
c     # fill in values at fixed grid points effected at time tc0
      do ng=1,mfgrids

      if ((xlowfg(ng).lt.xlowmbc+mx*dx.and.xhifg(ng).gt.xlowmbc).and.
     &      (ylowfg(ng).lt.ylowmbc+my*dy.and.yhifg(ng).gt.ylowmbc).and.
     &      (ilastoutfg(ng).lt.noutfg(ng).and.tcf.ge.tstartfg(ng))) then


         if (tlastoutfg(ng)+dtfg(ng).ge.tc0.and.
     &                        tlastoutfg(ng)+dtfg(ng).le.tcf) then
c        # fixedgrid ng has an output time within [tc0,tcf] interval
c        # and it overlaps this computational grid spatially
         i0=i0fg(ng) !# index into the ng grid in the work array
         call fgridinterp(fgridearly(i0),xlowfg(ng),ylowfg(ng),
     &    xhifg(ng),yhifg(ng),dxfg(ng),dyfg(ng),mxfg(ng),myfg(ng),
     &    tc0,mfgridvars(ng),q,nvar,mx,my,mbc,dx,dy,nvar,xlowmbc,
     &    ylowmbc,maux,aux,ioutarrivaltimes(ng),ioutsurfacemax(ng),0)
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
         if (ioutsurfacemax(ng)+ioutarrivaltimes(ng).gt.0) then
           i0=i0fg2(ng)
         call fgridinterp(fgridoften(i0),xlowfg(ng),ylowfg(ng),
     &    xhifg(ng),yhifg(ng),dxfg(ng),dyfg(ng),mxfg(ng),myfg(ng),
     &    tc0,mfgridvars2(ng),q,nvar,mx,my,mbc,dx,dy,nvar,xlowmbc,
     &    ylowmbc,maux,aux,ioutarrivaltimes(ng),ioutsurfacemax(ng),2)
         endif
         endif

      endif
      enddo
      tcfmax=max(tcfmax,tcf)
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c
c     # take one step on the conservation law:
c
      call step2(mbig,mx,my,nvar,maux,
     &           mbc,mx,my,
     &              q,aux,dx,dy,dt,cflgrid,
     &              fm,fp,gm,gp,
     &              work(i0faddm),work(i0faddp),
     &              work(i0gaddm),work(i0gaddp),
     &              work(i0q1d),work(i0dtdx1),work(i0dtdy1),
     &              work(i0aux1),work(i0aux2),work(i0aux3),
     &              work(i0next),mwork1,rpn2,rpt2)
c
c            
        mptr_level = node(nestlevel,mptr)
        write(outunit,811) mptr, mptr_level, cflgrid
 811    format(" Courant # of grid ",i5," level",i3," is ",d12.4)
c
!$OMP  CRITICAL (cflm)
        cflmax = dmax1(cflmax,cflgrid)
        cfl_level = dmax1(cfl_level,cflgrid)
!$OMP END CRITICAL (cflm)
c
c       # update q
        dtdx = dt/dx
        dtdy = dt/dy
        do 50 m=1,nvar
        do 50 i=mbc+1,mitot-mbc
        do 50 j=mbc+1,mjtot-mbc
         if (mcapa.eq.0) then
c
c            # no capa array.  Standard flux differencing:

           q(m,i,j) = q(m,i,j)
     &           - dtdx * (fm(m,i+1,j) - fp(m,i,j))
     &           - dtdy * (gm(m,i,j+1) - gp(m,i,j))
         else
c            # with capa array.
           q(m,i,j) = q(m,i,j)
     &           - (dtdx * (fm(m,i+1,j) - fp(m,i,j))
     &           +  dtdy * (gm(m,i,j+1) - gp(m,i,j))) / aux(mcapa,i,j)
         endif

 50      continue
c
c     # Copied here from b4step2 since need to do before saving to qc1d:
      do i=1,mitot
        do j=1,mjtot
          if (q(1,i,j).lt.drytolerance) then
             q(1,i,j) = max(q(1,i,j),0.d0)
             do m=2,nvar
                q(m,i,j)=0.d0
                enddo
             endif
        enddo
      enddo
c
      if (method(5).eq.1) then
c        # with source term:   use Godunov splitting
         call src2(mx,my,nvar,mbc,mx,my,xlowmbc,ylowmbc,dx,dy,
     &             q,maux,aux,time,dt)
         endif


c     ::::::::::::::::::::::::Fixed Grid data afterstep:::::::::::::::::::::::
c     # fill in values at fixed grid points effected at time tcf
      do ng=1,mfgrids
      if ((xlowfg(ng).lt.xlowmbc+mx*dx.and.xhifg(ng).gt.xlowmbc).and.
     &      (ylowfg(ng).lt.ylowmbc+my*dy.and.yhifg(ng).gt.ylowmbc).and.
     &      (ilastoutfg(ng).lt.noutfg(ng).and.tcf.ge.tstartfg(ng))) then

         if (tlastoutfg(ng)+dtfg(ng).ge.tc0
     &                     .and.tlastoutfg(ng)+dtfg(ng).le.tcf) then


c        # fixedgrid ng has an output time within [tc0,tcf] interval
c        # and it overlaps this computational grid spatially
         i0=i0fg(ng) !# index into the ng grid in the work array

         call fgridinterp(fgridlate(i0),xlowfg(ng),ylowfg(ng),
     &    xhifg(ng),yhifg(ng),dxfg(ng),dyfg(ng),mxfg(ng),myfg(ng),
     &    tcf,mfgridvars(ng),q,nvar,mx,my,mbc,dx,dy,nvar,xlowmbc,
     &    ylowmbc,maux,aux,ioutarrivaltimes(ng),ioutsurfacemax(ng),0)
c            # routine to interpolate solution
c            # at tcf to the fixed grid storage array,
c            #saving solution and tcf at every grid point

         endif

c        # fill in values for eta if they need to be saved for later checking max/mins
c        # check for arrival times
         if (ioutsurfacemax(ng)+ioutarrivaltimes(ng).gt.0) then
         i0=i0fg2(ng)
         call fgridinterp(fgridoften(i0),xlowfg(ng),ylowfg(ng),
     &    xhifg(ng),yhifg(ng),dxfg(ng),dyfg(ng),mxfg(ng),myfg(ng),
     &    tc0,mfgridvars2(ng),q,nvar,mx,my,mbc,dx,dy,nvar,xlowmbc,
     &    ylowmbc,maux,aux,ioutarrivaltimes(ng),ioutsurfacemax(ng),1)
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
         if (level .eq. levelcheck) then
         if (ioutsurfacemax(ng)+ioutarrivaltimes(ng).gt.0) then
           i0=i0fg2(ng)
         call fgridinterp(fgridoften(i0),xlowfg(ng),ylowfg(ng),
     &    xhifg(ng),yhifg(ng),dxfg(ng),dyfg(ng),mxfg(ng),myfg(ng),
     &    tc0,mfgridvars2(ng),q,nvar,mx,my,mbc,dx,dy,nvar,xlowmbc,
     &    ylowmbc,maux,aux,ioutarrivaltimes(ng),ioutsurfacemax(ng),2)
          endif
          endif


      endif
      enddo
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
         write(*,*)" at end of stepgrid: dumping grid ",mptr
         do i = 1, mitot
         do j = 1, mjtot
            write(outunit,545) i,j,(q(ivar,i,j),ivar=1,nvar)
c            write(*,545) i,j,(q(i,j,ivar),ivar=1,nvar)
         end do
         end do
      endif
c
      return
      end


