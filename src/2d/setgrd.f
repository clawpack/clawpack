c
c -----------------------------------------------------------------
c
      subroutine setgrd (nvar,cut,naux,dtinit,t0)
c
      use geoclaw_module
c
      implicit double precision (a-h,o-z)

      include  "call.i"
      dimension spoh(maxlv)

      integer verbosity_regrid
      logical  vtime
      data     vtime/.false./

c     # may as well not bother to calculate time step for error est.
c
c :::::::::::::::::::::::::::: SETGRD :::::::::::::::::::::::::::::::;
c  set up the entire tree/grid structure.  only at this time t = 0
c  can we take advantage of initialization routines.
c  remember that regridding/error estimation needs to have two
c  time steps of soln. values.
c  6/21/05: added dtinit arg. to allow for better choice of initial timestep
c   as discovered by advance/setgrd in first step.
c ::::::::::::::::::::::::::::::::::::::;::::::::::::::::::::::::::::;
c
      dtinit = possk(1)
      if (mxnest .eq. 1) go to 99
c
      levnew =  2
      time   = t0
      verbosity_regrid = method(4)
c
 10   if (levnew .gt. mxnest) go to 30
          levold = levnew - 1
          if (lstart(levold) .eq. 0) go to 30
          lbase  = levold
          lfnew  = lbase
c
c  set up level to be flagged. need a solution t=0,and t=dt.
c  error estimation makes next one at t=2*dt for Richardson.
c
         call advanc(levold,nvar,dtlev,vtime,naux)
         evol = evol + rvol
         rvol = 0.d0
         kfac = 1
         do i = 1, levold-1
           kfac = kfac * kratio(i)
         end do
         dtinit = min(dtinit, dtlev*kfac)
 
c        dont count it in real integration stats
         do 20 level=1,mxnest
 20         rvoll(level) = 0.d0
c
c  flag, cluster, and make new grids
c
         call grdfit(lbase,levold,nvar,naux,cut,time,t0)
         if (newstl(levnew) .ne. 0) lfnew = levnew
c
c  init new level. after each iteration. fix the data structure
c  also reinitalize coarser grids so fine grids can be advanced
c  and interpolate correctly for their bndry vals from coarser grids.
c
         call ginit(newstl(levnew),.true., nvar, naux, t0)
         lstart(levnew) = newstl(levnew)
         lfine = lfnew
         call ginit(lstart(levold),.false., nvar, naux, t0)
c
c count number of grids on newly created levels (needed for openmp
c parallelization). this is also  done in regridding.
c  set up numgrids now for new level, since advanc will use it for parallel execution
c 
         mptr = lstart(levnew)
         ngrids = 0
         ncells = 0
         do while (mptr .gt. 0)
            ngrids = ngrids + 1
            ncells = ncells + (node(ndihi,mptr)-node(ndilo,mptr)+1)
     .                      * (node(ndjhi,mptr)-node(ndjlo,mptr)+1)
            mptr = node(levelptr, mptr)
          end do
          numgrids(levnew) = ngrids
          numcells(levnew) = ncells
          if (verbosity_regrid .ge. levnew) then
             write(*,100) ngrids,ncells,levnew
 100         format("there are ",i4," grids with ",i8,
     &               " cells at level ", i3)
          endif 
c
      levnew = levnew + 1
      go to 10
 30   continue
c
c  switch location of old and new storage for soln. vals, 
c  and reset time to 0.0 (or initial time t0)
c
      if (mxnest .eq. 1) go to 99
c
      lev = 1
 40   if ((lev .eq. mxnest) .or. (lev .gt. lfine))  go to 60
        mptr = lstart(lev)
 50        itemp                = node(store1,mptr)
           node(store1,mptr)    = node(store2,mptr)
           node(store2,mptr)    = itemp
           rnode(timemult,mptr) = t0
           mptr                 = node(levelptr,mptr)
           if (mptr .ne. 0) go to 50
       lev = lev + 1
       go to 40
 60   continue
c
c initial updating so can do conservation check. can do before
c bndry flux arrays set, since dont need them for this
c
      do 65 level = 1, lfine-1
         call update(lfine-level,nvar,naux)
 65   continue
c
c set up boundary flux conservation arrays
c
      do 70 level = 1, lfine-1
         call prepf(level+1,nvar,naux)
         call prepc(level,nvar)
 70   continue
c
      if (.not. varRefTime) go to 99   ! keep consistent with gfixup and filval

         do level = 1, lfine   ! compute max speed for all grids at each level
           spoh(level) = 0.d0 ! to keep track of max wave speed for all new grids
           mptr = lstart(level)
 75        continue
             loc    = node(store1,mptr)
             locaux = node(storeaux,mptr)
             nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
             ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
             mitot  = nx + 2*nghost
             mjtot  = ny + 2*nghost

             sp_over_h = get_max_speed(alloc(loc),mitot,mjtot,nvar,
     &                         alloc(locaux),naux,nghost,
     &                         hxposs(level),hyposs(level))
c            write(*,*)"max wave speed for level ",level," grid, ",mptr,
c    .        " is", sp_over_h
             spoh(level) = max(spoh(level),sp_over_h)
           mptr = node(levelptr,mptr)
           if (mptr .ne. 0) go to 75
         end do
c
c check that level 1 has a time step that doesnt exceed cfl based on wave speed just calculated
c
         dtc = possk(1)
         if (dtc*spoh(1) .gt. cflv1) then
              write(*,*)" coarse level time step too big: should reset"
              write(*,*) " dtc ", dtc," gives cfl = ",dtc*spoh(1)
         endif
c
c        even if initial timestep was good could be too small
c        try automatically resetting 
c        RANDY: do you like this option
c
         possk(1) = cfl/spoh(1)
         write(*,*)"  automatically setting initial dt to ",possk(1)

c
c  set new time step and refinement ratios in time
c  at this initial time just setting initial conditions
c  important to do this from coarser to finer levels to make sure 
c  subcycling lines up right.
c  this sets ratios, and possk array.  
c  note that code below does not change coarsest level time step,
c  this is set during timestepping or by user initially
c
      do level = 2, lfine   
         dtc = possk(level-1)
         dtf = cfl/spoh(level)
         if (dtf .gt. dtc) then
            kratio(level-1) = 1  ! cant have larger timestep than parent level
            possk(level)    = dtc  ! cant have larger timestep than parent level
         else
            kratio(level-1) = ceiling(dtc/dtf)
            possk(level)    = possk(level-1)/kratio(level-1)
        endif
c       write(6,*)" setting ref. ratio in time for level ",level," to ",
c    .            kratio(level-1)
      end do
c
 99   continue
c
      return
      end
