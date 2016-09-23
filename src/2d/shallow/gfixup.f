c
c  -----------------------------------------------------------
c
      subroutine gfixup(lbase, lfnew, nvar, naux, newnumgrids,
     .                  maxnumnewgrids)
c
      use geoclaw_module
      use refinement_module, only: varRefTime
      use amr_module
      use topo_module, only: topo_finalized,aux_finalized
      implicit double precision (a-h,o-z)

      dimension spoh(maxlv)
      integer omp_get_thread_num, omp_get_max_threads
      integer mythread/0/, maxthreads/1/
      integer newnumgrids(maxlv),listnewgrids(maxnumnewgrids)

      integer clock_start, clock_finish, clock_rate
      integer wallclock_start, wallclock_finish
      integer mbad
c
c ::::::::::::::::::::::::: GFIXUP ::::::::::::::::::::::::::::::::;
c        interpolate initial values for the newly created grids.
c        the start of each level is located in newstl array.
c        since only levels greater than lbase were examined, start
c        looking there.
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;
c
c   reclaim old storage (position 8) and list space 15 and 16
c   before allocating new storage. remember, finest level grids
c  (if level = mxnest so that error never estimated) don't have
c  2 copies of solution values at old and new times.
c
c
      call putsp(lbase,lbase,nvar,naux)
      level = lbase + 1
 1    if (level .gt. lfine) go to 4
      call putsp(lbase,level,nvar,naux)
          mptr = lstart(level)
 2        if (mptr .eq. 0) go to 3
              nx = node(ndihi,mptr) - node(ndilo,mptr) + 1
              ny = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
              mitot = nx + 2*nghost
              mjtot = ny + 2*nghost
              nwords        = mitot*mjtot*nvar
              if (level .lt. mxnest) 
     .           call reclam(node(store2, mptr), nwords)
              node(store2, mptr) = 0
              mptr          = node(levelptr, mptr)
          go to 2
 3        level   = level + 1
          go to 1
c
 4    lcheck = lbase + 1

c
c update topo to current time, before setting all new grids at new levels
c
      time = rnode(timemult, lstart(lbase))
      if (.not. topo_finalized) then
          call topo_update(time)
      endif

 5    if (lcheck .gt. mxnest) go to 89
          hx = hxposs(lcheck)
          hy = hyposs(lcheck)
          spoh(lcheck) = 0.d0 ! to keep track of max wave speed for all new grids
c
c prepare for doing next loop over grids at a given level in parallel
c unlike other level loops, these are newly created grids, not yet merged in
c so take grids from newstl (NEWSTartOfLevel), not lstart. Dont yet know
c how many either.
       call prepnewgrids(listnewgrids,newnumgrids(lcheck),lcheck)
c
c  interpolate level lcheck
c   first get space, since cant do that part in parallel
       do  j = 1, newnumgrids(lcheck)
          mptr = listnewgrids(j)
            nx = node(ndihi,mptr) - node(ndilo,mptr) + 1
            ny = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
            mitot = nx + 2*nghost
            mjtot = ny + 2*nghost
            loc    = igetsp(mitot * mjtot * nvar)
            node(store1, mptr)  = loc
            if (naux .gt. 0) then
              locaux = igetsp(mitot * mjtot * naux)
             else
              locaux = 1
            endif
            node(storeaux, mptr)  = locaux
       end do

c   
c                 other reduction variables initialized in stst1
       this_spoh = 0.d0
       call system_clock(wallclock_start,clock_rate)  

!$OMP PARALLEL DO 
!$OMP&            PRIVATE(clock_start,clock_finish,clock_rate)
!$OMP&            PRIVATE(j,mptr,nx,ny,mitot,mjtot,corn1,corn2,loc)
!$OMP&            PRIVATE(locaux,time,mic,mjc,xl,xr,yb,yt,ilo,ihi)
!$OMP&            PRIVATE(jlo,jhi,sp_over_h)
!$OMP&            SHARED(newnumgrids,listnewgrids,nghost,node,hx,hy)
!$OMP&            SHARED(rnode,intratx,intraty,lcheck,nvar,alloc,naux)
!$OMP&            REDUCTION(MAX:this_spoh)
!$OMP&            REDUCTION(+:timeFilval)
!$OMP&            SCHEDULE(dynamic,1)
!$OMP&            DEFAULT(none)

      do  j = 1, newnumgrids(lcheck)
          mptr = listnewgrids(j)

c  changed to move setaux out of this loop. instead, copy aux in filval 
c  along with soln.involves changing intcopy to icall and making flag array
c  can only do this after topo stops moving
              nx = node(ndihi,mptr) - node(ndilo,mptr) + 1
              ny = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
              mitot = nx + 2*nghost
              mjtot = ny + 2*nghost
              corn1 = rnode(cornxlo,mptr)
              corn2 = rnode(cornylo,mptr)
              loc   =  node(store1, mptr)
              if (naux .gt. 0) then
                locaux =  node(storeaux, mptr)
              else
                locaux = 1
              endif

c
c      We now fill in the values for grid mptr using filval. It uses
c      piecewise linear interpolation to obtain values from the
c      (lcheck - 1) grid, then overwrites those with whatever (lcheck)
c      grids are available. We take advantage of the fact that the
c      (lcheck - 1) grids have already been set, and that the distance
c      between any point in mptr and a (lcheck - 1) and (lcheck - 2)
c      interface is at least one (lcheck - 1) cell wide.
c
 
c          # make a coarsened patch with ghost cells so can use
c          # grid interpolation routines, but only set "interior".
c          # extra 2 cells so that can use linear interp. on
c          # "interior" of coarser patch to fill fine grid.
           mic = nx/intratx(lcheck-1) + 2
           mjc = ny/intraty(lcheck-1) + 2
           xl  = rnode(cornxlo,mptr)
           xr  = rnode(cornxhi,mptr)
           yb  = rnode(cornylo,mptr)
           yt  = rnode(cornyhi,mptr)
           time  = rnode(timemult, mptr)
           ilo   = node(ndilo, mptr)
           ihi   = node(ndihi, mptr)
           jlo   = node(ndjlo, mptr)
           jhi   = node(ndjhi, mptr)
 
           call system_clock(clock_start,clock_rate)
           call filval(alloc(loc),mitot,mjtot,hx,hy,lcheck,time,
     1                 mic,mjc,
     2                 xl,xr,yb,yt,nvar,
     3                 mptr,ilo,ihi,jlo,jhi,
     4                 alloc(locaux),naux,
     5                 sp_over_h)
           call system_clock(clock_finish,clock_rate)
           timeFilval = timeFilval + clock_finish - clock_start
           this_spoh = max(this_spoh, sp_over_h)
         end do
!$OMP END PARALLEL DO

       call system_clock(wallclock_finish,clock_rate)  
       timeFilvalTot = timeFilvalTot + wallclock_finish-wallclock_start

       spoh(lcheck) = this_spoh
 
 
c
c  done filling new grids at level. move them into lstart from newstl
c  (so can use as source grids for filling next level). can also
c  get rid of loc. 7 storage for old level.
c
 80   mptr = lstart(lcheck)
 85   if (mptr .eq. 0) go to 90
          nx = node(ndihi,mptr) - node(ndilo,mptr) + 1
          ny = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
          mitot = nx + 2*nghost
          mjtot = ny + 2*nghost
          call reclam(node(store1,mptr),mitot*mjtot*nvar)
          if (naux .gt. 0) then
            call reclam(node(storeaux,mptr),mitot*mjtot*naux)
          endif
          mold   = mptr
          mptr   = node(levelptr,mptr)
          call putnod(mold)
          call freeBndryList(mold)
          go to 85
 90   lstart(lcheck) = newstl(lcheck)
      lcheck = lcheck + 1
      go to 5
c
 89   lfine = lfnew
c
c     initialize 2nd (old time) storage block for new grids not at top level
c
      levend = min(lfine,mxnest-1)
      do 110 level = lbase+1, levend
         mptr = lstart(level)
 105     if (mptr .eq. 0) go to 110
            nx = node(ndihi,mptr) - node(ndilo,mptr) + 1
            ny = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
            mitot = nx + 2*nghost
            mjtot = ny + 2*nghost
            nwords = mitot*mjtot*nvar
            node(store2,mptr) = igetsp(nwords)
         mptr = node(levelptr,mptr)
         go to 105
 110   continue

c
c -------------
c  grid structure now complete again. safe to print, etc. assuming
c  things initialized to zero in nodget.
c -------------
c
       if (.not. varRefTime) go to 99
c
c  set new time step and refinement ratios in time
c  important to do this from coarser to finer levels to make sure 
c  subcycling lines up right.
c  this sets ratios, and possk array.  
c  note that code below does not change coarsest level time step,
c  this is set during timestepping
c
      do level = lbase+1, lfine   
         dtc = possk(level-1)
C        Need this check in order to account for completely dry grids
         if (abs(spoh(level)) < tiny(1.d0)) then
           dtf = huge(1.d0)
         else
           dtf = cfl/(spoh(level))
         endif
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

 99   return
      end
c
c -----------------------------------------------------------------------------------------
c
c  use different routine since need to scan new grid list (newstl) not lstart
c  to make grids.  
c  could make one routine by passing in source of list, but this changed 4 other routines
c  so I didnt want to have to deal with it

       subroutine prepnewgrids(listnewgrids,num,level)

       use amr_module
       implicit double precision (a-h,o-z)
       integer listnewgrids(num)

       mptr = newstl(level)
       do j = 1, num
          listnewgrids(j) = mptr
          mptr = node(levelptr, mptr)
       end do

       if (mptr .ne. 0) then
         write(*,*)" Error in routine setting up grid array "
         stop
       endif

       return
       end
