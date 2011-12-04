c
c --------------------------------------------------------------
c
      subroutine advanc (level,nvar,dtlevnew,vtime,naux)
c
      implicit double precision (a-h,o-z)

      include  "call.i"

      logical vtime
      integer omp_get_thread_num, omp_get_max_threads
      integer mythread/0/, maxthreads/1/
      integer listgrids(numgrids(level))

c     maxgr is maximum number of grids  many things are
c     dimensioned at, so this is overall. only 1d array
c     though so should suffice. problem is
c     not being able to dimension at maxthreads

      integer locfp_save(maxgr)
      integer locgp_save(maxgr)


c
c  ::::::::::::::; ADVANC :::::::::::::::::::::::::::::::::::::::::::
c  integrate all grids at the input  'level' by one step of its delta(t)
c  this includes:  setting the ghost cells 
c                  advancing the solution on the grid
c                  adjusting fluxes for flux conservation step later
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      hx   = hxposs(level)
      hy   = hyposs(level)
      delt = possk(level)
c     this is linear alg.
      call prepgrids(listgrids,numgrids(level),level)
c     maxthreads initialized to 1 above in case no openmp
!$    maxthreads = omp_get_max_threads()

!$OMP PARALLEL DO PRIVATE(j,locnew, locaux, mptr,nx,ny,mitot
!$OMP&                    ,mjtot,time),
!$OMP&            SHARED(level,nvar,naux,alloc,intrat,delt,
!$OMP&                   nghost,node,rnode,numgrids,listgrids),
!$OMP&            SCHEDULE (dynamic,1)
ccccccc!$OMP&            DEFAULT(none),
      do  j = 1, numgrids(level)
c          mget is n^2 alg.
c          mptr   = mget(j,level)
          mptr = listgrids(j)
          nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
          ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
          mitot  = nx + 2*nghost
          mjtot  = ny + 2*nghost
          locnew = node(store1,mptr)
          locaux = node(storeaux,mptr)
          time   = rnode(timemult,mptr)
c
          call bound(time,nvar,nghost,alloc(locnew),mitot,mjtot,mptr,
     1               alloc(locaux),naux)

        end do
!$OMP END PARALLEL DO

c
c save coarse level values if there is a finer level for wave fixup
      if (level+1 .le. mxnest) then
         if (lstart(level+1) .ne. null) then
            call saveqc(level+1,nvar,naux)
         endif
      endif
c
      time = rnode(timemult,lstart(level))
      call fgrid_advance(time,delt) 
      
      dtlevnew = rinfinity
      cfl_level = 0.d0    !# to keep track of max cfl seen on each level
c 
c each grid needs to get extra storage for integration serially (for now at least) in case
c dynamic memory moves the alloc array during expansion.  so get it once and for all for
c biggest chunk it may need. That uses max1d on a side, which for parallel apps.
c is generally much smaller than serial, though this could be a waste.
c remember that max1d includes ghost cells
c
c if necessary can have each thread use max of grids it owns, but then
c cant use dynamic grid assignments.
c
c next loop will get enough storage, it wont be dimensioned as called for
c grids smaller than max1d on a side.
      do j = 1, maxthreads
        
        locfp_save(j) = igetsp(2*max1d*max1d*nvar)
        locgp_save(j) = igetsp(2*max1d*max1d*nvar)

      end do

!$OMP PARALLEL DO PRIVATE(j,locold, locnew, mptr,nx,ny,mitot,mjtot)  
!$OMP&            PRIVATE(locfp, locfm, locgp,locgm,ntot,xlow,ylow)
!$OMP&            PRIVATE(locaux,lenbc,locsvf,locsvq,locx1d,i)
!$OMP&            PRIVATE(dtnew,time,mythread)
!$OMP&            SHARED(rvol,rvoll,level,nvar,mxnest,alloc,intrat,delt)
!$OMP&            SHARED(nghost,intratx,intraty,hx,hy,naux,listsp)
!$OMP&            SHARED(node,rnode,dtlevnew,numgrids,listgrids)
!$OMP&            SHARED(locfp_save,locgp_save)
!$OMP&            SCHEDULE (DYNAMIC,1)
!$OMP&            DEFAULT(none)
      do  j = 1, numgrids(level)
c          mptr   = mget(j,level)
          mptr = listgrids(j)
          locold = node(store2, mptr)
          locnew = node(store1, mptr)
          nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
          ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
          time   = rnode(timemult,mptr)
c
          mitot  = nx + 2*nghost
          mjtot  = ny + 2*nghost

c         ::: get scratch storage for fluxes and slopes
!--          locfp = igetsp(mitot*mjtot*nvar)
!--          locfm = igetsp(mitot*mjtot*nvar)
!--          locgp = igetsp(mitot*mjtot*nvar)
!--          locgm = igetsp(mitot*mjtot*nvar)
c
c next way so dont call igetsp so much, less parallel bottleneck in critical section
!--           locfp = igetsp(2*mitot*mjtot*nvar)
!--           locfm = locfp + mitot*mjtot*nvar
!--           locgp = igetsp(2*mitot*mjtot*nvar)
!--           locgm = locgp + mitot*mjtot*nvar
c next way for dynamic memory enlargement safety
!$         mythread = omp_get_thread_num()
           locfp = locfp_save(mythread+1)
           locfm = locfp + mitot*mjtot*nvar
           locgp = locgp_save(mythread+1)
           locgm = locgp + mitot*mjtot*nvar

c
c  copy old soln. values into  next time step's soln. values
c  since integrator will overwrite it. only for grids not at
c  the finest level. finest level grids do not maintain copies
c  of old and new time solution values.
c
          if (level .lt. mxnest) then
             ntot   = mitot * mjtot * nvar
cdir$ ivdep
             do 10 i = 1, ntot
 10            alloc(locold + i - 1) = alloc(locnew + i - 1)
          endif
c
      xlow = rnode(cornxlo,mptr) - nghost*hx
      ylow = rnode(cornylo,mptr) - nghost*hy

!$OMP CRITICAL(rv)
      rvol = rvol + nx * ny
      rvoll(level) = rvoll(level) + nx * ny
!$OMP END CRITICAL(rv)


      locaux = node(storeaux,mptr)
c
      if (node(ffluxptr,mptr) .ne. 0) then
         lenbc  = 2*(nx/intratx(level-1)+ny/intraty(level-1))
         locsvf = node(ffluxptr,mptr)
         locsvq = locsvf + nvar*lenbc
         locx1d = locsvq + nvar*lenbc
         call qad(alloc(locnew),mitot,mjtot,nvar,
     1            alloc(locsvf),alloc(locsvq),lenbc,
     2            intratx(level-1),intraty(level-1),hx,hy,
     3            naux,alloc(locaux),alloc(locx1d),delt,mptr)
      endif

c        # see if the grid about to advanced has gauge data to output
c        # this corresponds to previous time step, but output done
c        # now to make linear interpolation easier, since grid
c        # now has boundary conditions filled in.
c        # no testing here for mgauges>0 so that do not
c        # need to use gauges.i. the only time advanc is
c        # called that isn't "real" is in the initial setting
c        # up of grids (setgrd), but source grids are 0 there so
c        # nothing will be output.

c    should change the way  dumpguage does io - right now is critical section
           call dumpgauge(alloc(locnew),alloc(locaux),xlow,ylow,
     .                    nvar,mitot,mjtot,naux,mptr)

c
      call stepgrid(alloc(locnew),alloc(locfm),alloc(locfp),
     1            alloc(locgm),alloc(locgp),
     2            mitot,mjtot,nghost,
     3            delt,dtnew,hx,hy,nvar,
     4            xlow,ylow,time,mptr,naux,alloc(locaux))

      if (node(cfluxptr,mptr) .ne. 0)
     1   call fluxsv(mptr,alloc(locfm),alloc(locfp),
     2               alloc(locgm),alloc(locgp),
     3               alloc(node(cfluxptr,mptr)),mitot,mjtot,
     4               nvar,listsp(level),delt,hx,hy)
      if (node(ffluxptr,mptr) .ne. 0) then
         lenbc = 2*(nx/intratx(level-1)+ny/intraty(level-1))
         locsvf = node(ffluxptr,mptr)
         call fluxad(alloc(locfm),alloc(locfp),
     1               alloc(locgm),alloc(locgp),
     2               alloc(locsvf),mptr,mitot,mjtot,nvar,
     4               lenbc,intratx(level-1),intraty(level-1),
     5               nghost,delt,hx,hy)
      endif
c
!--          call reclam(locfp, mitot*mjtot*nvar)
!--          call reclam(locfm, mitot*mjtot*nvar)
!--          call reclam(locgp, mitot*mjtot*nvar)
!--          call reclam(locgm, mitot*mjtot*nvar)
c         next way to reclaim was to minimize calls to
c         reclam, due to critical section and openmp
!--          call reclam(locfp, 2*mitot*mjtot*nvar)
!--          call reclam(locgp, 2*mitot*mjtot*nvar)

!$OMP CRITICAL (newdt)
          dtlevnew = dmin1(dtlevnew,dtnew)
!$OMP END CRITICAL (newdt)    
c
          rnode(timemult,mptr)  = rnode(timemult,mptr)+delt
      end do
!$OMP END PARALLEL DO
c
c     debug statement:
c     write(*,*)" from advanc: level ",level," dtlevnew ",dtlevnew

c new way to reclaim for safety with dynamic memory and openmp
      do j = 1, maxthreads
        call reclam(locfp_save(j),2*max1d*max1d*nvar)
        call reclam(locgp_save(j),2*max1d*max1d*nvar)
      end do      
c
      return
      end
c
c -------------------------------------------------------------
c
       subroutine prepgrids(listgrids,num, level)

       implicit double precision (a-h,o-z)
       include "call.i"
       integer listgrids(num)

       mptr = lstart(level)
       do j = 1, num
          listgrids(j) = mptr
          mptr = node(levelptr, mptr)
       end do

      if (mptr .ne. 0) then
         write(*,*)" Error in routine setting up grid array "
         stop
      endif

      return
      end
