c
c  -------------------------------------------------------------
c
      subroutine tick(nvar,iout,nstart,nstop,cut,vtime,time,ichkpt,
     &                naux,nout,tout,tchk,t0)
c
      implicit double precision (a-h,o-z)

      include  "call.i"

      logical    vtime, dumpout, dumpchk, varRefTime/.true./
      dimension dtnew(maxlv), ntogo(maxlv), tlevel(maxlv)
      dimension tout(nout)
      dimension tchk(maxout)

c
c :::::::::::::::::::::::::::: TICK :::::::::::::::::::::::::::::
c  main driver routine.  controls:
c        integration  of all grids.
c        error estimation / regridding
c        output counting
c        updating of fine to coarse grids

c  parameters:
c     nstop   = # of coarse grid time steps to be taken
c     iout    = output interval every 'iout' coarse time steps
c               (if 0, not used - set to inf.)
c     vtime   = true for variable timestep, calculated each coarse step
c
c  integration strategy is to advance a fine grid until it catches
c  up to the coarse grid. this strategy is applied recursively.
c  coarse grid goes first.
c
c  nsteps: used to count how number steps left for a level to be
c          integrated before it catches up with the next coarser level.
c  ncycle: counts number of coarse grid steps = # cycles.
c
c  icheck: counts the number of steps (incrementing by 1
c          each step) to keep track of when that level should
c          have its error estimated and finer levels should be regridded.
c ::::::::::::::::::::::::::::::::::::;::::::::::::::::::::::::::
c
      integer verbosity_regrid     
c     # Eventually add this to new input params?
c     # For now use verbosity, or set to 0 to suppress printing regrid info:
      verbosity_regrid = method(4) 

      ncycle         = nstart
      call setbestsrc()     ! need at very start of run, including restart
      if ( iout .eq. 0) iout  = iinfinity
      if (ichkpt .eq. 0) ichkpt = iinfinity
      nextout = 1
      if (nout .gt. 0) then
         tfinal = tout(nout)
c        if this is a restart, make sure output times start after restart time
         if (nstart .gt. 0) then
            do ii = 1, nout
              if (tout(ii) .gt. time) then
                nextout = ii
                go to 2
              endif
            end do
  2         continue
         endif
      else
         tfinal  = rinfinity
      endif

      nextchk = 1
      if (nstart .gt. 0 .and. ichkpt.lt.0) then
c        if this is a restart, make sure chkpt times start after restart time
         do ii = 1, -ichkpt
           if (tchk(ii) .gt. time) then
              nextchk = ii
              go to 3
              endif
           enddo
  3      continue
         endif

      tlevel(1)      = time
      do 5 i       = 2, mxnest
       tlevel(i) = tlevel(1)
 5     continue
c
c  ------ start of coarse grid integration loop. ------------------
c
 20   if (ncycle .ge. nstop .or. time .ge. tfinal) goto 999

      if (nextout  .le. nout) then
         outtime       = tout(nextout)
      else
         outtime       = rinfinity
      endif

      if (nextchk  .le. -ichkpt) then
         chktime       = tchk(nextchk)
      else
         chktime       = rinfinity
      endif

c     if (time.lt.outtime .and. time + possk(1) .ge. outtime) then
      if (time.lt.outtime .and. time+1.001*possk(1) .ge. outtime) then
c     apr 2010 mjb: modified so allow slightly larger timestep to
c     hit output time exactly, instead of taking minuscule timestep
c     should still be stable since increase dt in only 3rd digit.
c        ## adjust time step  to hit outtime exactly, and make output
         oldposs = possk(1)
         possk(1) = outtime - time
c        write(*,*)" old possk is ", possk(1)
         diffdt = oldposs - possk(1)  ! if positive new step is smaller

          if (method(4).ge.2) then    ! verbosity test to
            write(*,122) diffdt,outtime  ! notify of change
 122        format(" Adjusting timestep by ",e10.3,
     .             " to hit output time of ",e12.6)
c           write(*,*)" new possk is ", possk(1)
            if (diffdt .lt. 0.) then ! new step is slightly larger
              pctIncrease = -100.*diffdt/oldposs   ! minus sign to make whole expr. positive
              write(*,123) pctIncrease
 123          format(" New step is ",e8.2," % larger.",
     .               "  Should still be stable")
             endif
       endif
         do 12 i = 2, mxnest
 12         possk(i) = possk(i-1) / kratio(i-1)
         nextout = nextout + 1
        dumpout = .true.
      else
        dumpout = .false.
      endif

      if (time.lt.chktime .and. time + possk(1) .ge. chktime) then
c        ## adjust time step  to hit chktime exactly, and do checkpointing
         possk(1) = chktime - time
         do 13 i = 2, mxnest
 13         possk(i) = possk(i-1) / kratio(i-1)
         nextchk = nextchk + 1
        dumpchk = .true.
      else
        dumpchk = .false.
      endif

c  take output stuff from here - put it back at end.
c
          level        = 1
          ntogo(level) = 1
          do 10 i = 1, maxlv
 10          dtnew(i)  = rinfinity
c
c     ------------- regridding  time?  ---------
c
c check if either
c   (i)  this level should have its error estimated before being advanced
c   (ii) this level needs to provide boundary values for either of
c        next 2 finer levels to have their error estimated.
c        this only affects two grid levels higher, occurs because
c        previous time step needs boundary vals for giant step.
c  no error estimation on finest possible grid level
c
 60       continue
          if (icheck(level) .ge. kcheck) then
               lbase = level
          else if (level+1 .ge. mxnest) then
               go to 90
          else if (icheck(level+1) .ge. kcheck) then
               lbase = level+1
          else if (level+2 .ge. mxnest) then
               go to 90
          else if (icheck(level+2) .ge. kcheck) then
               lbase = level+2
          else
               go to 90
          endif
          if (lbase .eq. mxnest .or. lbase .gt. lfine) go to 70
c
c regrid level 'lbase+1' up to finest level.
c level 'lbase' stays fixed.
c
          if (rprint) write(outunit,101) lbase
101       format(8h  level ,i5,32h  stays fixed during regridding )
          call regrid(nvar,lbase,cut,naux,t0)
          call setbestsrc()     ! need at every grid change
c         call conck(1,nvar,time)
c         call outtre(lstart(lbase+1),.true.,nvar,naux)
c note negative time to signal regridding output in plots
c         call valout(lbase,lfine,-tlevel(lbase),nvar,naux)
c
c  maybe finest level in existence has changed. reset counters.
c
          if (rprint .and. lbase .lt. lfine) then
             call outtre(lstart(lbase+1),.false.,nvar,naux)
          endif
 70       continue
          do 80  i  = lbase, lfine
 80          icheck(i) = 0
          do 81  i  = lbase+1, lfine
 81          tlevel(i) = tlevel(lbase)
c
c          MJB: modified to check level where new grids start, which is lbase+1
c               initial code incorrectly used level in the following lines
          if (verbosity_regrid.ge.lbase+1) then
c             mptr    = lstart(level)
              mptr    = lstart(lbase+1)
              numgrids = 0
              numcells = 0
              do while (mptr .gt. 0) 
                  numgrids = numgrids + 1
                  numcells = numcells + (iregend(level)-iregst(level)+1)
     &                       * (jregend(level)-jregst(level)+1)
                  mptr = node(levelptr, mptr)
                  enddo
              if (lbase+1.gt.1) then
c                # don't bother on level 1, which always covers full domain
c                # with specified mx, my.
                 write(6,1007) lbase+1,numgrids, numcells
 1007            format('   Regridded Level',i3,' and higher ...  Now',
     &               i4,' grid(s) with',i9,' cells on this level')
                 do levnew = lbase+1,lfine
                     write(6,1006) intratx(levnew-1),intraty(levnew-1),
     &                             kratio(levnew-1),levnew
 1006                format('   Refinement ratios...  in x:', i3, 
     &                 '  in y:',i3,'  in t:',i3,' for level ',i4)
                 end do
                 endif
              endif

c  ------- done regridding --------------------
c
c integrate all grids at level 'level'.
c
 90       continue


          call advanc(level,nvar,dtlevnew,vtime,naux)

c         # rjl modified 6/17/05 to print out *after* advanc and print cfl
c         # rjl & mjb changed to cfl_level, 3/17/10

          timenew = tlevel(level)+possk(level)
          if (tprint) then
              write(outunit,100)level,cfl_level,possk(level),timenew
              endif
          if (method(4).ge.level) then
              write(6,100)level,cfl_level,possk(level),timenew
              endif
100       format(' AMRCLAW: level ',i2,'  CFL = ',e8.3,
     &           '  dt = ',e10.4,  '  final t = ',e12.6)


c        # to debug individual grid updates...
c        call valout(level,level,time,nvar,naux)
c
c done with a level of integration. update counts, decide who next.
c
          ntogo(level)  = ntogo(level) - 1
          dtnew(level)  = dmin1(dtnew(level),dtlevnew)
          tlevel(level) = tlevel(level) + possk(level)
          icheck(level) = icheck(level) + 1
c
          if (level .lt. lfine) then
             level = level + 1
c            #  check if should adjust finer grid time step to start wtih
             if (((possk(level-1) - dtnew(level-1))/dtnew(level-1)) .gt.
     .            .05) then
                dttemp = dtnew(level-1)/kratio(level-1)
                ntogo(level) = (tlevel(level-1)-tlevel(level))/dttemp+.9
              else
                ntogo(level) = kratio(level-1)
              endif
             possk(level) = possk(level-1)/ntogo(level)
             go to 60
          endif
c
 105      if (level .eq. 1) go to 110
              if (ntogo(level) .gt. 0) then
c                same level goes again. check for ok time step
 106             if ((possk(level)-dtnew(level))/dtnew(level)
     .                .gt. .05)  then
                    write(6,*)" ***adjusting timestep for level ",level
                 write(6,*)"   old ntogo dt ",ntogo(level),possk(level)
c                   adjust time steps for this and finer levels
                    ntogo(level) = ntogo(level) + 1
                    possk(level) = (tlevel(level-1)-tlevel(level))/
     .                             ntogo(level)
                    if (varRefTime) then
                       kratio(level-1) = ceiling(possk(level-1) /
     .                                           possk(level))
                    endif
                 write(6,*)"   new ntogo dt ",ntogo(level),possk(level)
                    go to 106
                 endif
                 go to 60
              else
                 level = level - 1
                 call update(level,nvar)
              endif
          go to 105
c
c  --------------one complete coarse grid integration cycle done. -----
c
c      time for output?  done with the whole thing?
c
 110      continue
          time    = time   + possk(1)
          ncycle  = ncycle + 1
          call conck(1,nvar,time)

       if ((ichkpt.gt.0 .and. mod(ncycle,ichkpt).eq.0) 
     &      .or. dumpchk) then
                call check(ncycle,time,nvar,naux)
       endif

       if ((mod(ncycle,iout).eq.0) .or. dumpout) then
         call valout(1,lfine,time,nvar,naux)
         if (printout) call outtre(mstart,.true.,nvar,naux)
       endif

      if ( .not. vtime) go to 201   
c
c      ## adjust time steps if variable time step and/or variable refinement ratios in time
c
        if (.not. varRefTime) then
c
c         find new dt for next cycle (passed back from integration routine).
           do 115 i = 2, lfine
             ii = lfine+1-i
             dtnew(ii) = min(dtnew(ii),dtnew(ii+1)*kratio(ii))
 115       continue
           possk(1) = dtnew(1)
           do 120 i = 2, mxnest
 120         possk(i) = possk(i-1) / kratio(i-1)

        else  ! since refinement ratio in time can change need to set new timesteps in different order
c             ! use same alg. as when setting refinement when first make new fine grids
          possk(1) = dtnew(1)
          do 125 i = 2, lfine
             if (dtnew(i)  .gt. possk(i-1)) then
               kratio(i-1) = 1  ! cant have larger timestep than parent level
               possk(i)    = possk(i-1)
            else
               kratio(i-1) = ceiling(possk(i-1)/dtnew(i))  ! round up for stable integer ratio
               possk(i)    = possk(i-1)/kratio(i-1)        ! set exact timestep on this level
           endif
 125    continue
        
      endif

 201  go to 20
c
999   continue

c
c  # computation is complete to final time or requested number of steps
c
       if (ncycle .ge. nstop .and. tfinal .lt. rinfinity) then
c         # warn the user that calculation finished prematurely
          write(outunit,102) nstop
          write(6,102) nstop
  102     format('*** Computation halted after nv(1) = ',i8,
     &           '  steps on coarse grid')
          endif
c
c  # final output (unless we just did it above)
c
      if (.not. dumpout) then
         if (mod(ncycle,iout).ne.0) then
           call valout(1,lfine,time,nvar,naux)
           if (printout) call outtre(mstart,.true.,nvar,naux)
         endif
      endif

c  # checkpoint everything for possible future restart
c  # (unless we just did it based on ichkpt or dumpchk)
c
c  # don't checkpoint at all if user set ichkpt=0

      if (ichkpt .gt. 0) then
         if ((ncycle/ichkpt)*ichkpt .ne. ncycle) then
           call check(ncycle,time,nvar,naux)
         endif
       endif

      if (ichkpt .lt. 0) then
         if (.not. dumpchk) then
           call check(ncycle,time,nvar,naux)
         endif
       endif

      return
      end
