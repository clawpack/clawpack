c
c ----------------------------------------------------------------
c
       program amr2ez
c
c  Use adaptive mesh refinement to solve the hyperbolic 2-d equation:
c
c              u  +  f(u)    + g(u)   = 0
c               t         x        y
c
c or the more general non-conservation law form:
c              u  +  A u     + B u    = 0
c               t         x        y
c
c  using the wave propagation method as in CLAWPACK in combination
c  with the locally uniform embedded grids of AMR.

c  Estimate error with Richardson extrap. (in errest.f)
c  + gradient checking (in errsp.f).  Initial conditions set
c  in (qinit.f), b.c.'s in (physbd.f).

c  Specify rectangular domain from
c           (xlower,ylower) to (xupper,yupper).
c
c  No rotated rectangles are used in this version.
c  Periodic b.c.'s finally implemented.
c
c =========================================================================
c  Copyright 1996,  Marsha J. Berger and Randall J. LeVeque
c
c  This software is made available for research and instructional use only.
c  You may copy and use this software without charge for these non-commercial
c  purposes, provided that the copyright notice and associated text is
c  reproduced on all copies.  For all other uses (including distribution of
c  modified versions), please contact the author at the address given below.
c
c  *** This software is made available "as is" without any assurance that it
c  *** will work for your purposes.  The software may in fact have defects, so
c  *** use the software at your own risk.
c
c  --------------------------------------
c    AMRCLAW Version 0.4,  June, 1999
c     compatible with CLAWPACK Version 4.0
c    Homepage: http://www.amath.washington.edu/~claw/
c  --------------------------------------
c
c   Authors:
c
c             Marsha J. Berger
c             Courant Institute of Mathematical Sciences
c             New York University
c             251 Mercer St.
c             New York, NY 10012
c             berger@cims.nyu.edu
c
c             Randall J. LeVeque
c             Applied Mathematics
c             Box 352420
c             University of Washington,
c             Seattle, WA 98195-2420
c             rjl@amath.washington.edu
c
c =========================================================================
c


c
c ----------------------------------------------------------------
c
      use geoclaw_module
      use topo_module

      implicit double precision (a-h,o-z)

      include "call.i"
      include "regions.i"
      include "gauges.i"

      common /combc2/ mthbc(4)
      common /comfine/ dxmin,dymin

      character * 12     pltfile,infile,outfile,dbugfile,matfile
      character * 20     parmfile
      character * 25 fname
      logical            vtime,rest
      dimension          tout(maxout)
      dimension          tchk(maxout)

      integer oldmode
c
c
c  you may want to turn this on for SUN workstation, or replace
c  set to signal on overflow, divide by zero, and illegal operation
c
c       oldmode = ieee_handler("set","common",SIGFPE_ABORT)
c       if (oldmode .ne. 0) then
c           write(outunit,*)' could not set ieee trapper '
c           write(*,*)      ' could not set ieee trapper '
c           stop
c        endif
c
      infile   = 'amr2ez.data'
      outfile  = 'fort.amr'
      pltfile  = 'fort.ncar'
      dbugfile = 'fort.debug'
      matfile  = 'fort.nplot'
      parmfile  = 'fort.parameters'

      open(outunit, file=outfile,status='unknown',form='formatted')
      open(dbugunit,file=dbugfile,status='unknown',form='formatted')
      open(parmunit,file=parmfile,status='unknown',form='formatted')
c
c     ## New in 4.4:
c     ## Open file and skip over leading lines with # comments:
      fname = 'amr2ez.data'
      call opendatafile(inunit,fname)
c
      open(10,file='fort.info',status='unknown',form='formatted')
c
c
c     Number of space dimensions:  ## New in 4.4
c     read(55,*) ndim  ## not yet in amrclaw
c
c     ## The remainder is unchanged from 4.3:

c
c     domain variables
      read(inunit,*) nx
      read(inunit,*) ny
      read(inunit,*) mxnest
      if (abs(mxnest) .gt. maxlv) then
         write(outunit,*)
     &    'Error ***   mxnest > max. allowable levels (maxlv) in common'
         write(*,*)
     &    'Error ***   mxnest > max. allowable levels (maxlv) in common'
         stop
      endif
      if (mxnest .lt. 0) then
c         # mxnext<0 flags anisotropic refinement - read in refinement
c         # ratios in x,y,t from next three lines:
          mxnest = -mxnest
          read(inunit,*) (intratx(i),i=1,max(1,mxnest-1))
          read(inunit,*) (intraty(i),i=1,max(1,mxnest-1))
          read(inunit,*) (kratio(i), i=1,max(1,mxnest-1))
      else
          read(inunit,*) (intratx(i),i=1,max(1,mxnest-1))
          do i=1,max(1,mxnest-1)
             intraty(i) = intratx(i)
             kratio(i)  = intratx(i)
          enddo
      endif

c      if (intrat(1) .lt. 0) then
c          # this flags the situation where we do not want to refine in t
c           write(6,*) '*** No refinement in time!'
c           intrat(1) = -intrat(1)
c           do i=1,max(1,mxnest-1)
c              kratio(i) = 1
c              enddo
c        else
c          # the normal situation is to refine dt by same factors as dx,dy
c           do i=1,max(1,mxnest-1)
c              kratio(i) = intrat(i)
c              enddo
c        endif


      read(inunit,*) nout
      if (nout .gt. maxout) then
         write(outunit,*) 'Error ***   nout > maxout in common'
         write(*,*)       'Error ***   nout > maxout in common'
         stop
      endif
      read(inunit,*) outstyle
        if (outstyle.eq.1) then
           read(inunit,*) tfinal
c          # array tout is set below after reading t0
           iout = 0
           endif
        if (outstyle.eq.2) then
           read(inunit,*) (tout(i), i=1,nout)
           iout = 0
           endif
        if (outstyle.eq.3) then
           read(inunit,*) iout,nstop
           nout = 0
           endif
      read(inunit,*) possk(1)
      read(inunit,*) dtv2
      read(inunit,*) cflv1
      read(inunit,*) cfl
      read(inunit,*) nv1
      if (outstyle.eq.1 .or. outstyle.eq.2) then
         nstop = nv1
      endif

      read(inunit,*) method(1)
      vtime = (method(1) .eq. 1)
      read(inunit,*) method(2)
      iorder = method(2)
      read(inunit,*) method(3)
      if (method(3) .lt. 0) then
         write(6,*) '*** ERROR ***  method(3) < 0'
         write(6,*) '    dimensional splitting not supported in amrclaw'
         stop
         endif

      read(inunit,*) method(4)
      read(inunit,*) method(5)
      read(inunit,*) mcapa1
      read(inunit,*) naux
      if (naux .gt. maxaux) then
         write(outunit,*) 'Error ***   naux > maxaux in common'
         write(*,*)       'Error ***   naux > maxaux in common'
         stop
      endif
      do iaux = 1, naux
         read(inunit,*) auxtype(iaux)
      end do

      read(inunit,*) nvar
      read(inunit,*) mwaves
      if (mwaves .gt. maxwave) then
         write(outunit,*) 'Error ***   mwaves > maxwave in common'
         write(*,*)       'Error ***   mwaves > maxwave in common'
         stop
      endif
      read(inunit,*) (mthlim(mw), mw=1,mwaves)

      read(inunit,*) t0
      read(inunit,*) xlower
      read(inunit,*) xupper
      read(inunit,*) ylower
      read(inunit,*) yupper

      read(inunit,*) nghost
      read(inunit,*) mthbc(1)
      read(inunit,*) mthbc(2)
      read(inunit,*) mthbc(3)
      read(inunit,*) mthbc(4)

c     1 = left, 2 = right 3 = bottom 4 = top boundary
      xperdom = (mthbc(1).eq.2 .and. mthbc(2).eq.2)
      yperdom =  (mthbc(3).eq.2 .and. mthbc(4).eq.2)
      spheredom =  (mthbc(3).eq.5 .and. mthbc(4).eq.5)

      if (spheredom) then
         write(6,*) '*** Error: spherical domain not yet fully '
         write(6,*) '    implemented in GeoClaw'
	 stop
	 endif

      if ((mthbc(1).eq.2 .and. mthbc(2).ne.2) .or.
     &    (mthbc(2).eq.2 .and. mthbc(1).ne.2)) then
         write(6,*) '*** ERROR ***  periodic boundary conditions: '
         write(6,*) '  mthbc(1) and mthbc(2) must BOTH be set to 2'
         stop
         endif

      if ((mthbc(3).eq.2 .and. mthbc(4).ne.2) .or.
     &    (mthbc(4).eq.2 .and. mthbc(3).ne.2)) then
         write(6,*) '*** ERROR ***  periodic boundary conditions: '
         write(6,*) '  mthbc(3) and mthbc(4) must BOTH be set to 2'
         stop
         endif

      if ((mthbc(3).eq.5 .and. mthbc(4).ne.5) .or.
     &    (mthbc(4).eq.5 .and. mthbc(3).ne.5)) then
         write(6,*) '*** ERROR ***  sphere bcs at top and bottom: '
         write(6,*) '  mthbc(3) and mthbc(4) must BOTH be set to 5'
         stop
         endif

      if (spheredom .and. .not. xperdom) then
         write(6,*)'*** ERROR ***  sphere bcs at top and bottom: '
         write(6,*)'must go with periodic bcs at left and right  '
         stop
         endif

      if (outstyle.eq.1) then
	   do i=1,nout
	      tout(i) = t0 + i*(tfinal-t0)/float(nout)
	      enddo
           endif


c     restart and checkpointing
      read(inunit,*) rest
      read(inunit,*) ichkpt
      if (ichkpt .lt. 0) then
         if (-ichkpt .gt. maxout) then
            write(6,*) 'Error -ichkpt can be no greater than maxout'
            stop
            endif
         read(inunit,*) (tchk(i), i=1,-ichkpt)
         endif
c
c     refinement variables
      read(inunit,*) tol
      read(inunit,*) tolsp
      read(inunit,*) kcheck
      read(inunit,*) ibuff
      read(inunit,*) cut
c
c     style of output
c
      read(inunit,*) printout
      read(inunit,*) ncarout
      read(inunit,*) matlabout

c
c
c     # read verbose/debugging flags
      read(inunit,*) dprint
      read(inunit,*) eprint
      read(inunit,*) edebug
      read(inunit,*) gprint
      read(inunit,*) nprint
      read(inunit,*) pprint
      read(inunit,*) rprint
      read(inunit,*) sprint
      read(inunit,*) tprint
      read(inunit,*) uprint
      close(inunit)

c     # look for capacity function via auxtypes:
      mcapa = 0
      do 20 iaux = 1, naux
        if (auxtype(iaux) .eq. "capacity") then
          if (mcapa .ne. 0) then
            write(*,*)" only 1 capacity allowed"
            stop
          else
            mcapa = iaux
          endif
        endif

c       # change to Version 4.1 terminology:
        if (auxtype(iaux) .eq. "leftface") auxtype(iaux) = "xleft"
        if (auxtype(iaux) .eq. "bottomface") auxtype(iaux) = "yleft"

        if (.not. (auxtype(iaux) .eq."xleft" .or.
     .             auxtype(iaux) .eq. "yleft".or.
     .             auxtype(iaux) .eq. "capacity".or.
     .             auxtype(iaux) .eq. "center"))  then
                  write(*,*)" unknown type for auxiliary variables"
                  write(*,*)" use  xleft/yleft/center/capacity"
                  stop
        endif
   20   continue

c     ::: error checking of input data :::::::::::::::::::::::

      if (mcapa .ne. mcapa1) then
         write(outunit,*) 'Error ***  mcapa does not agree with auxtype'
         write(*,*) 'Error ***  mcapa does not agree with auxtype'
         stop
      endif
      if (nvar .gt. maxvar) then
         write(outunit,*) 'Error ***   nvar > maxvar in common'
         write(*,*)       'Error ***   nvar > maxvar in common'
         stop
      endif
      if (2*nghost .gt. min(nx,ny) .and. ny.ne.1) then
         mindim = 2*nghost
         write(outunit,*) 'Error ***   need finer domain >',
     .         mindim, ' cells'
         write(*,*)       'Error ***   need finer domain >',
     .         mindim, ' cells'
c        stop
      endif
      if (mcapa .gt. naux) then
         write(outunit,*) 'Error ***   mcapa > naux in input file'
         write(*,*)       'Error ***   mcapa > naux in input file'
         stop
      endif

      if (.not. vtime .and. nout .ne. 0) then
         write(outunit,*)' cannot specify output times with fixed dt'
         write(*,*)      ' cannot specify output times with fixed dt'
         stop
      endif
c
c
      if (ncarout)
     .   open(pltunit1,file=pltfile,status='unknown',form='formatted')
c
c
c     # write out parameters to fort.parameters....  still to be added
c
      write(parmunit,*) ' '
      write(parmunit,*) 'Running GeoClaw2 with parameter values:'
      write(parmunit,*) ' '


      write(6,*) ' '
      write(6,*) 'Running GeoClaw2 ...  '
      write(6,*) ' '

c     # default values of parameters that may be reset if user's setprob
c     # routine calls settopo, setdtopo, setqinit, setregions or setgauges.
      mgauges = 0
      mtopofiles = 0
      mregions = 0
      idtopo = 0
      iqinit = 0
      icoordsys = 0
      coeffmanning = 0.d0
      frictiondepth = 0.d0
c
      matlabu   = 0
      hxposs(1) = (xupper - xlower) / nx
      hyposs(1) = (yupper - ylower) / ny
      cflmax = 0.d0
c
c
c
      if (rest) then
         call restrt(nsteps,time,nvar)
         nstart  = nsteps
         tstart = time
c        # call user routine to set up problem parameters:
         call setprob()
         write(6,*) ' '
         write(6,*) 'Restarting from previous run'
         write(6,*) '   at time = ',time
         write(6,*) ' '
      else
c        # call user routine to set up problem parameters:
         call setprob()
         lentot = 0
         lenmax = 0
         lendim = 0
         rvol   = 0.0d0
         do 8 i   = 1, mxnest
 8         rvoll(i) = 0.0d0
         evol   = 0.0d0
         call   stst1

c        # changed 4/24/09: store dxmin,dymin for setaux before
c        # grids are made, in order to average up from finest grid.
         dxmin = hxposs(mxnest)
         dymin = hyposs(mxnest)

         call   domain (nvar,vtime,nx,ny,naux,t0)
c        # hold off on gauges until grids are set. the fake call to advance at the very
c        # first timestep looks at the gauge array but it is not yet built
         mgaugeSave = mgauges
         mgauges = 0
         call   setgrd (nvar,cut,naux,dtinit,t0)
         mgauges = mgaugeSave

!--         if (possk(1) .gt. dtinit*cflv1/cfl .and. vtime) then
!--c        ## initial time step was too large. reset to dt from setgrd
!--              write(6,*) "*** Initial time step reset for desired cfl"
!--              possk(1) = dtinit
!--              do i = 2, mxnest-1
!--                 possk(i) = possk(i-1)*kratio(i-1)
!--              end do
!--         endif

         time = t0
         nstart = 0
      endif

      if (icoordsys .eq. 0) then
         write(6,*) 'ERROR:  icoordsys is not set properly'
         write(6,*) '        perhaps you neglected to call setgeo?'
         stop
         endif

      tstart = time

      write(parmunit,*) ' '
      write(parmunit,*) '--------------------------------------------'
      write(parmunit,*) ' '
      write(parmunit,*) '   rest = ', rest, '   (restart?)'
      write(parmunit,*) '   tstart = ',tstart
      write(parmunit,*) ' '
c
c  print out program parameters for this run
c
      write(outunit,107)tol,tolsp,iorder,kcheck,ibuff,nghost,cut,
     1            mxnest,ichkpt,cfl
      write(outunit,109) xupper,yupper,xlower,ylower,nx,ny
      write(outunit,139)(intratx(i),i=1,mxnest)
      write(outunit,139)(intraty(i),i=1,mxnest)
      write(outunit,119) naux
      write(outunit,129) (iaux,auxtype(iaux),iaux=1,naux)
      if (mcapa .gt. 0) write(outunit,149) mcapa
107   format(/
     *       ' amrclaw parameters:',//,
     *       ' error tol            ',e12.5,/,
     *       ' spatial error tol    ',e12.5,/,
     *       ' order of integrator     ',i9,/,
     *       ' error checking interval ',i9,/,
     *       ' buffer zone size        ',i9,/,
     *       ' nghost                  ',i9,/,
     *       ' volume ratio cutoff  ',e12.5,/,
     *       ' max. refinement level   ',i9,/,
     *       ' user sub. calling times ',i9,/,
     *       ' cfl # (if var. delt) ',e12.5,/)
109   format(' xupper(upper corner) ',e12.5,/,
     *       ' yupper(upper corner) ',e12.5,/,
     *       ' xlower(lower corner) ',e12.5,/,
     *       ' ylower(lower corner) ',e12.5,/,
     *       ' nx = no. cells in x dir.',i9,/,
     *       ' ny = no. cells in y dir.',i9,/,
     *       ' refinement ratios       ',6i5,/,/)
119   format(' no. auxiliary vars.     ',i9)
129   format('       var ',i5,' of type ', a10)
139   format(' refinement ratios:       ',6i5,/)
149   format(' capacity fn. is aux. var',i9)

      write(6,*) ' '
      write(6,*) 'Done reading data, starting computation ...  '
      write(6,*) ' '


      call outtre (mstart,printout,nvar,naux)
      write(outunit,*) "  original total mass ..."
      call conck(1,nvar,time)
      call valout(1,lfine,time,nvar,naux)
c
c     --------------------------------------------------------
c     # tick is the main routine which drives the computation:
c     --------------------------------------------------------
      call tick(nvar,iout,nstart,nstop,cut,vtime,time,ichkpt,naux,
     &          nout,tout,tchk,t0)
c     --------------------------------------------------------

c     # Done with computation, cleanup:


      lentotsave = lentot
      call cleanup(nvar,naux)
      if (lentot .ne. 0) then
        write(outunit,*) lentot," words not accounted for ",
     &                   "in memory cleanup"
        write(*,*)        lentot," words not accounted for ",
     &                   "in memory cleanup"
      endif
c
c report on statistics
c
      open(matunit,file=matfile,status='unknown',form='formatted')
      write(matunit,*) matlabu-1
      write(matunit,*) mxnest
      close(matunit)

      write(outunit,*)
      write(outunit,*)
      write(outunit,901) lentotsave
      write(outunit,902) lenmax
      write(outunit,903) lendim

      write(outunit,904) rvol
      do 60 level = 1,mxnest
 60     write(outunit,905) level, rvoll(level)

      write(outunit,906) evol
      if (evol+rvol .gt. 0.) then
         ratmet = rvol / (evol+rvol) * 100.0d0
      else
         ratmet = 0.0d0
      endif
      write(outunit,907) ratmet
      write(outunit,908) cflmax

 901  format('current  space usage = ',i12)
 902  format('maximum  space usage = ',i12)
 903  format('need space dimension = ',i12,/)
 904  format('number of cells advanced for time integration = ',f20.6)
 905  format(3x,'# cells advanced on level ',i4,' = ',f20.2)
 906  format('number of cells advanced for error estimation = ',f20.6,/)
 907  format(' percentage of cells advanced in time  = ', f10.2)
 908  format(' maximum Courant number seen = ', f10.2)
c
      write(outunit,909)
 909  format(//,' ------  end of AMRCLAW integration --------  ')
c
c     # Close output and debug files.
      close(outunit)
      close(dbugunit)
c
      stop
      end
