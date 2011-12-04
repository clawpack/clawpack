
c=========================================================================
      subroutine setqinit
c=========================================================================
c
c     Read initial perturbation from setqinit.data
c
c     For problems where h = q(i,j,1) is depth of fluid, initialized to 
c     sealevel in qinit.
c
c     If iqinit = 0, no initial perturbation to sealevel.
c
c     If iqinit = 1,2 or 3 initial perturbation is read in.
c     Then the qinit file should have 3 columns:
c              x   y  z  if icoordsys = 1 
c            lat lon  z if icoordsys = 2 

c     If iqinit = 4, file defines eta, and h =q(i,j,1)= max(eta-b,0).
c
c     Longitude and latitude advance in the standard GIS way from 
c     upper left corner across in x and then down in y.


      implicit double precision (a-h,o-z)
      character*25 fname
      character*150 qinitfname
      logical foundFile

      include "call.i"
      include "qinit.i"

      write(parmunit,*) ' '
      write(parmunit,*) '--------------------------------------------'
      write(parmunit,*) 'SETQINIT:'
      write(parmunit,*) '-------------'

      fname  = 'setqinit.data'
      inquire(file=fname,exist=foundFile)
      if (.not. foundFile) then
        write(*,*) 'You must provide a file ', fname
        stop
      endif

      iunit = 7
      call opendatafile(iunit, fname)

      read(7,*) iqinit

      if (iqinit.eq.0) then
c         # no perturbation specified
          write(parmunit,*)  '  iqinit = 0, no perturbation'
          write(6,*) '  iqinit = 0, no perturbation'
          return
          endif

      read(7,*) qinitfname
      read(7,*) minlevelqinit, maxlevelqinit
      close(unit=7)

      write(parmunit,*) '   minlevel, maxlevel, qinitfname:'
      write(parmunit,*)  minlevelqinit, maxlevelqinit, qinitfname
c      write(6,*) '   minlevel, maxlevel, qinitfname:'
c      write(6,*)  minlevelqinit, maxlevelqinit, qinitfname

      call readqinit(mxqinit,myqinit,dxqinit,dyqinit,
     &      xlowqinit,xhiqinit,ylowqinit,yhiqinit,qinitfname)

      mqinit=mxqinit*myqinit

c===================Check qinit Memory=============================
      if (mqinit.gt.maxqinitsize) then
          write(*,*) 'SETQINIT: not enough workspace'
          write(*,*) 'Increase maxqinitsize to', mqinit
          write(*,*) 'in qinit.i, then recompile'
          Stop
        else
          write(parmunit,*) 
     &        '  qinit space allocated:' ,maxqinitsize
          write(parmunit,*) '  qinit space used:',mqinit
        endif

      return
      end
