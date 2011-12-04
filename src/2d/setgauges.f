c=========================================================================
      subroutine setgauges
c=========================================================================

      implicit double precision (a-h,o-z)
      character*25, parameter :: fname = 'setgauges.data'
      logical foundFile

      include "gauges.i"
      include "call.i"

      write(parmunit,*) ' '
      write(parmunit,*) '--------------------------------------------'
      write(parmunit,*) 'SETGAUGES:'
      write(parmunit,*) '---------'

      inquire(file=fname,exist=foundFile)
      if (.not. foundFile) then
        write(*,*) 'You must provide a file ', fname 
        stop
      endif

      iunit = 7
      call opendatafile(iunit, fname)

      read(7,*) mgauges
      if (mgauges.gt.maxgauges) then
            write(*,*) 'ERROR in setgauges'
            write(*,*) 'mgauges = ',mgauges,'   maxgauges = ',maxgauges
            write(*,*) 'increase maxgauges in gauges.i'
            stop
            endif

      if (mgauges .eq. 0) then
         write(parmunit,*) 'No gauges specified'
         return
         endif


      write(parmunit,*) '  mgauges = ',  mgauges
      do i=1,mgauges
          read(7,*) igauge(i),xgauge(i),ygauge(i), 
     &              t1gauge(i), t2gauge(i)
          mbestsrc(i) = 0   ! initialize for starters
          enddo
      close(7)

c     # open file for output of gauge data
c     # all data is output in one binary file with format
c     # gauge #, time, height

      open(unit=OUTGAUGEUNIT,file='fort.gauge',status='unknown',
     .           form='formatted')

      return
      end
