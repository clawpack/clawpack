c=========================================================================
      subroutine setfixedgrids
c=========================================================================

      implicit double precision (a-h,o-z)
      character*25, parameter :: fname = 'setfixedgrids.data'
      logical foundFile

      include "fixedgrids.i"
      include "call.i"

      write(parmunit,*) ' '
      write(parmunit,*) '--------------------------------------------'
      write(parmunit,*) 'SETFIXEDGRIDS:'
      write(parmunit,*) '-----------'

      inquire(file=fname,exist=foundFile)
      if (.not. foundFile) then
        write(*,*) 'You must provide a file ', fname
        stop
      endif

      iunit = 7
      call opendatafile(iunit, fname)

      read(7,*) mfgrids
      if (mfgrids.gt.maxfgrids) then
           write(*,*) 'SETFIXEDGRIDS: ERROR mfgrids > maxfgrids'
           write(*,*) 'Decrease the number of fixed grids or'
           write(*,*) 'Increase maxfgrids in fixedgrids.i'
           stop
           endif

      if (mfgrids .eq. 0) then
         write(parmunit,*) '  No fixed grids specified for output'
         return
         endif

      do i=1,mfgrids
         read(7,*) tstartfg(i),tendfg(i),noutfg(i),xlowfg(i),xhifg(i),
     &             ylowfg(i),yhifg(i),mxfg(i),myfg(i), 
     &             ioutarrivaltimes(i), ioutsurfacemax(i)

      enddo
      close(7)

c     # set some parameters for each grid
      do i=1,mfgrids
c       # set dtfg (the timestep length between outputs) for each grid
        if (tendfg(i).le.tstartfg(i)) then
           if (noutfg(i).gt.1) then 
             write(*,*) 'SETFIXEDGRIDS: ERROR for fixed grid', i
             write(*,*) 'tstartfg=tendfg yet noutfg>1'
             write(*,*) 'set tendfg > tstartfg or set noutfg = 1'
             stop
           else
              dtfg(i)=0.d0
           endif
        else
           if (noutfg(i).lt.2) then
              write(*,*) 'SETFIXEDGRIDS: ERROR for fixed grid', i
              write(*,*) 'tendfg>tstartfg, yet noutfg=1'
              write(*,*) 'set noutfg > 2'
              stop
           else
              dtfg(i)=(tendfg(i)-tstartfg(i))/(noutfg(i)-1)
           endif
        endif

c       # initialize tlastoutfg and ilastoutfg
        tlastoutfg(i)= tstartfg(i)-dtfg(i)
        ilastoutfg(i)= 0

c       # set spatial intervals dx and dy on each grid
        if (mxfg(i).gt.1) then
           dxfg(i)= (xhifg(i)-xlowfg(i))/(mxfg(i)-1)
        elseif (mxfg(i).eq.1) then
           dxfg(i)=0.d0
        else
           write(*,*) 'SETFIXEDGRIDS: ERROR for fixed grid', i
           write(*,*) 'x grid points mxfg<=0, set mxfg>= 1'
        endif

        if (myfg(i).gt.1) then
           dyfg(i)= (yhifg(i)-ylowfg(i))/(myfg(i)-1)
        elseif (myfg(i).eq.1) then
           dyfg(i)=0.d0
        else
           write(*,*) 'SETFIXEDGRIDS: ERROR for fixed grid', i
           write(*,*) 'y grid points myfg<=0, set myfg>= 1'
        endif 
      enddo

c     # set the number of variables stored for each grid
c     # this should be (the number of variables you want to write out + 1)
      do i=1, mfgrids 
         mfgridvars(i)  = 6  
         mfgridvars2(i) = 3*ioutsurfacemax(i) + ioutarrivaltimes(i)
      enddo

c      # find entry point into work arrays for each fixed grid
c      # make sure enough space has been alotted for fixed grids in memory
       i0fg(1)=1
       i0fg2(1)=1
       do i=2,mfgrids
         i0fg(i)= i0fg(i-1)   + mfgridvars(i-1)*mxfg(i-1)*myfg(i-1)
         i0fg2(i)= i0fg2(i-1) + mfgridvars2(i-1)*mxfg(i-1)*myfg(i-1)
       enddo
       mspace=i0fg(mfgrids) + 
     &          mfgridvars(mfgrids)*mxfg(mfgrids)*myfg(mfgrids)
       mspace2=i0fg2(mfgrids) + 
     &          mfgridvars2(mfgrids)*mxfg(mfgrids)*myfg(mfgrids)
       mspace=mspace+mspace2
       if (mspace.gt.maxfgridsize) then
         write(*,*) 'SETFIXEDGRIDS: ERROR not enough memory allocated'
         write(*,*) 'Decrease the number and size of fixed grids or'
         write(*,*) 'set maxfgridsize in fixedgrids.i to:', mspace
         stop
       endif

c      #initialize fixed grid work arrays to NaN
c      #this will prevent non-filled values from being misinterpreted
       do k=1,maxfgridsize
          fgridearly(k)=d_nan()
          fgridlate(k)=d_nan()
          fgridoften(k)=d_nan()
       enddo
       tcfmax=-1.d16
      write(parmunit,*) '  mfgrids = ',mfgrids

      do i=1,mfgrids
         write(parmunit,701) 
  701    format(2i4,6d12.3)
         enddo

      return
      end

      real*8 function d_nan()
      real*8 dnan
      integer inan(2)
      equivalence (dnan,inan)
      inan(1)=2147483647
      inan(2)=2147483647
      d_nan=dnan
      return
      end
