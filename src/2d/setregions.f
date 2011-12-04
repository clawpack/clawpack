c=========================================================================
      subroutine setregions
c=========================================================================

      implicit double precision (a-h,o-z)
      character*25 fname
      logical foundFile

      include "regions.i"
      include "call.i"

      write(parmunit,*) ' '
      write(parmunit,*) '--------------------------------------------'
      write(parmunit,*) 'SETREGIONS:'
      write(parmunit,*) '-----------'

      fname  = 'setregions.data'
      inquire(file=fname,exist=foundFile)
      if (.not. foundFile) then
        write(*,*) 'You must provide a file ', fname
        stop
      endif

      iunit = 7
      call opendatafile(iunit, fname)

      read(7,*) mregions
      if (mregions.gt.maxregions) then
           write(*,*) 'SETREGIONS: ERROR mregions > maxregions'
           write(*,*) 'Decrease the number of regions or'
           write(*,*) 'Increase maxregions in regions.i'
           stop
           endif

      if (mregions .eq. 0) then
         write(parmunit,*) '  No regions specified for refinement'
         return
         endif

      do i=1,mregions
         read(7,*) minlevelregion(i),maxlevelregion(i),
     &         tlowregion(i),thiregion(i),
     &         xlowregion(i),xhiregion(i),ylowregion(i),yhiregion(i)
         enddo

      close(7)

      write(parmunit,*) '  mregions = ',mregions
      write(parmunit,*) 
     &  '  minlevel, maxleve, tlow, thi, xlow, xhi, ylow, yhigh values:'
c    &   ' min max tlow           thi           xlow            xhi',
c    &   '         ylow           yhi'

      do i=1,mregions
         write(parmunit,701) minlevelregion(i),maxlevelregion(i),
     &         tlowregion(i),thiregion(i),
     &         xlowregion(i),xhiregion(i),ylowregion(i),yhiregion(i)
  701    format(2i4,6d12.3)
         enddo

      return
      end
