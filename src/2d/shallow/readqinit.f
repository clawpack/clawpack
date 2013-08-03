c=========================================================================
      subroutine readqinit(mx,my,dx,dy,xlow,xhi,ylow,yhi,fname)
c=========================================================================

      use amr_module
      use qinit_module
      implicit double precision (a-h,o-z)

      character*150 fname

      logical foundFile


      inquire(file=fname,exist=foundFile)
      if (.not. foundFile) then
        write(*,*) 'You must provide a file ', fname
        stop
      endif

      write(6,*) '  '
      write(6,*) 'Reading qinit data from file  ', fname
      write(6,*) '  '

      write(parmunit,*) '  '
      write(parmunit,*) 'Reading qinit data from'
      write(parmunit,*) fname
      write(parmunit,*) '  '

      mpoints=0
      mx=0
      open(unit=19,
     &     file=fname,
     &     status='unknown',form='formatted')

c
c     # currently only supports one file type:
c     # x,y,z values, one per line in standard order from NW corner to SE
c     # z is perturbation from standard depth h,hu,hv set in qinit_geo,
c     # if iqinit = 1,2, or 3 respectively.
c     # if iqinit = 4, the z column corresponds to the definition of the 
c     # surface elevation eta. The depth is then set as q(i,j,1)=max(eta-b,0)

      read(19,end=20,fmt=*) x0,y0,z0
      mpoints=mpoints+1
      mx=mx+1
      qinitwork(mpoints)=z0
      x=x0

      do while (.true.)
           read(19,end=20,fmt=*) xnew,ynew,znew
           mpoints=mpoints+1
           qinitwork(mpoints)=znew
           if (xnew.gt.x) then
              mx=mx + 1
           else
              mx=1
           endif
        x=xnew
        y=ynew
        enddo

 20   continue


      close(unit=19)


      if (mpoints.gt.maxqinitsize) then
          write(*,*) 'SETQINIT: not enough workspace'
          write(*,*) 'Increase maxqinitsize to', mpoints
          write(*,*) 'in qinit.i, then recompile'
          Stop
      endif

      my= mpoints/mx
      xf=x
      yf=y
      dx=(xf-x0)/(mx-1)
      dy=-(yf-y0)/(my-1)

c      !Changes with new bilinear interp DLG 10/08:
c       below changed for new bilinear interpolation of auxiliary data
c       values are now noded based, no longer considered cell-centered
c      xlow = x0-.5d0*dx
c      yhi = y0+.5d0*dy
c      xhi = xf+.5d0*dx
c      ylow = yf-.5d0*dy

      xlow = x0
      yhi = y0
      xhi = xf
      ylow = yf

      return
      end
