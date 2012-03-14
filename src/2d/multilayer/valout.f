c
c -----------------------------------------------------
c
      subroutine valout (lst, lend, time, nvar, naux)
c
      use multilayer_module, only: layers,rho
      use amr_module
      
      implicit double precision (a-h,o-z)
      character*10  matname1, matname2

      dimension eta(layers+1),h(layers),hu(layers),hv(layers)
      
      double precision :: wind_x,wind_y,pressure

      iadd(ivar,i,j)  = loc + ivar - 1 + nvar*((j-1)*mitot+i-1)
      iaddaux(iaux,i,j) = locaux + iaux-1 + naux*(i-1) +
     .                                      naux*mitot*(j-1)



c ::::::::::::::::::::::::::: VALOUT ::::::::::::::::::::::::::::::::::;
c graphics output of soln values for contour or surface plots.
c modified for GeoClaw to output the surface level along with q.
c    surface = q(1,i,j) + aux(1,i,j)
c :::::::::::::::::::::::::::::::::::::;:::::::::::::::::::::::::::::::;

c
c     ### MATLAB graphics output
c
      if (matlabout) then
c        ###  make the file names and open output files
         matname1 = 'fort.qxxxx'
         matname2 = 'fort.txxxx'
         matunit1 = 50
         matunit2 = 60
         nstp     = matlabu
         do 55 ipos = 10, 7, -1
            idigit = mod(nstp,10)
            matname1(ipos:ipos) = char(ichar('0') + idigit)
            matname2(ipos:ipos) = char(ichar('0') + idigit)
            nstp = nstp / 10
 55      continue
         open(unit=matunit1,file=matname1,status='unknown',
     .       form='formatted')

         level = lst
         ngrids = 0
 65      if (level .gt. lfine) go to 90
            mptr = lstart(level)
 70         if (mptr .eq. 0) go to 80
              ngrids  = ngrids + 1
              nx      = node(ndihi,mptr) - node(ndilo,mptr) + 1
              ny      = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
              loc     = node(store1, mptr)
              locaux  = node(storeaux,mptr)
              mitot   = nx + 2*nghost
              mjtot   = ny + 2*nghost
              write(matunit1,1001) mptr, level, nx, ny
 1001 format(i5,'                 grid_number',/,
     &       i5,'                 AMR_level',/,
     &       i5,'                 mx',/,
     &       i5,'                 my')


c  old        xcorn = rnode(cornxlo,mptr) - .5d0*hxposs(level)
c  old        ycorn = rnode(cornylo,mptr) - .5d0*hyposs(level)
              xlow = rnode(cornxlo,mptr)
              ylow = rnode(cornylo,mptr)
              write(matunit1,1002)
     &              xlow,ylow,hxposs(level),hyposs(level)
 1002 format(e18.8,'    xlow', /,
     &       e18.8,'    ylow', /,
     &       e18.8,'    dx', /,
     &       e18.8,'    dy',/)

      do j = nghost+1, mjtot-nghost
         do i = nghost+1, mitot-nghost
            do ivar=1,nvar
               if (dabs(alloc(iadd(ivar,i,j))) .lt. 1d-90) then
                  alloc(iadd(ivar,i,j)) = 0.d0
               endif
            enddo

            if (layers > 1) then
                do k=1,layers
                    index = 3*(k-1)
                    h(k) = alloc(iadd(index+1,i,j)) / rho(k)
                    hu(k) = alloc(iadd(index+2,i,j)) / rho(k)
                    hv(k) = alloc(iadd(index+3,i,j)) / rho(k)
                enddo
                eta(2) = h(2) + alloc(iaddaux(i,j,1))
            else
                k = 1
                index = 3*(k-1)
                h(k) = alloc(iadd(index+1,i,j))
                hu(k) = alloc(iadd(index+2,i,j))
                hv(k) = alloc(iadd(index+3,i,j))
                eta(2) = alloc(iaddaux(1,i,j))
            endif
            eta(1) = h(1) + eta(2)
            wind_x = alloc(iaddaux(4,i,j))
            if (abs(wind_x) < 1d-90) then
                wind_x = 0.d0
            endif
            wind_y = alloc(iaddaux(5,i,j))
            if (abs(wind_y) < 1d-90) then
                wind_y = 0.d0
            endif
            pressure = alloc(iaddaux(6,i,j))
            if (abs(pressure) < 1d-90) then
                pressure = 0.d0
            endif
C             pressure_x = alloc(iaddaux(i,j,7))
C             if (abs(pressure_x) < 1d-90) then
C                 pressure_x = 0.d0
C             endif
C             pressure_y = alloc(iaddaux(i,j,8))
C             if (abs(pressure_y) < 1d-90) then
C                 pressure_y = 0.d0
C             endif
C             vorticity = alloc(iaddaux(i,j,9))
C             if (abs(vorticity) < 1d-90) then
C                 vorticity = 0.d0
C             endif

              write(matunit1,109) (h(k),hu(k),hv(k), k=1,layers),
     &              (eta(k),k=1,layers),wind_x,wind_y,pressure
C               write(matunit1,109) h(1),hu(1),hv(1),h(2),hu(2),hv(2),
C      &            eta(1),eta(2),wind_x,wind_y
C             write(matunit1,109) (h(k),hu(k),hv(k),eta(k), k=1,layers),
C      &          wind_x,wind_y        
            
C             write(matunit1,109) (alloc(iadd(i,j,ivar)), ivar=1,nvar),
C      &         surface, 
C      &         wind_x,wind_y,
C      &         pressure,
C      &         pressure_x,
C      &         pressure_y,
C      &         vorticity
         enddo
         write(matunit1,*) ' '
      enddo
  109       format(4e26.16)


            mptr = node(levelptr, mptr)
            go to 70
 80      level = level + 1
         go to 65

 90     continue

        open(unit=matunit2,file=matname2,status='unknown',
     .       form='formatted')

c     # nvar+6 variable printed since surface, wind, pressure, pressure_x and pressure_y also printed

      write(matunit2,1000) time,4*layers+3,ngrids,naux,2
 1000 format(e18.8,'    time', /,
     &       i5,'                 meqn'/,
     &       i5,'                 ngrids'/,
     &       i5,'                 maux'/,
     &       i5,'                 ndim'/,/)
c

      write(6,601) matlabu,time
  601 format('AMRCLAW: Frame ',i4,
     &       ' matlab plot files done at time t = ', d12.5,/)

      matlabu = matlabu + 1

      close(unit=matunit1)
      close(unit=matunit2)
      endif

      return
      end
