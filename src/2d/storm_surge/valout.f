c
c -----------------------------------------------------
c
      subroutine valout (lst, lend, time, nvar, naux)
c
      use amr_module
      use geoclaw_module, only: num_layers,rho
      use storm_module, only: wind_index, pressure_index
      use storm_module, only: output_storm_location
      implicit double precision (a-h,o-z)
      character*10  matname1, matname2

c     Work arrays
      dimension eta(num_layers),h(num_layers)
      dimension hu(num_layers),hv(num_layers)
      dimension storm_field(3)


c OLD INDEXING
c     iadd(i,j,ivar) = loc + i - 1 + mitot*((ivar-1)*mjtot+j-1)
c     iaddaux(i,j,iaux) = locaux + i - 1 + mitot*((iaux-1)*mjtot+j-1)
c NEW INDEXING ORDER SWITCHED
      iadd(ivar,i,j)  = loc + ivar - 1 + nvar*((j-1)*mitot+i-1)
      iaddaux(iaux,i,j) = locaux + iaux-1 + naux*(i-1) +
     .                                      naux*mitot*(j-1)



c ::::::::::::::::::::::::::: VALOUT ::::::::::::::::::::::::::::::::::;
c graphics output of soln values for contour or surface plots.
c modified for GeoClaw to output the surface level along with q.
c    surface = q(i,j,1) + aux(i,j,1)
c :::::::::::::::::::::::::::::::::::::;:::::::::::::::::::::::::::::::;

c
c     ### MATLAB graphics output
c
      if (matlabout) then
         ! Output the storm location at this time
         call output_storm_location(time)

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
               if (abs(alloc(iadd(ivar,i,j))) < 1d-90) then
                  alloc(iadd(ivar,i,j)) = 0.d0
               endif
            enddo
            
            ! Extract all but bottom layer depth and momenta
            do k=1,num_layers-1
                index = 3 * (k - 1)
                h(k) = alloc(iadd(index+1,i,j)) / rho(k)
                hu(k) = alloc(iadd(index+2,i,j)) / rho(k)
                hv(k) = alloc(iadd(index+3,i,j)) / rho(k)
            enddo
            index = 3 * (num_layers - 1)
            h(num_layers) = alloc(iadd(index+1,i,j)) / rho(num_layers)
            hu(num_layers) = alloc(iadd(index+2,i,j)) / rho(num_layers)
            hv(num_layers) = alloc(iadd(index+3,i,j)) / rho(num_layers)
            
            ! Calculate surfaces
            eta(num_layers) = h(num_layers) + alloc(iaddaux(1,i,j))
            do k=num_layers-1,1,-1
                eta(k) = h(k) + eta(k+1)
            enddo

            ! Storm fields and location
            storm_field(1) = alloc(iaddaux(wind_index,i,j))
            storm_field(2) = alloc(iaddaux(wind_index+1,i,j))
            storm_field(3) = alloc(iaddaux(pressure_index,i,j))
            forall (k=1:3,storm_field(k) < 1d-90)
              storm_field(k) = 0.d0
            end forall
            
            write(matunit1,109) (h(k),hu(k),hv(k), k=1,num_layers),
     &                          (eta(k),k=1,num_layers),
     &                          (storm_field(k),k=1,3)
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

c     # nvar+1 variable printed since surface also printed

      write(matunit2,1000) time,4*num_layers+3,ngrids,naux,2
 1000 format(e18.8,'    time', /,
     &       i5,'                 meqn'/,
     &       i5,'                 ngrids'/,
     &       i5,'                 naux'/,
     &       i5,'                 ndim'/,/)
c

      write(6,601) matlabu,time
  601 format('GeoClaw: Frame ',i4,
     &       ' output files done at time t = ', d12.6,/)

      matlabu = matlabu + 1

      close(unit=matunit1)
      close(unit=matunit2)
      endif

      return
      end
