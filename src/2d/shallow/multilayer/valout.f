c
c -----------------------------------------------------
c
      subroutine valout (lst, lend, time, nvar, naux)
c
      use amr_module
      use multilayer_module, only: num_layers, rho
      implicit double precision (a-h,o-z)
      character(len=10) :: fname1, fname2, fname3, fname4

c     # GeoClaw specific output....  add eta to q array before printing
c
c     # Write the results to the file fort.q<iframe>
c     # Use format required by matlab script  plotclaw2.m or Python tools
c
c     # set outaux = .true. to also output the aux arrays to fort.a<iframe>

      logical outaux
      integer output_aux_num 
      real(kind=8) :: h(num_layers), hu(num_layers), hv(num_layers)
      real(kind=8) :: eta(num_layers)
      real(kind=8), allocatable :: qeta(:)

      iadd(ivar,i,j)  = loc + ivar - 1 + nvar*((j-1)*mitot+i-1)
      iaddaux(iaux,i,j) = locaux + iaux-1 + naux*(i-1) +
     .                                      naux*mitot*(j-1)
      iaddqeta(ivar,i,j)  = 1 + ivar - 1+4*num_layers*((j-1)*mitot+i-1) 
c

c     # how many aux components requested?
      output_aux_num = 0
      do i=1,naux
         output_aux_num = output_aux_num + output_aux_components(i)
         enddo
        
c     # Currently outputs all aux components if any are requested!
      outaux = ((output_aux_num > 0) .and. 
     .         ((.not. output_aux_onlyonce) .or. 
     .          (abs(time - t0) < 1d-8)))

c     open(unit=77,file='fort.b',status='unknown',access='stream')


c     ### Python graphics output
c

c        ###  make the file names and open output files
         fname1 = 'fort.qxxxx'
         fname2 = 'fort.txxxx'
         fname3 = 'fort.axxxx'
         fname4 = 'fort.bxxxx'
         matunit1 = 50
         matunit2 = 60
         matunit3 = 70
         matunit4 = 71
         nstp     = matlabu
         do 55 ipos = 10, 7, -1
            idigit = mod(nstp,10)
            fname1(ipos:ipos) = char(ichar('0') + idigit)
            fname2(ipos:ipos) = char(ichar('0') + idigit)
            fname3(ipos:ipos) = char(ichar('0') + idigit)
            fname4(ipos:ipos) = char(ichar('0') + idigit)
            nstp = nstp / 10
 55      continue

         open(unit=matunit1,file=fname1,status='unknown',
     .       form='formatted')

         if (output_format == 3) then
c            # binary output          
             open(unit=matunit4,file=fname4,status='unknown',
     &               access='stream')
             endif

         level = lst
         ngrids = 0
c65      if (level .gt. lfine) go to 90
 65      if (level .gt. lend) go to 90
            mptr = lstart(level)
 70         if (mptr .eq. 0) go to 80
              ngrids  = ngrids + 1
              nx      = node(ndihi,mptr) - node(ndilo,mptr) + 1
              ny      = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
              loc     = node(store1, mptr)
              locaux  = node(storeaux,mptr)
              mitot   = nx + 2*nghost
              mjtot   = ny + 2*nghost
              if (ny.gt.1) then
                  write(matunit1,1001) mptr, level, nx, ny
                else
c                 # output in 1d format if ny=1:
                  write(matunit1,1003) mptr, level, nx
                endif
 1001 format(i5,'                 grid_number',/,
     &       i5,'                 AMR_level',/,
     &       i5,'                 mx',/,
     &       i5,'                 my')
 1003 format(i5,'                 grid_number',/,
     &       i5,'                 AMR_level',/,
     &       i5,'                 mx')


              xlow = rnode(cornxlo,mptr)
              ylow = rnode(cornylo,mptr)
              if (ny.gt.1) then
                  write(matunit1,1002)
     &              xlow,ylow,hxposs(level),hyposs(level)
                else
                  write(matunit1,1004)
     &              xlow,hxposs(level)
                endif
 1002 format(e18.8,'    xlow', /,
     &       e18.8,'    ylow', /,
     &       e18.8,'    dx', /,
     &       e18.8,'    dy',/)
 1004 format(e18.8,'    xlow', /,
     &       e18.8,'    dx', /)


        if (output_format == 1) then
             do j = nghost+1, mjtot-nghost
                do i = nghost+1, mitot-nghost
                   do ivar=1,nvar
                      if (abs(alloc(iadd(ivar,i,j))) < 1d-90) then
                         alloc(iadd(ivar,i,j)) = 0.d0
                      endif
                   enddo

                   ! Extract depth and momenta
                   do k=1,num_layers-1
                     index = 3 * (k - 1)
                     h(k) = alloc(iadd(index + 1,i,j)) / rho(k)
                     hu(k)= alloc(iadd(index + 2,i,j)) / rho(k)
                     hv(k) = alloc(iadd(index + 3,i,j)) / rho(k)
                   end do
                  ! Extract bottom layer
                  index = 3 * (num_layers - 1) 
                  h(num_layers) = alloc(iadd(index+1,i,j)) 
     &                                                 / rho(num_layers)
                  hu(num_layers) = alloc(iadd(index+2,i,j)) 
     &                                                 / rho(num_layers)
                  hv(num_layers) = alloc(iadd(index+3,i,j)) 
     &                                                 / rho(num_layers)

                  ! Calculate surfaces
                  eta(num_layers) = h(num_layers) 
     &                                           + alloc(iaddaux(1,i,j))
                  do k=num_layers-1,1,-1
                    eta(k) = h(k) + eta(k+1)
                    if (abs(eta(k)) < 1d-90) then
                      eta(k) = 0.d0
                    endif
                  enddo

                  write(matunit1,109) (h(k),hu(k),hv(k),k=1,num_layers),
     &                                (eta(k),k=1,num_layers)
                enddo
                write(matunit1,*) ' '
             enddo
  109        format(4e26.16)
         endif

         if (output_format == 3) then
C            # binary output          
           i1 = iadd(1,1,1)
           i2 = iadd(nvar,mitot,mjtot)
C            # Need to augment q with eta:
             allocate(qeta(4*num_layers*mitot*mjtot))
             do j=1,mjtot
                 do i=1,mitot

                   ! Extract depth and momenta
                   do k=1,num_layers-1
                     index = 3 * (k - 1)
                     h(k) = alloc(iadd(index + 1,i,j)) / rho(k)
                     hu(k)= alloc(iadd(index + 2,i,j)) / rho(k)
                     hv(k) = alloc(iadd(index + 3,i,j)) / rho(k)
                   end do                  
                   ! Extract bottom layer
                  index = 3 * (num_layers - 1) 
                  h(num_layers) = alloc(iadd(index+1,i,j)) 
     &                                                 / rho(num_layers)
                  hu(num_layers) = alloc(iadd(index+2,i,j)) 
     &                                                 / rho(num_layers)
                  hv(num_layers) = alloc(iadd(index+3,i,j)) 
     &                                                 / rho(num_layers)

                  ! Calculate surfaces
                  eta(num_layers) = h(num_layers) 
     &                                           + alloc(iaddaux(1,i,j))
                  do k=num_layers-1,1,-1
                    eta(k) = h(k) + eta(k+1)
                  enddo

                  do m=1,num_layers
                    index = 3*(m - 1)
                    qeta(iaddqeta(index+1,i,j)) = h(m)
                    qeta(iaddqeta(index+2,i,j)) = hu(m)
                    qeta(iaddqeta(index+3,i,j)) = hv(m)
                    qeta(iaddqeta(3*num_layers+m,i,j)) = eta(m)
C                     print *, iaddqeta(index+1,i,j)
C                     print *, iaddqeta(3*num_layers+m,i,j)
C                     print *, h(m),hu(m),hv(m),eta(m)
                  end do
C                  stop
                 enddo
             enddo

c            # NOTE: we are writing out ghost cell data also, unlike ascii
             write(matunit4) qeta

             deallocate(qeta)
             endif

            mptr = node(levelptr, mptr)
            go to 70
 80      level = level + 1
         go to 65

 90     continue

c       -------------------
c       # output aux arrays
c       -------------------

        if (outaux) then
c        # output aux array to fort.aXXXX

         level = lst
 165     if (level .gt. lend) go to 190
            mptr = lstart(level)
 170        if (mptr .eq. 0) go to 180
              nx      = node(ndihi,mptr) - node(ndilo,mptr) + 1
              ny      = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
              locaux  = node(storeaux,mptr)
              mitot   = nx + 2*nghost
              mjtot   = ny + 2*nghost


          if (output_format == 1) then
             open(unit=matunit3,file=fname3,status='unknown',
     .            form='formatted')
              if (ny.gt.1) then
                  write(matunit3,1001) mptr, level, nx, ny
                else
c                 # output in 1d format if ny=1:
                  write(matunit3,1003) mptr, level, nx
                endif
              xlow = rnode(cornxlo,mptr)
              ylow = rnode(cornylo,mptr)
              if (ny.gt.1) then
                  write(matunit3,1002)
     &              xlow,ylow,hxposs(level),hyposs(level)
                else
                  write(matunit3,1004)
     &              xlow,hxposs(level)
                endif

             do j = nghost+1, mjtot-nghost
                do i = nghost+1, mitot-nghost
                   do ivar=1,naux
                      if (abs(alloc(iaddaux(ivar,i,j))) .lt. 1d-90) 
     &                   alloc(iaddaux(ivar,i,j)) = 0.d0
                   enddo
                   write(matunit3,109) (alloc(iaddaux(ivar,i,j)), 
     &                              ivar=1,naux)
                enddo
                write(matunit3,*) ' '
             enddo
            endif
            
         if (output_format == 3) then
c            # binary output          
             open(unit=matunit3,file=fname3,status='unknown',
     &               access='stream')
             i1 = iaddaux(1,1,1)
             i2 = iaddaux(naux,mitot,mjtot)
c            # NOTE: we are writing out ghost cell data also, unlike ascii
             write(matunit3) alloc(i1:i2)
             endif


            mptr = node(levelptr, mptr)
            go to 170
 180     level = level + 1
         go to 165

 190    continue
        close(unit=matunit3)
        endif !# end outputting aux array


c     --------------
c     # fort.t file:
c     --------------

      open(unit=matunit2,file=fname2,status='unknown',
     .       form='formatted')
      if (ny.gt.1) then 
          ndim = 2
        else
c         # special case where 2d AMR is used for a 1d problem
c         # and we want to use 1d plotting routines
          ndim = 1
        endif

c     # NOTE: we need to print out nghost too in order to strip
c     #       ghost cells from q when reading in pyclaw.io.binary
c     # Print meqn = nvar+1 because eta added.
      write(matunit2,1000) time,nvar+num_layers,ngrids,naux,ndim,nghost
 1000 format(e18.8,'    time', /,
     &       i5,'                 meqn'/,
     &       i5,'                 ngrids'/,
     &       i5,'                 naux'/,
     &       i5,'                 ndim'/,
     &       i5,'                 nghost'/,/)
c

      write(6,601) matlabu,time
  601 format('AMRCLAW: Frame ',i4,
     &       ' output files done at time t = ', d12.6,/)

      matlabu = matlabu + 1

      close(unit=matunit1)
      close(unit=matunit2)
      if (output_format == 3) then
          close(unit=matunit4)
          endif

      return
      end