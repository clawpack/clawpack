! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! :::::     Routine to interpolate adjoint to given x,y point
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      subroutine interp_adjoint(lst, lend, nvar, naux, x, y, q, k)

      use adjoint_module, only: adjoints
      use amr_reload_module
      implicit none

      ! Function arguments
      integer, intent(in) :: lst, lend, naux, k
      integer :: ii_c, jj_c, ii_a, jj_a, nvar
      real(kind=8), intent(in) :: x,y
      integer :: nx,ny,loc,locaux,level,mptr,mitot,mjtot,
     .           ivar,iaux,i,j, iadd, iaddaux, iaddqeta
      real(kind=8) :: xlow, ylow, xhi, yhi, dx, dy,xm,ym,
     .                x_side, x_main, y_main, y_side
      real(kind=8) :: q(nvar+1), aux_a, aux_c,q_temp1(nvar+1),
     .                q_temp2(nvar+1), denom, aux_interp
      logical :: y_interp, yc_interp, ya_interp

      iadd(ivar,i,j)  = loc + ivar - 1 + nvar*((j-1)*mitot+i-1)
      iaddaux(iaux,i,j) = locaux + iaux-1 + naux*(i-1) + 
     .                                      naux*mitot*(j-1)
      iaddqeta(ivar,i,j)  = 1 + ivar - 1 + 4*((j-1)*mitot+i-1)

      do ivar=1,nvar+1
          q(ivar) = 0.0
      enddo

         level = lst
 65      if (level .gt. lend) go to 90
            mptr = adjoints(k)%lstart(level)
 70         if (mptr .eq. 0) go to 80

              ! Number of points in x and y (nx by ny grid)
              nx = adjoints(k)%node(ndihi,mptr) -
     &                  adjoints(k)%node(ndilo,mptr) + 1
              ny = adjoints(k)%node(ndjhi,mptr) -
     &                  adjoints(k)%node(ndjlo,mptr) + 1

              ! Finding x and y extreem values for grid
              xlow = adjoints(k)%rnode(cornxlo,mptr)
              ylow = adjoints(k)%rnode(cornylo,mptr)
              dx = adjoints(k)%hxposs(level)
              dy = adjoints(k)%hyposs(level)
              xhi = xlow + nx*dx
              yhi = ylow + ny*dy

              loc     = adjoints(k)%node(store1, mptr)
              locaux  = adjoints(k)%node(storeaux,mptr)

              ! Total number of points in x and y
              mitot = nx + 2*nghost
              mjtot = ny + 2*nghost

              if ((x < xlow) .or. (x > xhi) .or. (y < ylow) .or.  
     &            (y > yhi)) then
                    ! Skipping interpolation if the point of interest
                    ! is not in the current grid
                    continue
              else
                  xm = xlow - (nghost+0.5d0)*dx
                  ym = ylow - (nghost+0.5d0)*dy

                  ! Finding current cell in x (i) and y (j)
                  ii_c = int((x-xm)/dx)
                  jj_c = int((y-ym)/dy)

                  ! Finding correct cell to interpolate with
                  jj_a = int(((y-ym)/dy) + 0.5d0)
                  ii_a = int(((x-xm)/dx) + 0.5d0)

                  if (jj_c == jj_a .and. jj_a /= 0) then
                     jj_a = jj_a - 1
                  endif
                  if (ii_c == ii_a .and. ii_a /= 0) then
                      ii_a = ii_a - 1
                  endif
                  if (jj_a >= ny) then
                      jj_a = jj_a - 1
                  endif
                  if (ii_a >= nx) then
                      ii_a = ii_a - 1
                  endif

                  ! If the cell we should interpolate with is not
                  ! in the same wet/dry state, don't interpolate

                  ! Verify that we should interpolate along center y line
                  aux_c = adjoints(k)%alloc(iaddaux(1,ii_c,jj_c))
                  aux_a = adjoints(k)%alloc(iaddaux(1,ii_c,jj_a))
                  if(sign(aux_c, aux_a) .ne. sign(aux_c, aux_c))then
                      yc_interp = .false.
                  endif

                  ! Verify that we should interpolate along adjacent y edge
                  aux_a = adjoints(k)%alloc(iaddaux(1,ii_a,jj_a))
                  if(sign(aux_c, aux_a) .ne. sign(aux_c, aux_c))then
                      ya_interp = .false.
                  endif

                  ! Verify that we should interpolate along x
                  aux_a = adjoints(k)%alloc(iaddaux(1,ii_a,jj_c))
                  if(sign(aux_c, aux_a) .ne. sign(aux_c, aux_c)) then
                      ii_a = ii_c
                  endif

                  ! Interpolating in y
                  y_main = ym + (jj_c + 0.5d0)*dy
                  if (jj_c /= jj_a) then
                      y_interp = .true.
                      y_side = ym + (jj_a + 0.5d0)*dy
                      denom = y_side - y_main

                      if (yc_interp) then
                        do ivar=1,nvar
                          q_temp1(ivar) = ((y_side - y)/denom)
     &                      *adjoints(k)%alloc(iadd(ivar,ii_c,jj_c))
     &                      + ((y - y_main)/denom)
     &                      *adjoints(k)%alloc(iadd(ivar,ii_c,jj_a))
                        enddo
                          ! Adding eta
                          aux_interp = ((y_side - y)/denom)
     &                        *adjoints(k)%alloc(iaddaux(1,ii_c,jj_c))
     &                        + ((y - y_main)/denom)
     &                        *adjoints(k)%alloc(iaddaux(1,ii_c,jj_a))
                          q_temp1(nvar+1) = q_temp1(1) + aux_interp
                      else
                        do ivar=1,nvar
                          q_temp1(ivar) =
     &                      adjoints(k)%alloc(iadd(ivar,ii_c,jj_c))
                        enddo
                          ! Adding eta
                          q_temp1(nvar + 1) = q_temp1(1) +
     &                        adjoints(k)%alloc(iaddaux(1,ii_c,jj_c))
                      endif

                      if (ya_interp) then
                        do ivar=1,nvar
                          q_temp2(ivar) = ((y_side - y)/denom)
     &                      *adjoints(k)%alloc(iadd(ivar,ii_a,jj_c))
     &                      + ((y - y_main)/denom)
     &                      *adjoints(k)%alloc(iadd(ivar,ii_a,jj_a))
                        enddo
                          ! Adding eta
                          aux_interp = ((y_side - y)/denom)
     &                        *adjoints(k)%alloc(iaddaux(1,ii_a,jj_c))
     &                        + ((y - y_main)/denom)
     &                        *adjoints(k)%alloc(iaddaux(1,ii_a,jj_a))
                          q_temp2(nvar+1) = q_temp2(1) + aux_interp
                      else
                        do ivar=1,nvar
                          q_temp2(ivar) =
     &                      adjoints(k)%alloc(iadd(ivar,ii_a,jj_c))
                        enddo
                          ! Adding eta
                          q_temp2(nvar + 1) = q_temp2(1) +
     &                        adjoints(k)%alloc(iaddaux(1,ii_a,jj_c))
                      endif
                  else
                      y_interp = .false.
                  endif


                  ! Interpolating in x
                  x_main = xm + (ii_c + 0.5d0)*dx
                  if (ii_c /= ii_a) then
                      x_side = xm + (ii_a + 0.5d0)*dx
                      denom = x_side - x_main

                      if(y_interp) then
                        q = ((x_side - x)/denom)*q_temp1
     &                    + ((x - x_main)/denom)*q_temp2
                      else
                        do ivar=1,nvar
                          q(ivar) = ((x_side - x)/denom)
     &                      *adjoints(k)%alloc(iadd(ivar,ii_c,jj_c))
     &                      + ((x - x_main)/denom)
     &                      *adjoints(k)%alloc(iadd(ivar,ii_a,jj_c))
                        enddo
                        ! Adding eta
                        aux_interp = ((x_side - x)/denom)
     &                        *adjoints(k)%alloc(iaddaux(1,ii_c,jj_c))
     &                        + ((x - x_main)/denom)
     &                        *adjoints(k)%alloc(iaddaux(1,ii_a,jj_c))
                        q(nvar+1) = q(1) + aux_interp
                      endif
                  else
                      if (y_interp) then
                        q = q_temp1
                      else
                        do ivar=1,nvar
                          q(ivar) =
     &                      adjoints(k)%alloc(iadd(ivar,ii_c,jj_c))
                        enddo
                        ! Adding eta
                        q(nvar + 1) = q(1) +
     &                      adjoints(k)%alloc(iaddaux(1,ii_c,jj_c))
                      endif
                  endif
              endif

           mptr = adjoints(k)%node(levelptr, mptr)
           go to 70
  80       level = level + 1
           go to 65
  90   continue

      end subroutine interp_adjoint