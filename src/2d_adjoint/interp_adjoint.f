! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! :::::     Routine to interpolate adjoint to given x,y point
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      subroutine interp_adjoint(nvar, naux, x, y, q, k)

      use adjoint_module, only: adjoints
      implicit double precision (a-h,o-z)

      ! Function arguments
      integer, intent(in) :: naux, k
      integer :: ii_c, jj_c, ii_a, jj_a, nvar, z
      real(kind=8), intent(in) :: x,y
      integer :: nx,ny,loc,locaux,level,mptr,mitot,mjtot,
     .           ivar,iaux,i,j, iadd, iaddaux, iaddqeta
      real(kind=8) :: xlow, ylow, xhi, yhi, dx, dy,xm,ym,
     .                x_side, x_main, y_main, y_side
      real(kind=8) :: q(nvar), aux_a, aux_c,q_temp1(nvar),
     .                q_temp2(nvar), denom, aux_interp
      logical :: y_interp, yc_interp, ya_interp

      iadd(ivar,i,j)  = loc + ivar - 1 +
     .                     nvar*((j-1)*mitot+i-1)
      getaux(i,j) = adjoints(k)%alloc(iadd(4,i,j))
     .                     - adjoints(k)%alloc(iadd(1,i,j))

      do ivar=1,nvar
          q(ivar) = 0.0
      enddo

      do z = 1, adjoints(k)%ngrids
          mptr = adjoints(k)%gridpointer(z)
          level = adjoints(k)%gridlevel(mptr)

          ! Number of points in x and y (nx by ny grid)
          nx = adjoints(k)%ncellsx(mptr)
          ny = adjoints(k)%ncellsy(mptr)

          ! Finding x and y extreem values for grid
          xlow = adjoints(k)%xlowvals(mptr)
          ylow = adjoints(k)%ylowvals(mptr)
          dx = adjoints(k)%hxposs(level)
          dy = adjoints(k)%hyposs(level)
          xhi = xlow + nx*dx
          yhi = ylow + ny*dy

          loc     = adjoints(k)%loc(mptr)

          ! Total number of points in x and y
          mitot = nx + 2*adjoints(k)%nghost
          mjtot = ny + 2*adjoints(k)%nghost

          if ((x < xlow) .or. (x > xhi)
     .            .or. (y < ylow) .or. (y > yhi)) then
              ! Skipping interpolation if the point of interest
              ! is not in the current grid
              continue
          else
             xm = xlow - (adjoints(k)%nghost+0.5d0)*dx
             ym = ylow - (adjoints(k)%nghost+0.5d0)*dy

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
              aux_c = getaux(ii_c,jj_c)
              aux_a = getaux(ii_a,jj_c)

              if(sign(aux_c, aux_a) .ne. sign(aux_c, aux_c))then
                  yc_interp = .false.
              endif

              ! Verify that we should interpolate along adjacent y edge
              aux_c = getaux(ii_c,jj_a)
              aux_a = getaux(ii_a,jj_a)
              if(sign(aux_c, aux_a) .ne. sign(aux_c, aux_c))then
                  ya_interp = .false.
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
     &                     *adjoints(k)%alloc(iadd(ivar,ii_c,jj_c))
     &                     + ((y - y_main)/denom)
     &                     *adjoints(k)%alloc(iadd(ivar,ii_c,jj_a))
                    enddo
                  else
                      do ivar=1,nvar
                          q_temp1(ivar) =
     &                      adjoints(k)%alloc(iadd(ivar,ii_c,jj_c))
                      enddo
                  endif

                  if (ya_interp) then
                    do ivar=1,nvar
                        q_temp2(ivar) = ((y_side - y)/denom)
     &                     *adjoints(k)%alloc(iadd(ivar,ii_a,jj_c))
     &                     + ((y - y_main)/denom)
     &                     *adjoints(k)%alloc(iadd(ivar,ii_a,jj_a))
                    enddo
                  else
                    do ivar=1,nvar
                        q_temp2(ivar) =
     &                     adjoints(k)%alloc(iadd(ivar,ii_a,jj_c))
                    enddo
                  endif
              else
                  y_interp = .false.
              endif

              ! Verify that we should interpolate along x
              if (y_interp) then
                  aux_c = q_temp1(4) - q_temp1(1)
                  aux_a = q_temp2(4) - q_temp2(1)
              else
                  aux_c = getaux(ii_c,jj_c)
                  aux_a = getaux(ii_c,jj_a)
              endif
              if(sign(aux_c, aux_a) .ne. sign(aux_c, aux_c)) then
                  ii_a = ii_c
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
                  endif
              else
                  if (y_interp) then
                      q = q_temp1
                  else
                      do ivar=1,nvar
                          q(ivar) =
     &                      adjoints(k)%alloc(iadd(ivar,ii_c,jj_c))
                      enddo
                  endif
              endif
          endif

      enddo

      end subroutine interp_adjoint