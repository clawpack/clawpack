! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! :::::     Routine to interpolate adjoint to given x,y point
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

 subroutine interp_adjoint(nvar, r, q_interp, xlower_a, ylower_a, &
              dx_a, dy_a, mx_a, my_a, xlower_f, ylower_f, dx_f, &
              dy_f, mx_f, my_f, mask_adjoint, mptr_a, mask_forward)

      use adjoint_module, only: adjoints
      implicit none

      ! Function arguments
      integer, intent(in) :: nvar, r, mx_a, my_a
      logical, intent(in) :: mask_adjoint( &
                        1-adjoints(r)%nghost:mx_a+adjoints(r)%nghost, &
                        1-adjoints(r)%nghost:my_a+adjoints(r)%nghost)
      real(kind=8), intent(in) :: xlower_a, xlower_f, ylower_a, ylower_f
      integer, intent(in) :: mx_f, my_f, mptr_a
      real(kind=8), intent(in) :: dx_f, dx_a, dy_f, dy_a

      integer :: z,k, iz, jk, mitot
      integer :: ivar, i, j, iadd, loc
      real(kind=8) :: q_interp(nvar,mx_f,my_f), denom
      real(kind=8) :: x, xhigh_a, y, yhigh_a
      real(kind=8) :: dxz,dxzp,dyk,dykp, a, b, getaux
      logical :: mask_forward(mx_f,my_f)

      iadd(ivar,i,j)  = loc + ivar - 1 + nvar*((j-1)*mitot+i-1)
      getaux(i,j) = adjoints(r)%alloc(iadd(4,i,j)) &
                          - adjoints(r)%alloc(iadd(1,i,j))

      q_interp = 0.0
      xhigh_a  = xlower_a + mx_a*dx_a
      yhigh_a = ylower_a + my_a*dx_a
      loc    = adjoints(r)%loc(mptr_a)
      mitot = adjoints(r)%ncellsx(mptr_a) + 2*adjoints(r)%nghost

      do z = 1, mx_f
        do k = 1,my_f
          if (mask_forward(z,k)) then
            x = xlower_f + (z - 0.5d0)*dx_f
            y = ylower_f + (k - 0.5d0)*dy_f

            iz = int((x - xlower_a + 0.5d0*dx_a) / dx_a) + 1
            dxz = x - (xlower_a + (iz-0.5d0)*dx_a)
            dxzp = (xlower_a + ((iz+1)-0.5d0)*dx_a) - x
            jk = int((y - ylower_a + 0.5d0*dy_a) / dy_a) + 1
            dyk = y - (ylower_a + (jk-0.5d0)*dy_a)
            dykp = (ylower_a + ((jk+1)-0.5d0)*dy_a) - y

            ! Interpolate only if this cell is overlapping with grid
            if (mask_adjoint(iz,jk)) then
              do ivar=1,nvar

                ! Interpolate only if cells are in same wet/dry state
                if(getaux(iz+1,jk) <= 0.d0 .and. getaux(iz,jk) <= 0.d0) then
                  a = (dxzp*adjoints(r)%alloc(iadd(ivar,iz,jk)) + &
                      dxz*adjoints(r)%alloc(iadd(ivar,iz+1,jk)))/dx_a
                else
                  a = adjoints(r)%alloc(iadd(ivar,iz,jk))
                endif

                if(getaux(iz+1,jk+1) <= 0.d0 .and. getaux(iz,jk+1) <= 0.d0) then
                  b = (dxzp*adjoints(r)%alloc(iadd(ivar,iz,jk+1)) + &
                      dxz*adjoints(r)%alloc(iadd(ivar,iz+1,jk+1)))/dx_a
                else
                  b = adjoints(r)%alloc(iadd(ivar,iz,jk+1))
                endif

                if(getaux(iz,jk) <= 0.d0 .and. getaux(iz,jk+1) <= 0.d0) then
                  q_interp(ivar,z,k) = (dykp*a + dyk*b)/dy_a
                else
                  q_interp(ivar,z,k) = a
                endif

              enddo
            endif
          endif
        enddo
      enddo

      end subroutine interp_adjoint