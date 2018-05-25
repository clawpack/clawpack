! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! ::::::     Module to calculate the inner product with the adjoint
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module innerprod_module

contains

    subroutine calculate_innerproduct(t,q,k,mx_f,my_f,xlower_f, &
               ylower_f,dx_f,dy_f,meqn_f,mbc_f,aux1,innerprod)

        use adjoint_module
        use amr_module, only : flag_richardson, flag_gradient

        implicit none

        real(kind=8), intent(in) :: t,xlower_f,ylower_f,dx_f,dy_f
        integer :: k,mx_f,my_f,meqn_f,mbc_f
        real(kind=8), intent(in) :: q(meqn_f,1-mbc_f:mx_f+mbc_f,1-mbc_f:my_f+mbc_f)

        integer :: mx_a, my_a, mptr_a, mbc_a
        integer :: i, j, i1, i2, j1, j2, level, loc, z, m
        real(kind=8) :: dx_a, xlower_a, xupper_a, xupper_f
        real(kind=8) :: dy_a, ylower_a, yupper_a, yupper_f
        real(kind=8) :: x1, x2, y1, y2

        real(kind=8), intent(inout) :: innerprod(mx_f,my_f)
        real(kind=8) :: q_innerprod(mx_f,my_f)
        logical :: mask_forward(mx_f,my_f)
        real(kind=8) :: q_interp(adjoints(k)%meqn,mx_f,my_f), eta
        real(kind=8) :: aux1(1-mbc_f:mx_f+mbc_f,1-mbc_f:my_f+mbc_f)

        logical, allocatable :: mask_adjoint(:,:)

        xupper_f = xlower_f + mx_f*dx_f
        yupper_f = ylower_f + my_f*dy_f

        ! Loop over patches in adjoint solution
        do z = 1, adjoints(k)%ngrids
            mptr_a = adjoints(k)%gridpointer(z)
            level = adjoints(k)%gridlevel(mptr_a)

            ! Number of points in x and y
            mx_a = adjoints(k)%ncellsx(mptr_a)
            my_a = adjoints(k)%ncellsy(mptr_a)

            ! Finding extreem values for grid
            xlower_a = adjoints(k)%xlowvals(mptr_a)
            dx_a = adjoints(k)%hxposs(level)
            xupper_a = xlower_a + mx_a*dx_a
            ylower_a = adjoints(k)%ylowvals(mptr_a)
            dy_a = adjoints(k)%hyposs(level)
            yupper_a = ylower_a + my_a*dy_a

            loc = adjoints(k)%loc(mptr_a)

            ! Check if adjoint patch overlaps with forward patch
            x1 = max(xlower_f,xlower_a)
            x2 = min(xupper_f,xupper_a)
            y1 = max(ylower_f,ylower_a)
            y2 = min(yupper_f,yupper_a)

            if ((x1 > x2) .or. (y1 > y2)) then
                ! Skipping interpolation if grids don't overlap
                mask_forward = .false.
                continue
            else
                mbc_a = adjoints(k)%nghost
                allocate(mask_adjoint(1-mbc_a:mx_a+mbc_a, 1-mbc_a:my_a+mbc_a))

                ! Create a mask that is .true. only in part of patch intersecting forward patch:
                i1 = max(int((x1 - xlower_a + 0.5d0*dx_a) / dx_a), 0)
                i2 = min(int((x2 - xlower_a + 0.5d0*dx_a) / dx_a) + 1, mx_a+1)
                j1 = max(int((y1 - ylower_a + 0.5d0*dy_a) / dy_a), 0)
                j2 = min(int((y2 - ylower_a + 0.5d0*dy_a) / dy_a) + 1, my_a+1)

                forall (i=1-mbc_a:mx_a+mbc_a, j=1-mbc_a:my_a+mbc_a)
                    mask_adjoint(i,j) = ((i >= i1) .and. (i <= i2) .and. &
                            (j >= j1) .and. (j <= j2))
                end forall

                ! Create a mask that is .true. only in part of forward patch intersecting patch:

                i1 = max(int((x1 - xlower_f + 0.5d0*dx_f) / dx_f)+1, 0)
                i2 = min(int((x2 - xlower_f + 0.5d0*dx_f) / dx_f), mx_f)
                j1 = max(int((y1 - ylower_f + 0.5d0*dy_f) / dy_f)+1, 0)
                j2 = min(int((y2 - ylower_f + 0.5d0*dy_f) / dy_f), my_f)

                forall (i=1:mx_f, j=1:my_f)
                    mask_forward(i,j) = ((i >= i1) .and. (i <= i2) .and. &
                             (j >= j1) .and. (j <= j2))
                end forall

                ! Interpolate adjoint values to q_interp
                ! Note that values in q_interp will only be set properly where
                ! mask_adjoint == .true.
                call interp_adjoint( &
                    adjoints(k)%meqn, k, q_interp, xlower_a, ylower_a, &
                    dx_a, dy_a, mx_a, my_a, xlower_f, ylower_f, &
                    dx_f, dy_f, mx_f, my_f, &
                    mask_adjoint, mptr_a, mask_forward)

                q_innerprod = 0.d0

                ! Checking wet and dry states
                forall(i = 1:mx_f, j = 1:my_f, mask_forward(i,j))
                    mask_forward(i,j) = (mask_forward(i,j) .and. &
                          (aux1(i,j) < 0.d0) .and. &
                          (q_interp(4,i,j) - q_interp(1,i,j) < 0.d0))
                end forall

                ! For each valid point, calculate inner product
                if (flag_gradient) then
                  forall(i = 1:mx_f, j = 1:my_f, mask_forward(i,j))
                      q_innerprod(i,j) = abs( &
                        (q(1,i,j) + aux1(i,j)) * q_interp(4,i,j) &
                        + q(2,i,j) * q_interp(2,i,j) &
                        + q(3,i,j) * q_interp(3,i,j))
                  end forall
                endif

                if (flag_richardson) then
                  forall(i = 1:mx_f, j = 1:my_f, mask_forward(i,j))
                    q_innerprod(i,j) = abs( &
                        q(1,i,j) * q_interp(4,i,j) &
                        + q(2,i,j) * q_interp(2,i,j) &
                        + q(3,i,j) * q_interp(3,i,j))
                  end forall
                endif


                do i=1,mx_f
                    do j=1,my_f
                        if (q_innerprod(i,j) > innerprod(i,j)) then
                            innerprod(i,j) = q_innerprod(i,j)
                        endif
                    enddo
                enddo

                deallocate(mask_adjoint)
            endif
        enddo


    end subroutine calculate_innerproduct

end module innerprod_module
