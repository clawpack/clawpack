
subroutine fixedgrid_values(mx,my,meqn,mbc,maux,q,aux,dx,dy, &
                   xlower,ylower,mask_patch,values)

    ! Given a grid q (and aux if needed), set the elements of 
    !   values(mv,i,j)  for mv=1:FG_NUM_MAXVALS
    ! to the desired values that will be output and/or monitored on
    ! the fixed grid(s).
    !
    ! Only the elements for which mask_patch(i,j) == .true. need be set.

    use fixedgrid_module

    implicit none
    integer, intent(in) :: mx,my,meqn,mbc,maux
    real(kind=8), intent(in) :: q(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn)
    real(kind=8), intent(in) :: aux(1-mbc:mx+mbc, 1-mbc:my+mbc, maux)
    real(kind=8), intent(in) :: dx,dy,xlower,ylower
    logical, intent(in) :: mask_patch(1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: values(FG_NUM_VAL, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    integer :: i,j

    real(kind=8) :: h,s,hs,hss,s_dry_tol

    s_dry_tol = 1.d-2

    !forall (i=1-mbc:mx+mbc, j=1-mbc:my+mbc, mask_patch(i,j))
    !   values(1,i,j) = q(i,j,1)
    !   values(2,i,j) = q(i,j,1) + aux(i,j,1)
    !   end forall

    if (FG_NUM_VAL .ne. 4) then
        write(6,*) '*** Error FG_NUM_VAL in fixedgrid_module.f90 is ',FG_NUM_VAL
        write(6,*) '*** Does not agree with number of values set in fixedgrid_values.f90'
        stop
        endif 

    do i=1-mbc,mx+mbc
        do j=1-mbc,my+mbc
            if (mask_patch(i,j)) then
                h = q(i,j,1)  ! depth
                hs = sqrt(q(i,j,2)**2 + q(i,j,3)**2) 
                if (h > s_dry_tol) then
                    s = hs / h
                  else
                    s = 0.d0
                  endif
                hs = h*s
                hss = h*hs
                values(1,i,j) = h
                values(2,i,j) = s
                values(3,i,j) = hs
                values(4,i,j) = hss
                endif
            enddo
        enddo

end subroutine fixedgrid_values
