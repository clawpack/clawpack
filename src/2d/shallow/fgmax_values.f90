
subroutine fgmax_values(mx,my,meqn,mbc,maux,q,aux,dx,dy, &
                   xlower,ylower,mask_patch,values)

    ! Given a grid q (and aux if needed), set the elements of 
    !   values(mv,i,j)  for mv=1:FG_NUM_VAL
    ! to the desired values that will be output and/or monitored on
    ! the fixed grid(s).
    !
    ! Only the elements for which mask_patch(i,j) == .true. need be set.

    use fgmax_module
    use geoclaw_module, only: sea_level, dry_tolerance

    implicit none
    integer, intent(in) :: mx,my,meqn,mbc,maux
    real(kind=8), intent(in) :: q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(in) :: aux(maux, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(in) :: dx,dy,xlower,ylower
    logical, intent(in) :: mask_patch(1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: values(FG_NUM_VAL, 1-mbc:mx+mbc, 1-mbc:my+mbc)

    ! This version sets only 1 value to monitor, the value eta_tilde
    ! defined to be h+B where wet and sea_level where dry

    if (FG_NUM_VAL .ne. 1) then
        write(6,*) '*** Error FG_NUM_VAL in fgmax_module is ',FG_NUM_VAL
        write(6,*) '*** Does not agree with number expected in fgmax_values: 1'
        stop
        endif 

    where (mask_patch .and. (q(1,:,:) >= dry_tolerance))
        values(1,:,:) = q(1,:,:) + aux(1,:,:)
    endwhere

    where (mask_patch .and. (q(1,:,:) < dry_tolerance))
        values(1,:,:) = -1.d90   !# to indicate dry region
    endwhere

end subroutine fgmax_values
