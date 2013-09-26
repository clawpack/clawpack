
subroutine fgmax_finalize()

    ! Print out the maxval and aux arrays and de-allocate storage.

    use fgmax_module
    use amr_module, only: mxnest

    implicit none
    character(30) :: fname
    character(1) :: cfg,cma
    integer :: k,ifg,level,mv,ma
    type(fgrid), pointer :: fg

    do ifg=1,FG_num_fgrids

        fg => FG_fgrids(ifg)   

        cfg = char(ichar('0') + ifg)
        fname = 'fort.FG' // cfg // '.valuemax'
        print *, 'Writing to file ', fname
        open(unit=FG_UNIT,file=trim(fname),status='unknown',form='formatted')

        do k=1,fg%npts
            do mv=1,FG_NUM_VAL
                if (abs(fg%valuemax(mv,k)) .lt. 1.d-90) then
                    fg%valuemax(mv,k) = 0.d0
                    endif
                enddo
            write(FG_UNIT,111) fg%x(k),fg%y(k), fg%levelmax(k), &
                  (fg%valuemax(mv,k), mv=1,FG_NUM_VAL), &
                  (fg%tmax(mv,k), mv=1,FG_NUM_VAL), fg%arrival_time(k)
 111        format(2e17.8,i4,21e17.8)
            enddo

        close(FG_UNIT)


        do ma=1,FG_NUM_AUX
            cma = char(ichar('0') + ma)
            fname = 'fort.FG' // cfg // '.aux' // cma
            print *, 'Writing to file ', fname
            open(unit=FG_UNIT,file=trim(fname),status='unknown',form='formatted')

            do k=1,fg%npts
                write(FG_UNIT,112) fg%x(k),fg%y(k), &
                      (fg%aux(level,ma,k), level=1,mxnest)
 112            format(2e17.8,20e17.8)
                enddo

            close(FG_UNIT)
            enddo

        ! deallocate(fg%valuemax,fg%levelmax,fg%aux,fg%x,fg%y)
        enddo
    if (FG_DEBUG) then
        write(6,*) '+++ Fixed grid debugging written to fort.61 and fort.65'   
        endif

end subroutine fgmax_finalize
