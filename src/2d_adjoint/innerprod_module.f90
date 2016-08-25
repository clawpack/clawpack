! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! ::::::     Module to calculate the inner product with the adjoint
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module innerprod_module

contains

    function calculate_max_innerproduct(t,x_c,y_c,eta,q2,q3,aux1) result(max_innerprod)

        use adjoint_module

        real(kind=8), intent(in) :: t
        integer :: r
        real(kind=8) :: q_innerprod1, q_innerprod2, q_innerprod, max_innerprod
        real(kind=8) :: x_c,y_c,eta,q2,q3
        real(kind=8) :: aux1, t_nm

        max_innerprod = 0.d0

        aloop: do r=1, totnum_adjoints

            if (r .ne. 1) then
                t_nm = adjoints(r-1)%time
            else
                t_nm = 0.d0
            endif

            if (t < adjoints(r)%time .and. &
                (t +(trange_final - trange_start))>= t_nm) then

                q_innerprod1 = 0.d0
                q_innerprod2 = 0.d0

                q_innerprod1 = calculate_innerproduct(r,x_c,y_c,eta,q2,q3,aux1)
                if (r .ne. 1) then
                    q_innerprod2 = calculate_innerproduct(r-1,x_c,y_c,eta,q2,q3,aux1)
                endif

                q_innerprod = max(q_innerprod1, q_innerprod2)
                if (q_innerprod > max_innerprod) then
                    max_innerprod = q_innerprod
                endif
            endif

        enddo aloop

    end function calculate_max_innerproduct

    function calculate_innerproduct(r,x_c,y_c,eta,q2,q3,aux1) result(q_innerprod)

        use adjoint_module

        integer :: r
        real(kind=8) :: q_innerprod
        double precision, allocatable :: q_interp(:)
        real(kind=8) :: x_c,y_c,eta,q2,q3
        real(kind=8) :: aux_interp, aux1

        allocate(q_interp(nvar+1))

        ! If q is on land, don't computer the inner product
        if(aux1 > 0.d0) then
            q_innerprod = 0.d0
            return
        endif

        call interp_adjoint(1, adjoints(r)%lfine, nvar, &
            naux, x_c,y_c,q_interp, r)
        aux_interp = q_interp(4) - q_interp(1)

        ! If q_adjoint is on land, don't computer the inner product
        if(aux_interp > 0.d0) then
            q_innerprod = 0.d0
            return
        endif

        q_innerprod = abs( &
            eta * q_interp(4) &
            + q2 * q_interp(2) &
            + q3 * q_interp(3))

    end function calculate_innerproduct

end module innerprod_module
