! Module containing refinement flagging criteria
module refinement_module

    use geoclaw_module, only: GEO_PARM_UNIT

    implicit none
    save

    ! ========================================================================
    !  Refinement Criteria
    ! ========================================================================
    real(kind=8), private, allocatable :: wave_tolerance(:)
    real(kind=8), private, allocatable :: speed_tolerance(:)
    real(kind=8), private :: deep_depth
    integer, private :: max_level_deep
    
    ! ========================================================================
    !  Refinement Regions
    ! ========================================================================
    type region_type
        integer :: min_level,max_level
        real(kind=8) :: x_low,y_low,x_hi,y_hi,t_low,t_hi
    end type region_type
    
    integer, private :: num_regions
    type(region_type), private, allocatable :: regions(:)
    
    ! ========================================================================
    !  Flowgrades - Not updated yet, use at your own risk
    ! ========================================================================
    integer, private :: num_flowgrades
    real(kind=8), private, allocatable :: flowgradevalue(:)
    integer, private, allocatable :: iflowgradevariable(:), iflowgradetype(:)
    integer, private, allocatable :: iflowgrademinlevel(:)

contains

    ! ========================================
    !  logical function allowflag(x,y,t)
    !
    !  Indicate whether the grid point at (x,y,t) at this refinement level
    !  is allowed to be flagged for further refinement.
    ! 
    !  Modified for GeoClaw to check whether the point lies in any of
    !  the various regions specified in the data files.
    ! 
    !  This routine is called from routine flag2refine.
    ! 
    !  If Richardson error estimates are used (if tol>0) then this routine
    !  is also called from errf1.
    ! 
    !  TODO: Does allowflag work for t0 < 0.0? Maybe not for iqinit > 0
    !
    logical function allowflag(x,y,t,level)

        use geoclaw_module
        use topo_module
        use dtopo_module
        use qinit_module
          
        implicit none

        ! Function arguments
        real(kind=8), intent(in) :: x,y,t
        integer, intent(in) :: level
      
        ! Locals
        integer :: m

        allowflag = .false.

    !   following commented by dlg on 10/9/08.
    !   my understanding of maxleveldeep might be differnet
    !   still shouldn't be allowed if maxlevel allowed in a region is less
    !   than maxleveldeep
    !   might want to allow high levels of refinement in some deep regions
    !   but not others.
    !
    !   if (level .lt. maxleveldeep) then
    !      # allow refinement to next level in deep water
    !       allowflag = .true.
    !       go to 900  !# no need to check anything else
    !       endif

        do m=1,mtopofiles
            if (level < maxleveltopo(m)) then
                if (x > xlowtopo(m) .and. x < xhitopo(m) .and. &
                    y > ylowtopo(m) .and. y < yhitopo(m) .and. &
                    t > tlowtopo(m) .and. t < thitopo(m)) then

                    allowflag = .true.
                    return
                endif
            endif
        enddo

        do m=1,num_regions
            if (level < regions(m)%max_level) then 
                if (x > regions(m)%x_low .and. x <  regions(m)%x_hi.and. &
                    y > regions(m)%y_low .and. y <  regions(m)%y_hi.and. &
                    t > regions(m)%t_low .and. t <= regions(m)%t_hi) then

                    allowflag = .true.
                    return
                endif
            endif
        enddo

        do m=1,num_dtopo
            if (x >  xlowdtopo(m) .and. x < xhidtopo(m).and. &
                y >  ylowdtopo(m) .and. y < yhidtopo(m).and. &
                t >= t0dtopo(m)   .and. t <= tfdtopo(m)) then
             
                if (level < maxleveldtopo(m)) then
                    allowflag = .true.
                    return
                endif
            endif
        enddo

        ! TODO: Correct the assumption that t0 = 0.0
        if (t == 0.d0 .and. qinit_type > 0) then
            if (x > x_low_qinit .and. x < x_hi_qinit .and. &
                y > y_low_qinit .and. y < y_hi_qinit) then

                if (level < max_level_qinit) then
                    allowflag = .true.
                    return
                endif
            endif
        endif

    end function allowflag
    

    ! ::::::::::::::::::::: flag2refine ::::::::::::::::::::::::::::::::::
    !
    ! User routine to control flagging of points for refinement.
    !
    ! Specific for GeoClaw for tsunami applications and related problems
    !
    !
    ! The logical function allowflag(x,y,t) is called to
    ! check whether further refinement at this level is allowed in this cell
    ! at this time.
    !
    !    q   = grid values including ghost cells (bndry vals at specified
    !          time have already been set, so can use ghost cell values too)
    !
    !  aux   = aux array on this grid patch
    !
    ! amrflags  = array to be flagged with either the value
    !             DONTFLAG (no refinement needed)  or
    !             DOFLAG   (refinement desired)
    !
    !
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    subroutine flag2refine(mx,my,mbc,meqn,maux,xlower,ylower,dx,dy,t,level,tolsp, &
                           q,aux,amrflags,DONTFLAG,DOFLAG)

        use amr_module, only: mxnest
        use geoclaw_module, only:dry_tolerance,rho,eta_init,num_layers
    
        use topo_module, only: tlowtopo,thitopo,xlowtopo,xhitopo,ylowtopo,yhitopo
        use topo_module, only: minleveltopo,mtopofiles
    
        use dtopo_module, only: tfdtopo,xlowdtopo,xhidtopo,ylowdtopo,yhidtopo
        use dtopo_module, only: minleveldtopo,num_dtopo
    
        use qinit_module, only: x_low_qinit,x_hi_qinit,y_low_qinit,y_hi_qinit
        use qinit_module, only: min_level_qinit,qinit_type
 
        implicit none
    
        ! Subroutine arguments
        integer, intent(in) :: mx,my,mbc,meqn,maux,level
        real(kind=8), intent(in) :: xlower,ylower,dx,dy,t,tolsp
    
        real(kind=8), intent(in) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
        real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    
        ! Flagging
        real(kind=8),intent(inout) :: amrflags(1-mbc:mx+mbc,1-mbc:my+mbc)
        real(kind=8), intent(in) :: DONTFLAG
        real(kind=8), intent(in) :: DOFLAG
    
        logical :: allowflag
        external allowflag
    
        ! Generic locals
        integer :: i,j,m,k,layer_index
        real(kind=8) :: x_c,y_c,x_low,y_low,x_hi,y_hi
        real(kind=8) :: h,speed,eta,eta_below

        ! Initialize flags
        amrflags = DONTFLAG

        ! Loop over interior points on this grid
        ! (i,j) grid cell is [x_low,x_hi] x [y_low,y_hi], cell center at (x_c,y_c)
        y_loop: do j=1,my
            y_low = ylower + (j - 1) * dy
            y_c = ylower + (j - 0.5d0) * dy
            y_hi = ylower + j * dy
        
            x_loop: do i = 1,mx
                x_low = xlower + (i - 1) * dx
                x_c = xlower + (i - 0.5d0) * dx
                x_hi = xlower + i * dx

                ! The following conditions are only checked in the horizontal and
                ! override the allowflag routine

                ! Check to see if refinement is forced in any topography file region:
                do m=1,mtopofiles
                    if (level < minleveltopo(m) .and. t >= tlowtopo(m) .and. t <= thitopo(m)) then
                        if (  x_hi > xlowtopo(m) .and. x_low < xhitopo(m) .and. &
                              y_hi > ylowtopo(m) .and. y_low < yhitopo(m) ) then
                        
                            amrflags(i,j) = DOFLAG
                            cycle x_loop
                        endif
                    endif
                enddo

                ! Check to see if refinement is forced in any other region:
                do m=1,num_regions
                    if (level < regions(m)%min_level .and. &
                        t >= regions(m)%t_low .and. t <= regions(m)%t_hi) then
                        if (x_hi > regions(m)%x_low .and. x_low < regions(m)%x_hi .and. &
                            y_hi > regions(m)%y_low .and. y_low < regions(m)%y_hi ) then
                    
                            amrflags(i,j) = DOFLAG
                            cycle x_loop
                        endif
                    endif
                enddo

                ! Check if we're in the dtopo region and need to refine:
                ! force refinement to level minleveldtopo
                do m = 1,num_dtopo
                    if (level < minleveldtopo(m).and. &
                        t <= tfdtopo(m) .and. & !t.ge.t0dtopo(m).and.
                        x_hi > xlowdtopo(m) .and. x_low < xhidtopo(m).and. &
                        y_hi > ylowdtopo(m) .and. y_low < yhidtopo(m)) then
                    
                        amrflags(i,j) = DOFLAG    
                        cycle x_loop
                    endif
                enddo

                ! Check if we're in the region where initial perturbation is
                ! specified and need to force refinement:
                ! This assumes that t0 = 0.d0, should really be t0 but we do
                ! not have access to that parameter in this routine
                if (qinit_type > 0 .and. t == 0.d0) then 
                    if (level < min_level_qinit .and. & 
                        x_hi > x_low_qinit .and. x_low < x_hi_qinit .and. &
                        y_hi > y_low_qinit .and. y_low < y_hi_qinit) then
                    
                        amrflags(i,j) = DOFLAG
                        cycle x_loop
                    endif
                endif

                ! -----------------------------------------------------------------
                ! Refinement not forced, so check if it is allowed and if so,
                ! check if there is a reason to flag this point:
                if (allowflag(x_c,y_c,t,level)) then
                    ! These refinement criteria are checked per layer going backwards
                    ! The bottom layer is checked first and eta_below is set to the 
                    ! bathymetry
                    eta_below = aux(1,i,j)
                
                    do k=num_layers,1,-1
                        layer_index = 3 * (k - 1)
                    
                        ! Extract state
                        h = q(layer_index+1,i,j) / rho(k)
                        if (h > dry_tolerance(k)) then
                            eta = h + eta_below
                    
                            ! Check wave criteria
                            if (abs(eta - eta_init(k)) > wave_tolerance(k)) then
                                ! Check to see if we are near shore
                                if (h < deep_depth) then
                                    amrflags(i,j) = DOFLAG
                                    cycle x_loop
                                ! If we are not in too deep of water, also flag if
                                ! we are allowed to
                                else if (level < max_level_deep) then
                                    amrflags(i,j) = DOFLAG
                                    cycle x_loop
                                endif
                            endif
                    
                            ! Check speed criteria, note that it might be useful to 
                            ! also have a per layer criteria since this is not 
                            ! gradient based
                            speed = sqrt(q(layer_index+2,i,j)**2 &
                                       + q(layer_index+3,i,j)**2) &
                                       / q(layer_index+1,i,j)
                            do m=1,mxnest
                                if (speed > speed_tolerance(m) .and. level <= m) then
                                    amrflags(i,j) = DOFLAG
                                    cycle x_loop
                                endif
                            enddo
                        endif
                    enddo
                endif
            
            enddo x_loop
        enddo y_loop
    end subroutine flag2refine
    
    ! =========================================================================
    !  Reads in the refinement control parameters
    ! =========================================================================
    subroutine set_refinement(file_name)
        
        use amr_module, only: mxnest
        use geoclaw_module, only: num_layers
        
        implicit none
        
        ! Arguments
        character(len=*), optional, intent(in) :: file_name
        
        ! Locals
        integer, parameter :: unit = 127
        integer :: i

        write(GEO_PARM_UNIT,*) ' '
        write(GEO_PARM_UNIT,*) '--------------------------------------------'
        write(GEO_PARM_UNIT,*) 'Refinement Control Parameters:'
        write(GEO_PARM_UNIT,*) '------------------------------'

        if (present(file_name)) then
            call opendatafile(unit, file_name)
        else
            call opendatafile(unit, 'refinement.data')
        endif

        ! Basic criteria
        allocate(wave_tolerance(num_layers))
        read(unit,*) wave_tolerance
        allocate(speed_tolerance(mxnest))
        read(unit,*) (speed_tolerance(i),i=1,mxnest)
        read(unit,*) deep_depth
        read(unit,*) max_level_deep
        read(unit,*)
        
        ! Refinement region data
        read(unit,"(i2)") num_regions
        allocate(regions(num_regions))
        do i=1,num_regions
            read(unit,*) regions(i)%min_level, regions(i)%max_level, &
                         regions(i)%t_low, regions(i)%t_hi, &
                         regions(i)%x_low, regions(i)%x_hi, &
                         regions(i)%y_low, regions(i)%y_hi
        enddo
        close(unit)
        
        ! Write out data to parameter file
        write(GEO_PARM_UNIT,*) '   wave_tolerance:',wave_tolerance
        write(GEO_PARM_UNIT,*) '   speed_tolerance:',speed_tolerance
        write(GEO_PARM_UNIT,*) '   maxleveldeep:', max_level_deep
        write(GEO_PARM_UNIT,*) '   depthdeep:', deep_depth
        write(GEO_PARM_UNIT,*) ''
        write(GEO_PARM_UNIT,*) '  num_regions = ',num_regions
        write(GEO_PARM_UNIT,*) '  minlevel, maxlevel, tlow, thi, xlow, xhi, ylow, yhigh values:'
        do i=1,num_regions
            write(GEO_PARM_UNIT,*) regions(i)%min_level, regions(i)%max_level, &
                                   regions(i)%t_low, regions(i)%t_hi, &
                                   regions(i)%x_low, regions(i)%x_hi, &
                                   regions(i)%y_low, regions(i)%y_hi
        enddo
        
    end subroutine set_refinement
    
    
    ! =========================================================================
    ! TODO: This needs to be updated for the new module
    ! =========================================================================
    subroutine set_flow_grades(file_name)

        implicit none

        ! Input arguments
        character(len=*), optional, intent(in) :: file_name

        ! Locals
        integer, parameter :: iunit = 127
        integer :: i

        write(GEO_PARM_UNIT,*) ' '
        write(GEO_PARM_UNIT,*) '--------------------------------------------'
        write(GEO_PARM_UNIT,*) 'SET FLOW GRADES:'
        write(GEO_PARM_UNIT,*) '------------'

        ! Read user parameters from setflowgrades.data
        if (present(file_name)) then
            call opendatafile(iunit, file_name)
        else
            call opendatafile(iunit, 'setflowgrades.data')
        endif
        
        read(iunit,*) num_flowgrades

        if (num_flowgrades == 0) then
            write(GEO_PARM_UNIT,*) '  No flow grades specified'
            return
        endif

        ! Allocate arrays
        allocate(flowgradevalue(num_flowgrades),iflowgradevariable(num_flowgrades))
        allocate(iflowgradetype(num_flowgrades),iflowgrademinlevel(num_flowgrades))

        do i=1,num_flowgrades
            read(iunit,*) flowgradevalue(i),iflowgradevariable(i), &
                iflowgradetype(i),iflowgrademinlevel(i)
        enddo

        close(iunit)

        write(GEO_PARM_UNIT,*) '   mflowgrades:',  num_flowgrades

        do i=1,num_flowgrades
            write(GEO_PARM_UNIT,"(d12.3,3i4)") flowgradevalue(i), &
                iflowgradevariable(i),iflowgradetype(i),iflowgrademinlevel(i)

        enddo

    end subroutine set_flow_grades
    
end module refinement_module