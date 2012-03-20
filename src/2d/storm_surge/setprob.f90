subroutine setprob()

    use hurricane_module
    use multilayer_module
    use geoclaw_module
    use topo_module
    use dtopo_module

    implicit none
    
    call set_geo           ! sets basic parameters g and coord system
    call set_shallow       ! sets parameters specific to tsunamis
    call set_topo          ! specifies topography (bathymetry) files
    call set_dtopo         ! specifies file with dtopo from earthquake
    call setqinit         ! specifies file with dh if this used instead
    call setregions       ! specifies where refinement is allowed/forced
    call setgauges        ! locations of measuring gauges
    call setfixedgrids    ! specifies output on arbitrary uniform fixed grids 
    call set_hurricane_params ! Set hurricane parameters
    call set_multilayer_params ! Set multilayer parameters
    
end subroutine setprob
