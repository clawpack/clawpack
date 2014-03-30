subroutine setprob()

    use geoclaw_module
    use topo_module
    use qinit_module
    use fixedgrids_module
    use refinement_module

    use friction_module
    use multilayer_module
    
    implicit none

    integer :: total_aux

    call set_geo()                    ! sets basic parameters g and coord system
    call set_refinement()             ! sets refinement control parameters
    call read_dtopo_settings()        ! specifies file with dtopo from earthquake
    call read_topo_settings()         ! specifies topography (bathymetry) files
    call set_qinit()                  ! specifies file with dh if this used instead
    call set_fixed_grids()            ! Fixed grid settings

    call setup_variable_friction()    ! Variable friction parameter
    call set_multilayer()             ! Set multilayer parameters

end subroutine setprob
