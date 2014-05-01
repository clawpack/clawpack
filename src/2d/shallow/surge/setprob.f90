subroutine setprob()

    use regions_module, only: set_regions
    use gauges_module, only: set_gauges
    use geoclaw_module, only: set_geo
    use topo_module, only: read_topo_settings, read_dtopo_settings
    use qinit_module, only: set_qinit
    use fixedgrids_module, only: set_fixed_grids
    use refinement_module, only: set_refinement
    use storm_module, only: set_storm
    use friction_module, only: setup_variable_friction
    
    implicit none

    call set_geo()                    ! sets basic parameters g and coord system
    call set_refinement()             ! sets refinement control parameters
    call read_topo_settings()         ! specifies topography (bathymetry) files
    call read_dtopo_settings()        ! specifies file with dtopo from earthquake
    call set_qinit()                  ! specifies file with dh if this used instead
    call set_fixed_grids()            ! Fixed grid settings

    call set_storm()                  ! Set storm parameters
    call setup_variable_friction()    ! Set variable friction parameters
    
end subroutine setprob
