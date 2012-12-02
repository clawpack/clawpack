subroutine setprob()

    use regions_module, only: set_regions
    use gauges_module, only: set_gauges
    use geoclaw_module, only: set_geo
    use topo_module, only: read_topo_settings
    use dtopo_module, only: read_dtopo_settings
    use qinit_module, only: set_qinit
    use fixedgrids_module, only: set_fixed_grids
    use refinement_module, only: set_refinement
    use storm_module, only: set_storm
    
    implicit none

    call set_geo()                    ! sets basic parameters g and coord system
    call set_refinement()             ! sets refinement control parameters
    call read_topo_settings()         ! specifies topography (bathymetry) files
    call read_dtopo_settings()        ! specifies file with dtopo from earthquake
    call set_qinit()                  ! specifies file with dh if this used instead
    call set_fixed_grids()            ! Fixed grid settings

    call set_storm()                  ! Set storm parameters
    
end subroutine setprob
