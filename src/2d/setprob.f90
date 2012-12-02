subroutine setprob()

    use geoclaw_module
    use topo_module
    use dtopo_module
    use qinit_module
    use fixedgrids_module
    use refinement_module
    
    implicit none

    call set_geo()                    !# sets basic parameters g and coord system
    call read_multilayer_data()       !# Specifies parameters for multiple layers
    call set_refinement()             !# sets refinement control parameters
    call read_topo_settings()         !# specifies topography (bathymetry) files
    call read_dtopo_settings()        !# specifies file with dtopo from earthquake
    call set_qinit()                  !# specifies file with dh if this used instead
    call set_fixed_grids()            !# Fixed grid settings

end subroutine setprob
