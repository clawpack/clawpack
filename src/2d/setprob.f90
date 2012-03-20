subroutine setprob()

    use geoclaw_module
    use topo_module
    use dtopo_module
    use qinit_module

    implicit none

    call set_geo()          !# sets basic parameters g and coord system
    call set_multilayer()   !# Specifies parameters for multiple layers
    call set_shallow()      !# sets parameters specific to shallow water flows
    call set_topo()         !# specifies topography (bathymetry) files
    call set_dtopo()        !# specifies file with dtopo from earthquake
    call set_qinit()        !# specifies file with dh if this used instead
    call setregions()       !# specifies where refinement is allowed/forced
    call setgauges()        !# locations of measuring gauges
    call setfixedgrids()    !# specifies output on arbitrary uniform fixed grids

end subroutine setprob
