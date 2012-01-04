
module fixedgrids_module

!     # arrays and common blocks for fixed output grids

      implicit none
      save

      integer, parameter :: maxfgrids=1
      integer, parameter :: maxfgridsize=6
      integer :: mfgrids
      integer :: ioutarrivaltimes(maxfgrids)
      integer :: ioutsurfacemax(maxfgrids)
      integer :: mxfg(maxfgrids)
      integer :: myfg(maxfgrids)
      integer :: i0fg(maxfgrids)
      integer :: i0fg2(maxfgrids)
      integer :: noutfg(maxfgrids)
      integer :: ilastoutfg(maxfgrids)
      integer :: mfgridvars(maxfgrids)
      integer :: mfgridvars2(maxfgrids)
      
      real(kind=8) fgridearly(maxfgridsize)
      real(kind=8) fgridlate(maxfgridsize)
      real(kind=8) fgridoften(maxfgridsize)
     
      real(kind=8) xlowfg(maxfgrids)
      real(kind=8) xhifg(maxfgrids)
      real(kind=8) ylowfg(maxfgrids)
      real(kind=8) yhifg(maxfgrids)

      
      real(kind=8) tlastoutfg(maxfgrids)
      real(kind=8) tstartfg(maxfgrids)
      real(kind=8) tendfg(maxfgrids)
    
      real(kind=8) dtfg(maxfgrids)
      real(kind=8) dxfg(maxfgrids)
      real(kind=8) dyfg(maxfgrids)
      
      real(kind=8) tcfmax

! OLD ==================================================================
!
!     common /fgrid1/ fgridearly
!     common /fgrid2/ fgridlate
!     common /fgrid3/ fgridoften
!     common /fgridparams/
!    &       ylowfg, xhifg, yhifg, tstartfg, tendfg, dtfg, dxfg, dyfg,
!    &       tlastoutfg, tcfmax, xlowfg, i0fg, mxfg, myfg, noutfg, 
!    &       ilastoutfg, ioutarrivaltimes, ioutsurfacemax, i0fg2, 
!    &       mfgridvars, mfgridvars2, mfgrids

end module
