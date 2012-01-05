module gauges_module

      implicit none
      save

!     # common blocks for tide gauge data

      integer, parameter :: OUTGAUGEUNIT = 1000
      integer, parameter :: maxgauges=1000
      integer :: mgauges
      integer :: igauge(maxgauges)
      integer :: mbestsrc(maxgauges), mbestorder(maxgauges)
      
      real(kind=8) :: xgauge(maxgauges), ygauge(maxgauges)
      real(kind=8) :: t1gauge(maxgauges), t2gauge(maxgauges)
      
! -OLD -------------------------     
!     common /gauges/ xgauge,ygauge,igauge,mbestsrc,mbestorder,mgauges
!     common /gauget/ t1gauge,t2gauge

end module
