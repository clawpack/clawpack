module regions_module

      implicit none
      save
      
      integer, parameter :: maxregions = 20
      integer mregions
      integer :: minlevelregion(maxregions)
      integer :: maxlevelregion(maxregions)
      
      
      real(kind=8) :: xlowregion(maxregions)
      real(kind=8) :: ylowregion(maxregions)
      real(kind=8) :: xhiregion(maxregions)
      real(kind=8) :: yhiregion(maxregions)
      real(kind=8) :: tlowregion(maxregions)
      real(kind=8) :: thiregion(maxregions)
      

!--OLD---------------------------------------------
!      common /regionparams/
!     &       xlowregion,ylowregion,xhiregion,yhiregion,
!     &       tlowregion,thiregion,minlevelregion,maxlevelregion,
!     &       mregions

end module
