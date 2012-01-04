
module qinit_module

      implicit none
      save
      
!     # arrays and common blocks for qinit perturbation
!     # specified by setqinit and used in qinit.

      integer, parameter :: maxqinitsize  = 4400000
      integer :: iqinit
      
      real(kind=8) qinitwork(maxqinitsize)
      real(kind=8) xlowqinit
      real(kind=8) ylowqinit
      real(kind=8) tlowqinit
      real(kind=8) xhiqinit
      real(kind=8) yhiqinit
      real(kind=8) thiqinit
      real(kind=8) dxqinit
      real(kind=8) dyqinit
      integer mxqinit
      integer myqinit
      integer minlevelqinit
      integer maxlevelqinit

! -OLD -------------------------
!     common /qinitw/ qinitwork
!     common /qinitparams/
!    &	      xlowqinit,ylowqinit,xhiqinit,yhiqinit,dxqinit,dyqinit,
!    &	      tlowqinit,thiqinit,
!    &	      mxqinit,myqinit,minlevelqinit,
!    &	      maxlevelqinit,mqinitsize,iqinit

end module
