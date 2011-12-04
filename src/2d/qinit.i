
c     # arrays and common blocks for qinit perturbation
c     # specified by setqinit and used in qinit.

      parameter (maxqinitsize=4400000)

      dimension qinitwork(maxqinitsize)

      integer iqinit

      real*8 xlowqinit
      real*8 ylowqinit
      real*8 tlowqinit
      real*8 xhiqinit
      real*8 yhiqinit
      real*8 thiqinit
      real*8 dxqinit
      real*8 dyqinit

      real*8 mxqinit
      real*8 myqinit

      real*8 minlevelqinit
      real*8 maxlevelqinit


      common /qinitw/ qinitwork
      common /qinitparams/
     &	      xlowqinit,ylowqinit,xhiqinit,yhiqinit,dxqinit,dyqinit,
     &	      tlowqinit,thiqinit,
     &	      mxqinit,myqinit,minlevelqinit,
     &	      maxlevelqinit,mqinitsize,iqinit

