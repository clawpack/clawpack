!======================================================================
subroutine topo_update(t)
!======================================================================
!
!  called to update topography grids to the current time t
!  topo arrays are modified directly from topo0 and dtopo
!
!  this routine replaces the old method implemented in movetopo
!  where dtopo modified aux arrays directly. Now aux values are
!  always determined from topography grids.
!
!                   David George, dgeorge@usgs.gov, December 20 2013


   use topo_module

   implicit none

   !arguments
   real(kind=8), intent(in) ::t

   !locals
   integer :: i,j,m,mt
   integer :: ij,ij0,irank,idtopo1,idtopo2,jdtopo1,jdtopo2
   integer :: ijll,ijlr,ijul,ijur
   double precision :: x,y,xl,xr,yu,yl,zll,zlr,zul,zur,dz12,dz1,dz2

   if (t<minval(t0dtopo).or.topo_finalized.eqv..true.) then
      return
   endif
   if (minval(topotime)>=maxval(tfdtopo)) then
      deallocate(topo0work)
      topo_finalized = .true.
      return
   endif

   !first find t related values to avoid calculation for every i,j
   do m=1,num_dtopo
      !find t indices
      kdtopo1(m) = int(floor((t-t0dtopo(m))/dtdtopo(m)))+1
      kdtopo2(m) = int(ceiling((t-t0dtopo(m))/dtdtopo(m)))+1
      kdtopo1(m) = min(kdtopo1(m),mtdtopo(m))
      kdtopo2(m) = min(kdtopo2(m),mtdtopo(m))
      kdtopo1(m) = max(kdtopo1(m),1)
      kdtopo2(m) = max(kdtopo2(m),1)
      tdtopo1(m) = t0dtopo(m)+ dtdtopo(m)*real(kdtopo1(m)-1,kind=8) ! tdtopo1<= t
      tdtopo2(m) = t0dtopo(m)+ dtdtopo(m)*real(kdtopo2(m)-1,kind=8) ! tdtopo2>= t
      taudtopo(m) = 1.d0-max(0.d0,((t-tdtopo1(m))/dtdtopo(m)))
      taudtopo(m) = max(taudtopo(m),0.d0)
      index0_dtopowork1(m) = i0dtopo(m) + (kdtopo1(m)-1)*mxdtopo(m)*mydtopo(m)
      index0_dtopowork2(m) = i0dtopo(m) + (kdtopo2(m)-1)*mxdtopo(m)*mydtopo(m)
   enddo


   !first set topofiles aligned exactly with corresponding dtopo
   do i= mtopofiles - num_dtopo + 1, mtopofiles !topofile
      m = i - mtopofiles + num_dtopo !corresponding dtopofile
      !interpolate in time directly for matching nodes
      topowork(i0topo(i):i0topo(i) + mtopo(i)-1) = &
               topo0work(i0topo0(i):i0topo0(i) + mtopo(i)-1) &!initial topo
               + taudtopo(m)*dtopowork(index0_dtopowork1(m):index0_dtopowork1(m) + mtopo(i)-1) &
               + (1.0-taudtopo(m))*dtopowork(index0_dtopowork2(m):index0_dtopowork2(m) + mtopo(i)-1)
      !set time-stamp
      topotime(i) = t
   enddo

   !set non-dtopo associated topofiles
   do mt=1,mtopofiles-num_dtopo
      if (topo0save(mt)<=0) then
         !no intersection or dtopo area already covered by finer topo files
         !shouldn't ever need to alter this topofile
         topotime(mt)=t
         cycle
      endif
      if (topotime(mt)==t) then
         !topofile is already at correct time
         cycle
      endif

      do j=1,mytopo(mt)
         y = yhitopo(mt) - real(j-1,kind=8)*dytopo(mt)
         do i=1,mxtopo(mt)
            ij = i0topo(mt) + (j-1)*mxtopo(mt) + i -1
            ij0 = i0topo0(mt) + (j-1)*mxtopo(mt) + i -1
            x = xlowtopo(mt) +  real(i-1,kind=8)*dxtopo(mt)
            do irank = 1,num_dtopo
               m = mdtopoorder(irank)
               if ( (x>xhidtopo(m)).or.(x<xlowdtopo(m)).or. &
                          (y>yhidtopo(m)).or.(y<ylowdtopo(m))) then
                     !no intersection of point with this dtopo
                     cycle
               endif

               !find indices for bilinear dtopo
               !dtopo arrays are in form of DEM...high y values first
               !note for xy points lying on nodes all indices will be equal
               idtopo1 = int(floor((x-xlowdtopo(m))/dxdtopo(m)))+1
               idtopo2 = int(ceiling((x-xlowdtopo(m))/dxdtopo(m)))+1
               jdtopo1 = int(floor((yhidtopo(m)-y)/dydtopo(m))) + 1
               jdtopo2 = int(ceiling((yhidtopo(m)-y)/dydtopo(m))) + 1
               !indices for dtopo work array
               ijll = index0_dtopowork1(m) + (jdtopo2-1)*mxdtopo(m) + idtopo1 -1
               ijlr = index0_dtopowork1(m) + (jdtopo2-1)*mxdtopo(m) + idtopo2 -1
               ijul = index0_dtopowork1(m) + (jdtopo1-1)*mxdtopo(m) + idtopo1 -1
               ijur = index0_dtopowork1(m) + (jdtopo1-1)*mxdtopo(m) + idtopo2 -1

               !find x,y,z values for bilinear
               !z may be from only 1 or 2 nodes for coincidently aligned grids
               !bilinear should still evaluate correctly
               zll = dtopowork(ijll)
               zlr = dtopowork(ijlr)
               zul = dtopowork(ijul)
               zur = dtopowork(ijur)
               xl = xlowdtopo(m) + real(idtopo1-1,kind=8)*dxdtopo(m)
               xr = xl + dxdtopo(m)
               yu = yhidtopo(m) - real(jdtopo1-1,kind=8)*dydtopo(m)
               yl = yu - dydtopo(m)

               !bilinear value at (x,y) of dtopo cell at t=tdtopo1
               dz1 = zll*(xr-x)*(yu-y) + zlr*(x-xl)*(yu-y) + zul*(xr-x)*(y-yl) + zur*(x-xl)*(y-yl)

               !indices for work array at later time
               ijll = index0_dtopowork2(m) + (jdtopo2-1)*mxdtopo(m) + idtopo1 -1
               ijlr = index0_dtopowork2(m) + (jdtopo2-1)*mxdtopo(m) + idtopo2 -1
               ijul = index0_dtopowork2(m) + (jdtopo1-1)*mxdtopo(m) + idtopo1 -1
               ijur = index0_dtopowork2(m) + (jdtopo1-1)*mxdtopo(m) + idtopo2 -1
               zll = dtopowork(ijll)
               zlr = dtopowork(ijlr)
               zul = dtopowork(ijul)
               zur = dtopowork(ijur)
               !bilinear value of at (x,y) of dtopo cell at t=tdtopo2
               dz2 = zll*(xr-x)*(yu-y) + zlr*(x-xl)*(yu-y) + zul*(xr-x)*(y-yl) + zur*(x-xl)*(y-yl)
               dz12 = taudtopo(m)*dz1 + (1.0-taudtopo(m))*dz2
               dz12 = dz12/(dxdtopo(m)*dydtopo(m))
               !set topo value
               topowork(ij) = topo0work(ij0) + dz12
               !found value from finest dtopo
               exit
            enddo

         enddo
      enddo
      topotime(mt) = t
   enddo


end subroutine topo_update