
c
c -----------------------------------------------------------------------------
c
      double precision function get_max_speed(val,mitot,mjtot,nvar,
     &                                        aux,naux,nghost,hx,hy)

      use geoclaw_module
      implicit double precision (a-h,o-z)
      dimension   val(nvar,mitot,mjtot), aux(naux,mitot,mjtot)


      sp_over_h = 0.d0   ! compute max speed over h, since dx may not equal dy
      if (icoordsys .eq. 2) then
        do j = nghost+1, mjtot-nghost
         ymetric = Rearth*pi/180.d0
         hyphys = ymetric*hy

         do i = nghost+1, mitot-nghost
            xmetric = cos(aux(3,i,j))*Rearth*pi/180.d0
            hxphys = xmetric * hx
            h  = val(1,i,j)
            if (h .gt. drytolerance) then
              u  = val(2,i,j)/h
              v  = val(3,i,j)/h
            else
              u = 0.d0
              v = 0.d0
            endif
            sig = sqrt(grav*h)
            sp_over_h = max((abs(u)+sig)/hxphys,(abs(v)+sig)/hyphys,
     &                      sp_over_h)
         end do
         end do
      else  ! speeds in cartesian coords, no metrics needed
         do j = nghost+1, mjtot-nghost
         do i = nghost+1, mitot-nghost
            h  = val(1,i,j)
            if (h .gt. drytolerance) then
              u  = val(2,i,j)/h
              v  = val(3,i,j)/h
            else
              u = 0.d0
              v = 0.d0
            endif
            sig = sqrt(grav*h)
            sp_over_h = max((abs(u)+sig)/hx,(abs(v)+sig)/hy,sp_over_h)
         end do
         end do
      endif

       get_max_speed = sp_over_h

       return
       end
