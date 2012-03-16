subroutine qinit(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
    
    use qinit_module
    use geoclaw_module
    
    implicit none
    
    ! Subroutine arguments
    integer, intent(in) :: maxmx,maxmy,meqn,mbc,mx,my,maux
    double precision(in) :: xlower,ylower,dx,dy
    double precision(inout) :: q(meqn,1-mbc:maxmx+mbc,1-mbc:maxmy+mbc)
    double precision(inout) :: aux(maux,1-mbc:maxmx+mbc,1-mbc:maxmy+mbc)
    
    ! Locals
    integer :: i,j
    
    ! Set flat state based on eta_init
    do i=1,mx
        do j=1,my
            ! This loop only executes if layers > 1, multiply by rho then
            do m=1,layers-1
                layer_index = 3*(m-1) + 1
                q(layer_index,i,j) = max(0.d0,eta_init(m)-eta_init(m+1)) * rho(m)
                q(layer_index+1,i,j) = 0.d0
                q(layer_index+2,i,j) = 0.d0
            enddo
            ! Take care of last layer
            bottom_layer = 3*layers-2
            q(bottom_layer,i,j) = max(0.d0,eta_init(m) - aux(1,i,j))
            q(bottom_layer+1,i,j) = 0.d0
            q(bottom_layer+2,i,j) = 0.d0
            
            ! Multiply by rho if this is a multi-layer run
            if (layers > 1) then
                q(bottom_layer,i,j) = q(bottom_layer,i,j) * rho(layers)
            endif
        enddo
    enddo
    
    ! Add perturbation to surfaces (only handles layers == 1 case)
    if (iqinit > 0) then
        do i=1-mbc,mx+mbc
            x = xlower + (i-0.5d0)*dx
            xim = x - 0.5d0*dx
            xip = x + 0.5d0*dx
            do j=1-mbc,my+mbc
                y = ylower + (j-0.5d0)*dy
                yjm = y - 0.5d0*dy
                yjp = y + 0.5d0*dy


                if ((xip > xlowqinit).and.(xim < xhiqinit).and.  &
                    (yjp > ylowqinit).and.(yjm < yhiqinit)) then

                    xipc=min(xip,xhiqinit)
                    ximc=max(xim,xlowqinit)
                    xc=0.5d0*(xipc+ximc)

                    yjpc=min(yjp,yhiqinit)
                    yjmc=max(yjm,ylowqinit)
                    yc=0.5d0*(yjmc+yjpc)

                    dq = topointegral(ximc,xc,xipc,yjmc,yc,yjpc,xlowqinit,ylowqinit,dxqinit,dyqinit,mxqinit,myqinit,qinitwork,1)
                    dq = dq / ((xipc-ximc)*(yjpc-yjmc)*aux(2,i,j))

                    if (iqinit < 4) then 
                        if (aux(1,i,j) <= 0.d0) then
                            q(iqinit,i,j) = q(iqinit,i,j) + dq
                        endif
                    else if (iqinit == 4) then
                      q(1,i,j) = max(dq-aux(1,i,j),0.d0)
                    endif
                endif
            enddo
        enddo
    endif
    
end subroutine qinit