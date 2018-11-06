!
! :::::::::::::::::::::::::: UPDATE :::::::::::::::::::::::::::::::::
! update - update all grids at level 'level'.
!          this routine assumes cell centered variables.
!          the update is done from 1 level finer meshes under it.
! input parameter:
!    level  - ptr to the only level to be updated. levels coarser than
!             this will be at a diffeent time.
! rewritten in July 2016 from the original update.f
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!

subroutine update (level, nvar, naux)

    use geoclaw_module, only: sea_level, rho
    use amr_module

    use multilayer_module, only: num_layers, dry_tolerance, eta_init

    implicit none
    ! modified for shallow water on topography to use surface level eta
    ! rather than depth h = q(i,j,1) 
    ! eta - q(i,j,1) + aux(i,j,1)

    ! inputs
    integer, intent(in) :: nvar, naux, level

    integer :: ng, levSt, mptr, loc, locaux, nx, ny, mitot, mjtot
    integer :: ilo, jlo, ihi, jhi, mkid, iclo, jclo, ichi, jchi
    integer :: mi, mj, locf, locfaux, iplo, jplo, iphi, jphi
    integer :: iff, jff, nwet(num_layers), ico, jco, i, j, ivar, loccaux, lget
    integer :: listgrids(numgrids(level)), layer, i_layer
    real(kind=8) :: dt, totrat, bc, etasum(num_layers), hsum(num_layers)
    real(kind=8) :: husum(num_layers), hvsum(num_layers)
    real(kind=8) :: hf, bf, huf, hvf, etaf, hav(num_layers), hc(num_layers)
    real(kind=8) :: huc(num_layers), hvc(num_layers), capa, etaav(num_layers)
    real(kind=8) :: capac
    character(len=80) :: String

    lget = level

    String = "(19h    updating level ,i5)"
    
    if (uprint) then
        write(outunit, String) lget
    endif

    dt = possk(lget)

    ! need to set up data structure for parallel distribution of grids
!     call prepgrids(listgrids, numgrids(level), level)


!$OMP PARALLEL DO PRIVATE(ng,mptr,loc,loccaux,nx,ny,mitot,mjtot, &
!$OMP                    ilo,jlo,ihi,jhi,mkid,iclo,ichi, &
!$OMP                    jclo,jchi,mi,mj,locf,locfaux,iplo,iphi, &
!$OMP                    jplo,jphi,iff,jff,totrat,i,j,ivar,capac, &
!$OMP                    capa,bc,etasum,hsum,husum,hvsum, &
!$OMP                    levSt,ico,jco,hf,bf,huf,hvf, &
!$OMP                    etaf,etaav,hav,nwet,hc,huc,hvc, String), &
!$OMP             SHARED(numgrids,listOfGrids,level,intratx,intraty, &
!$OMP                   nghost,uprint,nvar,naux,mcapa,node,listsp, &
!$OMP                   alloc,lstart,dry_tolerance,listStart,lget, &
!$OMP                   num_layers, rho, eta_init), &
!$OMP            DEFAULT(none)


    do ng = 1, numgrids(lget)
        levSt   = listStart(lget)
        mptr    = listOfGrids(levSt+ng-1)
        loc     = node(store1,mptr)
        loccaux = node(storeaux,mptr)
        nx      = node(ndihi,mptr) - node(ndilo,mptr) + 1
        ny      = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
        mitot   = nx + 2*nghost
        mjtot   = ny + 2*nghost
        ilo     = node(ndilo,mptr)
        jlo     = node(ndjlo,mptr)
        ihi     = node(ndihi,mptr)
        jhi     = node(ndjhi,mptr)

        if (node(cfluxptr, mptr) /= 0) then
            call upbnd(alloc(node(cfluxptr,mptr)),alloc(loc),nvar,naux,mitot, &
                    mjtot,listsp(lget),mptr)
        endif

        mkid = lstart(lget+1)
        do while (mkid /= 0)
            iclo   = node(ndilo,mkid)/intratx(lget)
            jclo   = node(ndjlo,mkid)/intraty(lget)
            ichi   = node(ndihi,mkid)/intratx(lget)
            jchi   = node(ndjhi,mkid)/intraty(lget)
            mi      = node(ndihi,mkid)-node(ndilo,mkid) + 1 + 2*nghost
            mj      = node(ndjhi,mkid)-node(ndjlo,mkid) + 1 + 2*nghost
            locf    = node(store1,mkid)
            locfaux = node(storeaux,mkid)

            ! calculate the starting adn ending indices for coarse grid update, if overlap
            iplo = max(ilo,iclo)
            jplo = max(jlo,jclo)
            iphi = min(ihi,ichi)
            jphi = min(jhi,jchi)        

            if (iplo <= iphi .and. jplo <= jphi) then
                !calculate the starting index for the fine grid source pts.
                iff    = iplo*intratx(lget) - node(ndilo,mkid) + nghost + 1
                jff    = jplo*intraty(lget) - node(ndjlo,mkid) + nghost + 1
                totrat = intratx(lget) * intraty(lget)

                do i = iplo-ilo+nghost+1, iphi-ilo+nghost+1
                    do j = jplo-jlo+nghost+1, jphi-jlo+nghost+1
                        if (uprint) then
                            String = "(' updating pt. ',2i4,' of grid ',i3,' using ',2i4,' of grid ',i4)"
                            write(outunit,String) i,j,mptr,iff,jff,mkid                
                            
                            String = "(' old vals: ',4e25.15)"
                            write(outunit,String)(alloc(iadd(ivar,i,j,loc,mitot)),ivar=1,nvar)
                        endif
                        
                        if (mcapa == 0) then
                            capac = 1.0d0
                        else
                            capac = alloc(iaddcaux(i,j,loccaux,mcapa,mitot))
                        endif

                        bc = alloc(iaddctopo(i,j,loccaux,mitot))

                        etasum = 0.d0
                        hsum = 0.d0
                        husum = 0.d0
                        hvsum = 0.d0

                        nwet = 0

                        do jco = 1, intraty(lget)
                            do ico = 1, intratx(lget)
                                    
                                if (mcapa == 0) then
                                    capa = 1.0d0
                                else
                                    capa = alloc(iaddfaux(iff+ico-1,jff+jco-1,locfaux,mcapa,mi))
                                endif
                                
                                bf = alloc(iaddftopo(iff+ico-1,jff+jco-1,locfaux,mi))*capa
                                
                                do layer = num_layers, 1, -1

                                    hf = alloc(iaddf(3*layer-2,iff+ico-1,jff+jco-1,locf,mi))*capa / rho(layer)
                                    huf= alloc(iaddf(3*layer-1,iff+ico-1,jff+jco-1,locf,mi))*capa 
                                    hvf= alloc(iaddf(3*layer,iff+ico-1,jff+jco-1,locf,mi))*capa 

                                    if ( alloc(iaddf(3*layer-2,iff+ico-1,jff+jco-1,locf,mi))/rho(layer) &
                                        > dry_tolerance(layer)) then
                                        etaf = hf + bf
                                        nwet(layer) = nwet(layer) + 1
                                    else
                                        ! Set to zero for correct averages below.
                                        etaf = 0.d0
                                        huf=0.d0
                                        hvf=0.d0
                                    endif

                                    hsum(layer)   = hsum(layer) + hf
                                    husum(layer)  = husum(layer) + huf
                                    hvsum(layer)  = hvsum(layer) + hvf
                                    etasum(layer) = etasum(layer) + etaf  

                                    bf = bf + hf   
                                enddo
                            enddo
                        enddo

                        do layer = num_layers, 1, -1
                            if (nwet(layer) > 0) then
                                etaav(layer) = etasum(layer)/dble(nwet(layer))
                                hav(layer) = hsum(layer)/dble(nwet(layer))
                                hc(layer) = min(hav(layer), (max(etaav(layer)-bc*capac, 0.0d0)))
                                huc(layer) = (min(hav(layer), hc(layer)) / hsum(layer) * husum(layer))
                                hvc(layer) = (min(hav(layer), hc(layer)) / hsum(layer) * hvsum(layer))
                            else
                                hc(layer) = 0.0d0
                                huc(layer) = 0.0d0
                                hvc(layer) = 0.0d0
                            endif

                            alloc(iadd(3*layer-2,i,j,loc,mitot)) = hc(layer) / capac * rho(layer)
                            alloc(iadd(3*layer-1,i,j,loc,mitot)) = huc(layer) / capac 
                            alloc(iadd(3*layer,i,j,loc,mitot)) = hvc(layer) / capac 

                            bc = bc + hc(layer)
                        enddo

                        if (uprint) then
                            String = "(' new vals: ',4e25.15)"
                            write(outunit, String)(alloc(iadd(ivar, i, j,loc,mitot)), ivar=1, nvar)
                        endif

                        jff = jff + intraty(lget)
                    enddo

                    iff = iff + intratx(lget)
                    jff = jplo * intraty(lget) - node(ndjlo, mkid) + nghost + 1
                enddo

                
            endif
            mkid = node(levelptr, mkid)

        enddo
    enddo
!$OMP END PARALLEL DO
    return


contains

    integer pure function iadd(ivar,i,j,loc,mitot)
        integer, intent(in) :: i, j, ivar, loc, mitot
        iadd = loc + ivar-1 + nvar*((j-1)*mitot+i-1)
    end function iadd

    integer pure function iaddf(ivar,i,j,locf,mi)
        integer, intent(in) :: i, j, ivar, locf, mi
        iaddf = locf   + ivar-1 + nvar*((j-1)*mi+i-1)
    end function iaddf

    integer pure function iaddfaux(i,j,locfaux,mcapa,mi)
        integer, intent(in) :: i, j, locfaux, mcapa, mi
        iaddfaux = locfaux + mcapa-1 + naux*((j-1)*mi + (i-1))
    end function iaddfaux

    integer pure function iaddcaux(i,j,loccaux,mcapa,mitot)
        integer, intent(in) :: i, j, loccaux, mcapa, mitot
        iaddcaux = loccaux + mcapa-1 + naux*((j-1)*mitot+(i-1))
    end function iaddcaux

    integer pure function iaddftopo(i,j,locfaux,mi)
        integer, intent(in) :: i, j, locfaux, mi
        iaddftopo = locfaux + naux*((j-1)*mi + (i-1))
    end function iaddftopo

    integer pure function iaddctopo(i,j,loccaux, mitot)
        integer, intent(in) :: i, j, loccaux, mitot
        iaddctopo = loccaux + naux*((j-1)*mitot+(i-1))
    end function iaddctopo

end subroutine
