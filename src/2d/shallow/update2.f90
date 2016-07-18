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

    use geoclaw_module, only: dry_tolerance
    use amr_module
    implicit none
    ! modified for shallow water on topography to use surface level eta
    ! rather than depth h = q(i,j,1) 
    ! eta - q(i,j,1) + aux(i,j,1)

    ! inputs
    integer intent(in) :: nvar, naux
    real(kind=8) intent(in) :: level

    integer :: ng, levSt, mptr, loc, locaux, nx, ny, mitot, mjtot
    integer :: ilo, jlo, ihi, jhi, mkid
    integer :: listgrids(numgrids(level))
    real(kind=8) :: lget, dt

    lget = level

    character(len=80) :: String
    String = "(19h    updating level ,i5)"
    
    if (uprint)
        write(outunit, String) leget
    endif

    dt = possk(leget)

!$OMP PARALLEL DO PRIVATE(ng,mptr,loc,loccaux,nx,ny,mitot,mjtot,
!$OMP&                    ilo,jlo,ihi,jhi,locuse,mkid,iclo,ichi,
!$OMP&                    jclo,jchi,mi,mj,locf,locfaux,iplo,iphi,
!$OMP&                    jplo,jphi,iff,jff,totrat,i,j,ivar,capac,
!$OMP&                    capa,bc,etasum,hsum,husum,hvsum,drytol,
!$OMP&                    newt,levSt,ico,jco,hf,bf,huf,hvf,
!$OMP&                    etaf,etaav,hav,nwet,hc,huc,hvc),
!$OMP&             SHARED(numgrids,listOfGrids,level,intratx,intraty,
!$OMP&                   nghost,uprint,nvar,naux,mcapa,node,listsp,
!$OMP&                   alloc,lstart,dry_tolerance,listStart,lget),
!$OMP&            DEFAULT(none)

    ! need to set up data structure for parallel distribution of grids
    ! call prepgrids(listgrids, numgrids(level), level)

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

        if (node(cfluxptr, mptr) == 0) 
            mkid = lstart(lget+1)
        endif

        call upbnd(alloc(nod(cfluxptr,mprt),alloc(loc),nvar,naux,mitot, &
                    mjtot, listsp(leget),mptr), mitot, mjtot, listsp(leget), &
                    alloc(locuse),mptr)

        































    integer pure function iadd(ivar,i,j)
        integer, intent(in) :: i,j
        iadd = loc + ivar-1 + nvar*((j-1)*mitot+i-1)
    end function iadd

    integer pure function iaddf(ivar,i,j) 
        integer intent(in) :: i,j
        iaddf = locf   + ivar-1 + nvar*((j-1)*mi+i-1)
    end function iaddf

    integer pure function iaddfaux(i,j)
        integer intent(in) :: i,j
        iaddfaux = locfaux + mcapa-1 + naux*((j-1)*mi + (i-1))
    end function iaddfaux

    integer pure function iaddcaux(i,j)
        integer intent(in) :: i,j
        iaddcaux = loccaux + mcapa-1 + naux*((j-1)*mitot+(i-1))
    end function iaddcaux

    integer pure funtion iaddftopo(i,j)
        integer intent(in) :: i,j
        iaddftopo = locfaux +  naux*((j-1)*mi + (i-1))
    end function iaddftopo

    integer pure funtion iaddctopo(i,j)
        integer intent(in) :: i,j
        iaddctopo = loccaux +  naux*((j-1)*mitot+(i-1))
    end function iaddctopo
