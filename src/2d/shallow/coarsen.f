c
c ---------------------------------------------------
c Modified from coarsen.f in amrclaw
c Uses code taken from update.f to take into account
c wet and dry cells.
c ---------------------------------------------------
c
       subroutine coarsen(valdub,midub,mjdub,auxdub,
     &                    valbgc,mi2tot,mj2tot,auxbgc,nvar,naux)


       use geoclaw_module, only: dry_tolerance
c
       use amr_module, only: mcapa,nghost
       implicit none

       integer, intent(in) :: midub,mjdub,mi2tot,mj2tot,nvar,naux
       real(kind=8), intent(in) :: valdub(nvar,midub, mjdub)
       real(kind=8), intent(in) :: auxdub(naux,midub, mjdub)
       real(kind=8), intent(inout) :: valbgc(nvar,mi2tot,mj2tot)
       real(kind=8), intent(inout) :: auxbgc(naux,mi2tot,mj2tot)

c      # Local variables
       integer :: j,jfine,i,ifine,jco,ico,nwet
       real(kind=8) :: capac,bc,etasum,hsum,husum,hvsum,capa
       real(kind=8) :: hf,huf,hvf,bf,etaf,etaav,hav,hc,huc,hvc

c :::::::::::::::::::::::: COARSEN ::::::::::::::::::::::::::::::::
c coarsen = coarsen the fine grid data (with double the usual
c           number of ghost cells to prepare coarsened
c           grid for error estimation.
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
       do j = 1, mj2tot

         jfine = 2*(j-1) + 1

         do i = 1, mi2tot
           ifine = 2*(i-1) + 1

c          # Getting coarse capacity function and bathymetry
           if (mcapa .eq. 0) then
               capac=1.0d0
           else
               capac=auxbgc(2,i,j)
           endif

           bc = auxbgc(1,i,j)

c          # Calculating which values to use for coarsen
           etasum = 0.d0
           hsum = 0.d0
           husum = 0.d0
           hvsum = 0.d0

           nwet=0

           do jco = 1, 2
            do ico = 1, 2
               if (mcapa .eq. 0) then
                   capa=1.0d0
               else
                   capa=auxdub(2,ifine+ico-1,jfine+jco-1)
               endif

               hf =valdub(1,ifine+ico-1,jfine+jco-1)*capa
               huf=valdub(2,ifine+ico-1,jfine+jco-1)*capa
               hvf=valdub(3,ifine+ico-1,jfine+jco-1)*capa

               bf =auxdub(1,ifine+ico-1,jfine+jco-1)*capa

               if (hf > dry_tolerance) then
                   etaf = hf+bf
                   nwet=nwet+1
               else
                   etaf = 0.d0
                   huf=0.d0
                   hvf=0.d0
                   bf = 0.d0
               endif

               hsum   = hsum + hf
               husum  = husum + huf
               hvsum  = hvsum + hvf
               etasum = etasum + etaf
             enddo
           enddo

c          # Getting coarse averages
           if (nwet.gt.0) then
               etaav=etasum/dble(nwet)
               hav= hsum/dble(nwet)

c              hc=max(etaav-bc*capac,0.d0) !tsunamiclaw method
               hc=min(hav,(max(etaav-bc*capac,0.d0)))
               if (hsum.ne.0) then
                   huc=(min(hav,hc)/hsum)*husum
                   hvc=(min(hav,hc)/hsum)*hvsum
               else
                   huc= 0.d0
                   hvc= 0.d0
               endif
           else
               hc=0.d0
               huc=0.d0
               hvc=0.d0
           endif

c          # set h on coarse grid based on surface, not conservative near shoreline
           valbgc(1,i,j) = hc / capac
           valbgc(2,i,j) = huc / capac
           valbgc(3,i,j) = hvc / capac

         enddo
       enddo

       return
       end


