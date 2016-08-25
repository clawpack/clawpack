c
c ----------------------------------------------------------------
c Modified from auxcoarse.f in amrclaw
c Only coarsens topography
c ----------------------------------------------------------------
c
       subroutine auxcoarsen(auxdub,midub,mjdub,auxbgc,
     1                       mi2tot,mj2tot,naux,auxtype)

       use amr_module, only: mcapa
       implicit double precision (a-h, o-z)

       dimension     auxdub(naux,midub, mjdub)
       dimension     auxbgc(naux,mi2tot,mj2tot)
       character*10  auxtype(naux)

c :::::::::::::::::::::::: COARSEN ::::::::::::::::::::::::::::::::
c coarsen = coarsen the fine grid auxiliary data (with double the usual
c           number of ghost cells to prepare coarsened data
c           for error estimation.
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

       do 20 j = 1, mj2tot
         jfine = 2*(j-1) + 1
         do 20 i = 1, mi2tot
           ifine = 2*(i-1) + 1

           if (mcapa .eq. 0) then
             capac=1.0d0
           else
             capac=auxbgc(2,i,j)
           endif

           toposum = 0.d0
           do jco = 1,2
             do ico = 1,2
               if (mcapa .eq. 0) then
                 capa=1.0d0
               else
                 capa=auxdub(2,ifine+ico-1,jfine+jco-1)
               endif

               topof=auxdub(1,ifine+ico-1,jfine+jco-1)*capa

               toposum = toposum + topof
             enddo
           enddo


           auxbgc(1,i,j) = toposum/(4.d0 * capac)
20       continue


       return
       end
