c------------------------------------------------------------------------------
        subroutine fgrid_advance(time,dt)
        
c------------------------------------------------------------------------------
        
c:::::::::::::::FGRID_ADVANCE:::::::::::::::::::::::::::::::::::::::::::::::::
c   advance (output) all fgrids at all times that have not yet been output
c   but that have been bracketed by computational times.

      implicit double precision (a-h,o-z)
        
      include  "call.i"
      include  "fixedgrids.i"
      

      tc0=time !# start of computational step
      tcf=tc0+dt !# end of computational step

c     # see if any f-grids should be written out
      do ng=1,mfgrids
        if (tc0.gt.tstartfg(ng).and.ilastoutfg(ng).lt.noutfg(ng)) then
c     # fgrid ng may need to be written out
c     # find the first output number that has not been written out and
c     # find the first output number on a fixed grid that is >= tc0
c     # which will not be written out
           if (dtfg(ng).gt.0.d0) then
             ioutfgend= 1+max(0,nint((tc0-tstartfg(ng))/dtfg(ng)))
           else
             ioutfgend=1
           endif
           ioutfgend=min(ioutfgend,noutfg(ng))
           ioutfgstart=ilastoutfg(ng)+1
c     # write-out fgrid times that are less than tc0, and have not been written yet
c     # these should be the most accurate values at any given point in the fgrid
c     # since tc0> output time
           do ioutfg=ioutfgstart,ioutfgend
             toutfg=tstartfg(ng)+(ioutfg-1)*dtfg(ng)
             if (toutfg.lt.tc0) then
c               # write out the solution for fixed grid ng
                i0=i0fg(ng)
                i02=i0fg2(ng)
c               # test if arrival times should be output
                ioutflag = ioutarrivaltimes(ng)*
     &                         (noutfg(ng)-ilastoutfg(ng))

                call fgridout(fgridearly(i0),fgridlate(i0),
     &              fgridoften(i02),xlowfg(ng),xhifg(ng),ylowfg(ng),
     &              yhifg(ng),mxfg(ng),myfg(ng),
     &              mfgridvars(ng),mfgridvars2(ng),toutfg,
     &              ioutfg,ng,ioutarrivaltimes(ng),ioutflag)

                tlastoutfg(ng)=toutfg
                ilastoutfg(ng)=ilastoutfg(ng)+1
             endif
           enddo

        endif
      enddo
      
      return
      end
