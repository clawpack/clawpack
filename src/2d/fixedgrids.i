
c     # arrays and common blocks for fixed output grids

      parameter (maxfgrids=3)
      parameter (maxfgridsize=404012)

      dimension fgridearly(maxfgridsize)
      dimension fgridlate(maxfgridsize)
      dimension fgridoften(maxfgridsize)

      dimension xlowfg(maxfgrids)
      dimension xhifg(maxfgrids)
      dimension ylowfg(maxfgrids)
      dimension yhifg(maxfgrids)
      dimension mxfg(maxfgrids)
      dimension myfg(maxfgrids)

      dimension i0fg(maxfgrids)
      dimension i0fg2(maxfgrids)
      dimension noutfg(maxfgrids)
      dimension ilastoutfg(maxfgrids)
      dimension mfgridvars(maxfgrids)
      dimension mfgridvars2(maxfgrids)
      dimension ioutarrivaltimes(maxfgrids)
      dimension ioutsurfacemax(maxfgrids)

      dimension tlastoutfg(maxfgrids)
      dimension tstartfg(maxfgrids)
      dimension tendfg(maxfgrids)

      dimension dtfg(maxfgrids)
      dimension dxfg(maxfgrids)
      dimension dyfg(maxfgrids)

c===================================================================

      common /fgrid1/ fgridearly
      common /fgrid2/ fgridlate
      common /fgrid3/ fgridoften
      common /fgridparams/
     &       ylowfg, xhifg, yhifg, tstartfg, tendfg, dtfg, dxfg, dyfg,
     &       tlastoutfg, tcfmax, xlowfg, i0fg, mxfg, myfg, noutfg, 
     &       ilastoutfg, ioutarrivaltimes, ioutsurfacemax, i0fg2, 
     &       mfgridvars, mfgridvars2, mfgrids
