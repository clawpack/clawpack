
      parameter (maxregions=20)

      dimension xlowregion(maxregions)
      dimension ylowregion(maxregions)
      dimension xhiregion(maxregions)
      dimension yhiregion(maxregions)
      dimension minlevelregion(maxregions)
      dimension maxlevelregion(maxregions)
      dimension tlowregion(maxregions)
      dimension thiregion(maxregions)

      common /regionparams/
     &       xlowregion,ylowregion,xhiregion,yhiregion,
     &       tlowregion,thiregion,minlevelregion,maxlevelregion,
     &       mregions

