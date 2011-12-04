c     # common blocks for tide gauge data

      integer OUTGAUGEUNIT
      parameter (maxgauges=1000)
      parameter (OUTGAUGEUNIT=89)
      double precision xgauge(maxgauges), ygauge(maxgauges)
      integer igauge(maxgauges)
      integer mbestsrc(maxgauges), mbestorder(maxgauges)
      double precision t1gauge(maxgauges), t2gauge(maxgauges)
      common /gauges/ xgauge,ygauge,igauge,mbestsrc,mbestorder,mgauges
      common /gauget/ t1gauge,t2gauge

