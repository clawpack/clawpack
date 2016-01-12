#get the directory of this makefile
GEOLIB:=$(dir $(lastword $(MAKEFILE_LIST)))

#import AMR library dependencies
AMRLIB = $(CLAW)/amrclaw/src/2d
include $(AMRLIB)/Makefile.amr_2d_geoclaw

#list of common sources for amr 2d codes
COMMON_MODULES += \
 $(GEOLIB)/utility_module.f90 \
 $(GEOLIB)/geoclaw_module.f90 \
 $(GEOLIB)/gauges_module.f90 \
 $(GEOLIB)/topo_module.f90 \
 $(GEOLIB)/qinit_module.f90 \
 $(GEOLIB)/refinement_module.f90 \
 $(GEOLIB)/fixedgrids_module.f90 \
 $(GEOLIB)/fgmax_module.f90 \
 $(GEOLIB)/surge/holland_storm_module.f90 \
 $(GEOLIB)/surge/stommel_storm_module.f90 \
 $(GEOLIB)/surge/constant_storm_module.f90 \
 $(GEOLIB)/surge/storm_module.f90 \
 $(GEOLIB)/friction_module.f90

COMMON_SOURCES += \
  $(GEOLIB)/setprob.f90 \
  $(GEOLIB)/qinit.f90 \
  $(GEOLIB)/topo_update.f90 \
  $(GEOLIB)/cellgridintegrate2.f \
  $(GEOLIB)/topointegral.f \
  $(GEOLIB)/bilinearintegral.f \
  $(GEOLIB)/stepgrid.f \
  $(GEOLIB)/src2.f90 \
  $(GEOLIB)/src1d.f90 \
  $(GEOLIB)/step2.f90 \
  $(GEOLIB)/flux2fw.f \
  $(GEOLIB)/qad.f \
  $(GEOLIB)/valout.f \
  $(GEOLIB)/filval.f90 \
  $(GEOLIB)/filpatch.f90 \
  $(GEOLIB)/bc2amr.f \
  $(GEOLIB)/update.f \
  $(GEOLIB)/setaux.f90 \
  $(GEOLIB)/flag2refine2.f90  \
  $(GEOLIB)/allowflag.f90  \
  $(GEOLIB)/b4step2.f90 \
  $(GEOLIB)/upbnd.f  \
  $(GEOLIB)/tick.f \
  $(GEOLIB)/setgrd.f \
  $(GEOLIB)/gfixup.f \
  $(GEOLIB)/ginit.f \
  $(GEOLIB)/getmaxspeed.f90 \
  $(GEOLIB)/advanc.f \
  $(GEOLIB)/amr2.f90 \
  $(GEOLIB)/fgmax_read.f90 \
  $(GEOLIB)/fgmax_frompatch.f90 \
  $(GEOLIB)/fgmax_interpolate.f90 \
  $(GEOLIB)/fgmax_values.f90 \
  $(GEOLIB)/fgmax_finalize.f90 \
