
# Makefile for Clawpack code in this directory.
# This version only sets the local files and frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
CLAWMAKE = $(CLAWUTIL)/src/Makefile.common

# See the above file for details and a list of make options, or type
#   make .help
# at the unix prompt.


# Adjust these variables if desired:
# ----------------------------------

CLAW_PKG = geoclaw                  # Clawpack package to use
EXE = xgeoclaw                 # Executable to create
SETRUN_FILE = setrun.py        # File containing function to make data
OUTDIR = _output               # Directory for output
SETPLOT_FILE = setplot.py      # File containing function to set plots
PLOTDIR = _plots               # Directory for plots

# Environment variable FC should be set to fortran compiler, e.g. gfortran
FC ?= gfortran
FFLAGS ?=
MAKELEVEL ?= 0
ifeq ($(MAKELEVEL),0)
	# If you want to include NetCDF output add these flags
	# Note that the NETCDF_INC should point to the location of the hearder files
	# for NetCDF
	# NETCDF_INC ?= /usr/include
	# FFLAGS += -I$(NETCDF_INC) -lnetcdf -lnetcdff

	# Threading flags
	# FFLAGS += $(OMP_FLAG)
endif

# ---------------------------------
# List of sources for this program:
# ---------------------------------

AMRLIB = $(AMRCLAW)/src/2d
GEOLIB = $(GEOCLAW)/src/2d

MODULES = \
 $(GEOLIB)/utility_module.f90 \
 $(AMRLIB)/amr_module.f90 \
 $(AMRLIB)/regions_module.f90 \
 $(AMRLIB)/gauges_module.f90 \
 $(GEOLIB)/geoclaw_module.f90 \
 $(GEOLIB)/topo_module.f90 \
 $(GEOLIB)/dtopo_module.f90 \
 $(GEOLIB)/qinit_module.f90 \
 $(GEOLIB)/refinement_module.f90 \
 $(GEOLIB)/fixedgrids_module.f90 \
 $(GEOLIB)/storm_surge/holland_storm_module.f90 \
 $(GEOLIB)/storm_surge/constant_storm_module.f90 \
 $(GEOLIB)/storm_surge/stommel_storm_module.f90 \
 $(GEOLIB)/storm_surge/storm_module.f90

# Switch source file grdfit2.f and grdfit.f to change flagging strategy
SOURCES = \
  $(GEOLIB)/storm_surge/src2.f90 \
  $(GEOLIB)/storm_surge/src1d.f90 \
  $(GEOLIB)/storm_surge/setaux.f90 \
  $(GEOLIB)/storm_surge/setprob.f90 \
  $(GEOLIB)/storm_surge/qinit.f90 \
  $(GEOLIB)/storm_surge/flag2refine.f90 \
  $(GEOLIB)/storm_surge/flag2refine2.f90 \
  $(GEOLIB)/storm_surge/b4step2.f90 \
  $(GEOLIB)/storm_surge/tick.f \
  $(GEOLIB)/storm_surge/valout.f \
  $(GEOLIB)/bc2amr.f \
  $(GEOLIB)/movetopo.f \
  $(GEOLIB)/allowflag.f90 \
  $(GEOLIB)/cellgridintegrate.f \
  $(GEOLIB)/topointegral.f \
  $(GEOLIB)/bilinearintegral.f \
  $(GEOLIB)/stepgrid.f \
  $(GEOLIB)/rpn2.f \
  $(GEOLIB)/riemannsolvers.f \
  $(GEOLIB)/rpt2.f \
  $(GEOLIB)/step2.f \
  $(GEOLIB)/flux2fw.f \
  $(GEOLIB)/qad.f \
  $(GEOLIB)/filval_new.f90 \
  $(GEOLIB)/filpatch_new.f90 \
  $(GEOLIB)/update.f \
  $(GEOLIB)/dumpgauge.f \
  $(GEOLIB)/upbnd.f  \
  $(GEOLIB)/setgrd.f \
  $(GEOLIB)/gfixup.f \
  $(GEOLIB)/ginit.f \
  $(GEOLIB)/getmaxspeed.f90 \
  $(GEOLIB)/advanc.f \
  $(AMRLIB)/amr2.f90 \
  $(AMRLIB)/prefilp_new.f90 \
  $(AMRLIB)/quick_sort1.f \
  $(AMRLIB)/estdt.f \
  $(AMRLIB)/bound_new.f90 \
  $(AMRLIB)/check4nans.f \
  $(AMRLIB)/spest.f  \
  $(AMRLIB)/spest2.f  \
  $(AMRLIB)/init_iflags.f \
  $(AMRLIB)/igetsp.f \
  $(AMRLIB)/reclam.f \
  $(AMRLIB)/birect.f \
  $(AMRLIB)/cleanup.f \
  $(AMRLIB)/colate.f \
  $(AMRLIB)/colate2.f \
  $(AMRLIB)/errest.f \
  $(AMRLIB)/bufnst.f \
  $(AMRLIB)/bufnst2.f \
  $(AMRLIB)/errf1.f \
  $(AMRLIB)/fixcapaq.f \
  $(AMRLIB)/flglvl.f \
  $(AMRLIB)/flglvl2.f \
  $(AMRLIB)/fluxad.f \
  $(AMRLIB)/fluxsv.f \
  $(AMRLIB)/grdfit2.f \
  $(AMRLIB)/intfil_new.f90 \
  $(AMRLIB)/moment.f \
  $(AMRLIB)/nestck.f \
  $(AMRLIB)/nestck2.f \
  $(AMRLIB)/prepc.f \
  $(AMRLIB)/prepf.f \
  $(AMRLIB)/projec.f \
  $(AMRLIB)/projec2.f \
  $(AMRLIB)/signs.f \
  $(AMRLIB)/findcut.f \
  $(AMRLIB)/smartbis.f \
  $(AMRLIB)/putnod.f \
  $(AMRLIB)/putsp.f \
  $(AMRLIB)/regrid.f \
  $(AMRLIB)/setuse.f \
  $(AMRLIB)/stst1.f \
  $(AMRLIB)/trimbd_new.f90 \
  $(AMRLIB)/nodget.f \
  $(AMRLIB)/basic.f \
  $(AMRLIB)/outval.f \
  $(AMRLIB)/copysol.f \
  $(AMRLIB)/outvar.f \
  $(AMRLIB)/outmsh.f \
  $(AMRLIB)/outtre.f \
  $(AMRLIB)/domain.f  \
  $(AMRLIB)/cellave.f \
  $(AMRLIB)/fdisc.f \
  $(AMRLIB)/fss.f \
  $(AMRLIB)/zeroin.f \
  $(AMRLIB)/setflags.f \
  $(AMRLIB)/shiftset.f \
  $(AMRLIB)/shiftset2.f \
  $(AMRLIB)/conck.f \
  $(AMRLIB)/domshrink.f \
  $(AMRLIB)/domprep.f \
  $(AMRLIB)/domup.f \
  $(AMRLIB)/domcopy.f \
  $(AMRLIB)/coarsen.f \
  $(AMRLIB)/setdomflags.f \
  $(AMRLIB)/setIndices.f \
  $(AMRLIB)/coarseGridFlagSet.f \
  $(AMRLIB)/setPhysBndryFlags.f \
  $(AMRLIB)/griddomup.f \
  $(AMRLIB)/griddomcopy.f \
  $(AMRLIB)/griddomshrink.f \
  $(AMRLIB)/domgrid.f \
  $(AMRLIB)/addflags.f \
  $(AMRLIB)/baseCheck.f \
  $(AMRLIB)/drivesort.f \
  $(AMRLIB)/flagcheck.f \
  $(AMRLIB)/auxcoarsen.f \
  $(AMRLIB)/intcopy.f \
  $(AMRLIB)/preintcopy.f \
  $(AMRLIB)/icall.f \
  $(AMRLIB)/preicall.f \
  $(AMRLIB)/inlinelimiter.f \
  $(AMRLIB)/cstore.f \
  $(AMRLIB)/saveqc.f \
  $(AMRLIB)/check.f \
  $(AMRLIB)/restrt.f \
  $(AMRLIB)/opendatafile.f \
  $(AMRLIB)/init_alloc.f90 \
  $(AMRLIB)/restrt_alloc.f90 \
  $(AMRLIB)/resize_alloc.f90

#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)

### DO NOT remove this line - make depends on it ###

