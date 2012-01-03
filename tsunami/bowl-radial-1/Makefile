
# Makefile for Clawpack code in this directory.
# This version only sets the local files and frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
CLAWMAKE = $(CLAW)/clawutil/src/Makefile-4.x

# See the above file for details and a list of make options, or type
#   make .help
# at the unix prompt.


# Adjust these variables if desired:
# ----------------------------------

CLAW_PKG = geoclaw                  # Clawpack package to use
CLAW_EXE = xgeoclaw                 # Executable to create
CLAW_setrun_file = setrun.py        # File containing function to make data
CLAW_OUTDIR = _output               # Directory for output
CLAW_setplot_file = setplot.py      # File containing function to set plots
CLAW_PLOTDIR = _plots               # Directory for plots

# Environment variable FC should be set to fortran compiler, e.g. gfortran
FC ?= gfortran   # default if not set as environment variable
# Add any desired compiler flags such as -g here:
FFLAGS ?=  -O2 

# If you want to include NetCDF output add these flags
# FFLAGS ?=-I/usr/include  -lnetcdf -lnetcdff


# ---------------------------------
# List of sources for this program:
# ---------------------------------


CLAW_SOURCES = \

# Clawpack libraries to be used:
CLAW_LIB = $(CLAW)/amrclaw/2d/lib
GEOLIB = $(CLAW)/geoclaw/2d/lib

CLAW_LIBSOURCES = \
  $(GEOLIB)/geoclaw_mod.f90 \
  $(GEOLIB)/topo_mod.f90 \
  $(GEOLIB)/dtopo_mod.f90 \
  $(GEOLIB)/setprob_geo.f \
  $(GEOLIB)/qinit_geo.f \
  $(GEOLIB)/readqinit_geo.f \
  $(GEOLIB)/movetopo_geo.f \
  $(GEOLIB)/cellgridintegrate_geo.f \
  $(GEOLIB)/topointegral_geo.f \
  $(GEOLIB)/bilinearintegral_geo.f \
  $(GEOLIB)/stepgrid_geo.f \
  $(GEOLIB)/fgridinterp_geo.f \
  $(GEOLIB)/fgridout_geo.f \
  $(GEOLIB)/rpn2ez_geo.f \
  $(GEOLIB)/riemannsolvers_geo.f \
  $(GEOLIB)/rpt2_geo.f \
  $(GEOLIB)/src2_geo.f \
  $(GEOLIB)/src1d_geo.f \
  $(GEOLIB)/step2_geo.f \
  $(GEOLIB)/flux2fw_geo.f \
  $(GEOLIB)/qad_geo.f \
  $(GEOLIB)/valout_geo.f \
  $(GEOLIB)/filval_geo.f \
  $(GEOLIB)/filpatch_geo.f \
  $(GEOLIB)/bc2amr_geo.f \
  $(GEOLIB)/update_geo.f \
  $(GEOLIB)/amr2ez_geo.f \
  $(GEOLIB)/setaux_geo.f \
  $(GEOLIB)/setqinit_geo.f \
  $(GEOLIB)/setregions_geo.f \
  $(GEOLIB)/setgauges_geo.f \
  $(GEOLIB)/setfixedgrids_geo.f \
  $(GEOLIB)/b4step2_geo.f \
  $(GEOLIB)/flag2refine_geo.f  \
  $(GEOLIB)/allowflag_geo.f  \
  $(GEOLIB)/dumpgauge_geo.f \
  $(GEOLIB)/upbnd_geo.f  \
  $(GEOLIB)/tick_geo.f \
  $(GEOLIB)/setgrd_geo.f \
  $(GEOLIB)/gfixup_geo.f \
  $(GEOLIB)/ginit_geo.f \
  $(GEOLIB)/getmaxspeed_geo.f \
  $(CLAW_LIB)/prefilp.f \
  $(CLAW_LIB)/quick_sort1.f \
  $(CLAW_LIB)/advanc.f \
  $(CLAW_LIB)/estdt.f \
  $(CLAW_LIB)/bound.f \
  $(CLAW_LIB)/check4nans.f \
  $(CLAW_LIB)/spest.f  \
  $(CLAW_LIB)/igetsp.f \
  $(CLAW_LIB)/reclam.f \
  $(CLAW_LIB)/birect.f \
  $(CLAW_LIB)/cleanup.f \
  $(CLAW_LIB)/colate.f \
  $(CLAW_LIB)/errest.f \
  $(CLAW_LIB)/bufnst.f \
  $(CLAW_LIB)/errf1.f \
  $(CLAW_LIB)/fixcapaq.f \
  $(CLAW_LIB)/flglvl.f \
  $(CLAW_LIB)/fluxad.f \
  $(CLAW_LIB)/fluxsv.f \
  $(CLAW_LIB)/grdfit.f \
  $(CLAW_LIB)/intfil.f \
  $(CLAW_LIB)/moment.f \
  $(CLAW_LIB)/nestck.f \
  $(CLAW_LIB)/prepc.f \
  $(CLAW_LIB)/prepf.f \
  $(CLAW_LIB)/projec.f \
  $(CLAW_LIB)/signs.f \
  $(CLAW_LIB)/findcut.f \
  $(CLAW_LIB)/smartbis.f \
  $(CLAW_LIB)/putnod.f \
  $(CLAW_LIB)/putsp.f \
  $(CLAW_LIB)/regrid.f \
  $(CLAW_LIB)/setuse.f \
  $(CLAW_LIB)/stst1.f \
  $(CLAW_LIB)/trimbd.f \
  $(CLAW_LIB)/nodget.f \
  $(CLAW_LIB)/basic.f \
  $(CLAW_LIB)/outval.f \
  $(CLAW_LIB)/copysol.f \
  $(CLAW_LIB)/outvar.f \
  $(CLAW_LIB)/outmsh.f \
  $(CLAW_LIB)/outtre.f \
  $(CLAW_LIB)/domain.f  \
  $(CLAW_LIB)/cellave.f \
  $(CLAW_LIB)/fdisc.f \
  $(CLAW_LIB)/fss.f \
  $(CLAW_LIB)/zeroin.f \
  $(CLAW_LIB)/setflags.f \
  $(CLAW_LIB)/shiftset.f \
  $(CLAW_LIB)/conck.f \
  $(CLAW_LIB)/domshrink.f \
  $(CLAW_LIB)/domprep.f \
  $(CLAW_LIB)/domup.f \
  $(CLAW_LIB)/domcopy.f \
  $(CLAW_LIB)/coarsen.f \
  $(CLAW_LIB)/auxcoarsen.f \
  $(CLAW_LIB)/intcopy.f \
  $(CLAW_LIB)/preintcopy.f \
  $(CLAW_LIB)/icall.f \
  $(CLAW_LIB)/preicall.f \
  $(CLAW_LIB)/inlinelimiter.f \
  $(CLAW_LIB)/cstore.f \
  $(CLAW_LIB)/saveqc.f \
  $(CLAW_LIB)/check.f \
  $(CLAW_LIB)/restrt.f \
  $(CLAW_LIB)/opendatafile.f \
  $(CLAW_LIB)/init_alloc.f90 \
  $(CLAW_LIB)/restrt_alloc.f90 \
  $(CLAW_LIB)/resize_alloc.f90 \
# $(CLAW_LIB)/resize_alloc_static.f90 \



#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)

topo: .data
	python maketopo.py


### DO NOT remove this line - make depends on it ###

