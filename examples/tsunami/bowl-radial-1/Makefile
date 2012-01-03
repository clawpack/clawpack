
# Makefile for Clawpack code in this directory.
# This version only sets the local files and frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
CLAWMAKE = $(CLAW)/clawutil/src/Makefile.common

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
FFLAGS ?=  -g

# If you want to include NetCDF output add these flags
# FFLAGS ?=-I/usr/include  -lnetcdf -lnetcdff


# ---------------------------------
# List of sources for this program:
# ---------------------------------

AMRLIB = $(CLAW)/amrclaw/src/2d
GEOLIB = $(CLAW)/geoclaw/src/2d

MODULES = \
 $(AMRLIB)/amr_module.f \
 $(GEOLIB)/geoclaw_module.f \
 $(GEOLIB)/topo_module.f \
 $(GEOLIB)/dtopo_module.f \
 $(GEOLIB)/qinit_module.f \
 $(GEOLIB)/gauges_module.f \
 $(GEOLIB)/regions_module.f \
 $(GEOLIB)/fixedgrids_module.f \

SOURCES = \
  $(GEOLIB)/setprob.f \
  $(GEOLIB)/qinit.f \
  $(GEOLIB)/readqinit.f \
  $(GEOLIB)/movetopo.f \
  $(GEOLIB)/cellgridintegrate.f \
  $(GEOLIB)/topointegral.f \
  $(GEOLIB)/bilinearintegral.f \
  $(GEOLIB)/stepgrid.f \
  $(GEOLIB)/fgridinterp.f \
  $(GEOLIB)/fgridout.f \
  $(GEOLIB)/rpn2ez.f \
  $(GEOLIB)/riemannsolvers.f \
  $(GEOLIB)/rpt2.f \
  $(GEOLIB)/src2.f \
  $(GEOLIB)/src1d.f \
  $(GEOLIB)/step2.f \
  $(GEOLIB)/flux2fw.f \
  $(GEOLIB)/qad.f \
  $(GEOLIB)/valout.f \
  $(GEOLIB)/filval.f \
  $(GEOLIB)/filpatch.f \
  $(GEOLIB)/bc2amr.f \
  $(GEOLIB)/update.f \
  $(GEOLIB)/amr2ez.f \
  $(GEOLIB)/setaux.f \
  $(GEOLIB)/setqinit.f \
  $(GEOLIB)/setregions.f \
  $(GEOLIB)/setgauges.f \
  $(GEOLIB)/setfixedgrids.f \
  $(GEOLIB)/b4step2.f \
  $(GEOLIB)/flag2refine.f  \
  $(GEOLIB)/allowflag.f  \
  $(GEOLIB)/dumpgauge.f \
  $(GEOLIB)/upbnd.f  \
  $(GEOLIB)/tick.f \
  $(GEOLIB)/setgrd.f \
  $(GEOLIB)/gfixup.f \
  $(GEOLIB)/ginit.f \
  $(GEOLIB)/getmaxspeed.f \
  $(AMRLIB)/prefilp.f \
  $(AMRLIB)/quick_sort1.f \
  $(AMRLIB)/advanc.f \
  $(AMRLIB)/estdt.f \
  $(AMRLIB)/bound.f \
  $(AMRLIB)/check4nans.f \
  $(AMRLIB)/spest.f  \
  $(AMRLIB)/igetsp.f \
  $(AMRLIB)/reclam.f \
  $(AMRLIB)/birect.f \
  $(AMRLIB)/cleanup.f \
  $(AMRLIB)/colate.f \
  $(AMRLIB)/errest.f \
  $(AMRLIB)/bufnst.f \
  $(AMRLIB)/errf1.f \
  $(AMRLIB)/fixcapaq.f \
  $(AMRLIB)/flglvl.f \
  $(AMRLIB)/fluxad.f \
  $(AMRLIB)/fluxsv.f \
  $(AMRLIB)/grdfit.f \
  $(AMRLIB)/intfil.f \
  $(AMRLIB)/moment.f \
  $(AMRLIB)/nestck.f \
  $(AMRLIB)/prepc.f \
  $(AMRLIB)/prepf.f \
  $(AMRLIB)/projec.f \
  $(AMRLIB)/signs.f \
  $(AMRLIB)/findcut.f \
  $(AMRLIB)/smartbis.f \
  $(AMRLIB)/putnod.f \
  $(AMRLIB)/putsp.f \
  $(AMRLIB)/regrid.f \
  $(AMRLIB)/setuse.f \
  $(AMRLIB)/stst1.f \
  $(AMRLIB)/trimbd.f \
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
  $(AMRLIB)/conck.f \
  $(AMRLIB)/domshrink.f \
  $(AMRLIB)/domprep.f \
  $(AMRLIB)/domup.f \
  $(AMRLIB)/domcopy.f \
  $(AMRLIB)/coarsen.f \
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
  $(AMRLIB)/resize_alloc.f90 \
# $(AMRLIB)/resize_alloc_static.f90 \



#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)

topo: .data
	python maketopo.py


### DO NOT remove this line - make depends on it ###

