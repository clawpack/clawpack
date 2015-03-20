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
EXE = xgeoclaw                 # Executable to create
SETRUN_FILE = setrun.py        # File containing function to make data
OUTDIR = _output               # Directory for output
SETPLOT_FILE = setplot.py      # File containing function to set plots
PLOTDIR = _plots               # Directory for plots

RESTART = False

# Environment variable FC should be set to fortran compiler, e.g. gfortran

# Compiler flags can be specified here or set as an environment variable
FFLAGS ?=

# Include flags for LAPACK linking
UNAME = $(shell uname)
ifeq ($(UNAME), Darwin)
	ifndef LFLAGS
		LFLAGS = -framework accelerate $(FFLAGS)
	else
		LFLAGS += -framework accelerate $(FFLAGS)
	endif
else ifeq ($(shell hostname), sandy.ices.utexas.edu)
	# Created using Intel MKL Link Line Advisor
	# http://software.intel.com/en-us/articles/intel-mkl-link-line-advisor/
	FFLAGS += -I$(MKLROOT)/include
	ifndef LFLAGS                                                                                                                                    
	    LFLAGS =  -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm $(FFLAGS)
	else
	    LFLAGS +=  -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_ilp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm $(FFLAGS)
	endif
else ifeq ($(UNAME), Linux)
	ifndef LFLAGS
		LFLAGS = -llapack $(FFLAGS)
	else
		LFLAGS += -llapack $(FFLAGS)
	endif
endif

# ---------------------------------
# List of sources for this program:
# ---------------------------------

AMRLIB = $(CLAW)/amrclaw/src/2d
GEOLIB = $(CLAW)/geoclaw/src/2d/shallow

MODULES = \
 $(GEOLIB)/utility_module.f90 \
 $(AMRLIB)/amr_module.f90 \
 $(GEOLIB)/geoclaw_module.f90 \
 $(AMRLIB)/regions_module.f90 \
 $(AMRLIB)/gauges_module.f90 \
 $(GEOLIB)/topo_module.f90 \
 $(GEOLIB)/refinement_module.f90 \
 $(GEOLIB)/fixedgrids_module.f90 \
 $(GEOLIB)/fgmax_module.f90 \
 $(GEOLIB)/holland_storm_module.f90 \
 $(GEOLIB)/constant_storm_module.f90 \
 $(GEOLIB)/stommel_storm_module.f90 \
 $(GEOLIB)/storm_module.f90 \
 $(GEOLIB)/friction_module.f90 \
 $(GEOLIB)/multilayer/multilayer_module.f90 \
 ./qinit_module.f90

SOURCES = \
  $(GEOLIB)/multilayer/setprob.f90 \
  ./qinit.f90 \
  $(GEOLIB)/multilayer/setaux.f90 \
  $(GEOLIB)/multilayer/b4step2.f90 \
  $(GEOLIB)/multilayer/dumpgauge.f \
  $(GEOLIB)/topo_update.f90 \
  $(GEOLIB)/multilayer/getmaxspeed.f90 \
  $(GEOLIB)/multilayer/stepgrid.f \
  $(GEOLIB)/multilayer/valout.f \
  $(GEOLIB)/multilayer/src2.f90 \
  $(GEOLIB)/movetopo.f \
  $(GEOLIB)/cellgridintegrate2.f \
  $(GEOLIB)/topointegral.f \
  $(GEOLIB)/bilinearintegral.f \
  $(CLAW)/riemann/src/rpn2_layered_shallow_water.f90 \
  $(CLAW)/riemann/src/rpt2_layered_shallow_water.f90 \
  $(CLAW)/riemann/src/geoclaw_riemann_utils.f \
  $(GEOLIB)/src1d.f90 \
  $(GEOLIB)/step2.f90 \
  $(GEOLIB)/flux2fw.f \
  $(GEOLIB)/qad.f \
  $(GEOLIB)/filval.f90 \
  $(GEOLIB)/filpatch.f90 \
  $(AMRLIB)/prefilp.f90 \
  $(AMRLIB)/trimbd.f \
  $(AMRLIB)/bound.f \
  $(AMRLIB)/intfil.f \
  $(GEOLIB)/bc2amr.f \
  $(GEOLIB)/update.f \
  $(GEOLIB)/flag2refine2.f90  \
  $(AMRLIB)/flagregions2.f90  \
  $(GEOLIB)/allowflag.f90  \
  $(GEOLIB)/upbnd.f  \
  $(GEOLIB)/tick.f \
  $(GEOLIB)/setgrd.f \
  $(GEOLIB)/gfixup.f \
  $(GEOLIB)/ginit.f \
  $(GEOLIB)/advanc.f \
  $(GEOLIB)/amr2.f90 \
  $(GEOLIB)/fgmax_read.f90 \
  $(GEOLIB)/fgmax_frompatch.f90 \
  $(GEOLIB)/fgmax_interpolate.f90 \
  $(GEOLIB)/fgmax_values.f90 \
  $(GEOLIB)/fgmax_finalize.f90 \
  $(AMRLIB)/quick_sort1.f \
  $(AMRLIB)/estdt.f \
  $(AMRLIB)/check4nans.f90 \
  $(AMRLIB)/spest2.f  \
  $(AMRLIB)/init_iflags.f \
  $(AMRLIB)/igetsp.f \
  $(AMRLIB)/reclam.f \
  $(AMRLIB)/birect.f \
  $(AMRLIB)/cleanup.f \
  $(AMRLIB)/colate2.f \
  $(AMRLIB)/errest.f \
  $(AMRLIB)/prepregstep.f \
  $(AMRLIB)/prepbigstep.f \
  $(AMRLIB)/bufnst2.f \
  $(AMRLIB)/flagger.f \
  $(AMRLIB)/flagregions.f \
  $(AMRLIB)/errf1.f \
  $(AMRLIB)/fixcapaq.f \
  $(AMRLIB)/flglvl2.f \
  $(AMRLIB)/fluxad.f \
  $(AMRLIB)/fluxsv.f \
  $(AMRLIB)/grdfit2.f \
  $(AMRLIB)/moment.f \
  $(AMRLIB)/nestck2.f \
  $(AMRLIB)/prepc.f \
  $(AMRLIB)/prepf.f \
  $(AMRLIB)/projec2.f \
  $(AMRLIB)/signs.f \
  $(AMRLIB)/findcut.f \
  $(AMRLIB)/smartbis.f \
  $(AMRLIB)/putnod.f \
  $(AMRLIB)/putsp.f \
  $(AMRLIB)/regrid.f \
  $(AMRLIB)/setuse.f \
  $(AMRLIB)/stst1.f \
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
  $(AMRLIB)/shiftset2.f \
  $(AMRLIB)/conck.f \
  $(AMRLIB)/domshrink.f \
  $(AMRLIB)/domprep.f \
  $(AMRLIB)/domup.f \
  $(AMRLIB)/domcopy.f \
  $(AMRLIB)/setdomflags.f \
  $(AMRLIB)/setIndices.f \
  $(AMRLIB)/coarseGridFlagSet.f \
  $(AMRLIB)/addflags.f \
  $(AMRLIB)/baseCheck.f \
  $(AMRLIB)/drivesort.f \
  $(AMRLIB)/flagcheck.f \
  $(AMRLIB)/domgrid.f \
  $(AMRLIB)/setPhysBndryFlags.f \
  $(AMRLIB)/griddomup.f \
  $(AMRLIB)/griddomcopy.f \
  $(AMRLIB)/griddomshrink.f \
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
  $(AMRLIB)/resize_alloc.f90

#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)
