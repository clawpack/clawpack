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
else ifeq ($(UNAME), Linux)
	ifeq ($(FC), ifort)
		# Created using Intel MKL Link Line Advisor
		# Please refer to 
		# http://software.intel.com/en-us/articles/intel-mkl-link-line-advisor/
		# for other installation types
		FFLAGS += -I$(MKLROOT)/include
		ifndef LFLAGS                                                                                                                                    
		    LFLAGS =  -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm $(FFLAGS)
		else
		    LFLAGS +=  -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_ilp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm $(FFLAGS)
		endif
	else ifeq ($(FC), gfortran)
		ifndef LFLAGS
			LFLAGS = -llapack $(FFLAGS)
		else
			LFLAGS += -llapack $(FFLAGS)
		endif
	endif
endif

# ---------------------------------
# package sources for this program:
# ---------------------------------

GEOLIB = $(CLAW)/geoclaw/src/2d/shallow
include $(GEOLIB)/Makefile.geoclaw

# ---------------------------------------
# package sources specifically to exclude
# (i.e. if a custom replacement source 
#  under a different name is provided)
# ---------------------------------------

EXCLUDE_MODULES = \

EXCLUDE_SOURCES = \

# ----------------------------------------
# List of custom sources for this program:
# ----------------------------------------

MODULES = \
 $(GEOLIB)/multilayer/gauges_module.f90 \
 ./qinit_module.f90

SOURCES = \
  ./qinit.f90 \
  $(GEOLIB)/multilayer/setaux.f90 \
  $(GEOLIB)/multilayer/b4step2.f90 \
  $(GEOLIB)/multilayer/getmaxspeed.f90 \
  $(GEOLIB)/multilayer/stepgrid.f \
  $(GEOLIB)/multilayer/valout.f \
  $(GEOLIB)/multilayer/src2.f90 \
  $(CLAW)/riemann/src/rpn2_layered_shallow_water.f90 \
  $(CLAW)/riemann/src/rpt2_layered_shallow_water.f90 \
  $(CLAW)/riemann/src/geoclaw_riemann_utils.f

#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)
