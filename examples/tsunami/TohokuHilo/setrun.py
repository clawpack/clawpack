"""
Module to set up run time parameters for Clawpack.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.

"""

import os
import numpy as np


#------------------------------
def setrun(claw_pkg='geoclaw'):
#------------------------------

    """
    Define the parameters used for running Clawpack.

    INPUT:
        claw_pkg expected to be "geoclaw" for this setrun.

    OUTPUT:
        rundata - object of class ClawRunData

    """

    from clawpack.clawutil import clawdata

    assert claw_pkg.lower() == 'geoclaw',  "Expected claw_pkg = 'geoclaw'"

    num_dim = 2
    rundata = clawdata.ClawRunData(claw_pkg, num_dim)

    #------------------------------------------------------------------
    # Problem-specific parameters to be written to setprob.data:
    #------------------------------------------------------------------
    
    #probdata = rundata.new_UserData(name='probdata',fname='setprob.data')

    #------------------------------------------------------------------
    # Standard Clawpack parameters to be written to claw.data:
    #   (or to amr2ez.data for AMR)
    #------------------------------------------------------------------

    clawdata = rundata.clawdata  # initialized when rundata instantiated


    # Set single grid parameters first.
    # See below for AMR parameters.


    # ---------------
    # Spatial domain:
    # ---------------

    # Number of space dimensions:
    clawdata.num_dim = num_dim

    # Lower and upper edge of computational domain:
    # Shift by dx/2 of finest bathy resolution (1/3 second) so
    # that cell centers align with bathy points
    dx2 = 0.5/(3*3600.)
    clawdata.lower[0] = 132. - dx2   # west longitude
    clawdata.upper[0] = 222. - dx2  # east longitude

    clawdata.lower[1] = 9. - dx2     # south latitude
    clawdata.upper[1] = 53. - dx2     # north latitude



    # Number of grid cells: Coarsest grid
    clawdata.num_cells[0] = 45
    clawdata.num_cells[1] = 22


    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    clawdata.num_eqn = 3

    # Number of auxiliary variables in the aux array (initialized in setaux)
    clawdata.num_aux = 4

    # Index of aux array corresponding to capacity function, if there is one:
    clawdata.capa_index = 2

    
    
    # -------------
    # Initial time:
    # -------------

    clawdata.t0 = 0.0


    # Restart from checkpoint file of a previous run?
    # Note: If restarting, you must also change the Makefile to set:
    #    RESTART = True
    # If restarting, t0 above should be from original run, and the
    # restart_file 'fort.chkNNNNN' specified below should be in 
    # the OUTDIR indicated in Makefile.

    clawdata.restart = False               # True to restart from prior results
    clawdata.restart_file = 'fort.chk00006'  # File to use for restart data

    # -------------
    # Output times:
    #--------------

    # Specify at what times the results should be written to fort.q files.
    # Note that the time integration stops after the final output time.
    # The solution at initial time t0 is always written in addition.

    clawdata.output_style = 2

    if clawdata.output_style==1:
        # Output nout frames at equally spaced times up to tfinal:
        clawdata.num_output_times = 13
        clawdata.tfinal = 13.*3600.0
        clawdata.output_t0 = True  # output at initial (or restart) time?

    elif clawdata.output_style == 2:
        # Specify a list of output times.
        from numpy import linspace
        times = []
        times = linspace(7.75,8.75,9) * 3600.
        clawdata.output_times =  list(times) #  + [10*3600.]

    elif clawdata.output_style == 3:
        # Output every iout timesteps with a total of ntot time steps:
        clawdata.output_step_interval = 2
        clawdata.total_steps = 6
        clawdata.output_t0 = True
        

    clawdata.output_format == 'ascii'      # 'ascii' or 'netcdf' 

    clawdata.output_q_components = 'all'   # could be list such as [True,True]
    clawdata.output_aux_components = 'none'  # could be list
    clawdata.output_aux_onlyonce = True    # output aux arrays only at t0



    # ---------------------------------------------------
    # Verbosity of messages to screen during integration:
    # ---------------------------------------------------

    # The current t, dt, and cfl will be printed every time step
    # at AMR levels <= verbosity.  Set verbosity = 0 for no printing.
    #   (E.g. verbosity == 2 means print only on levels 1 and 2.)
    clawdata.verbosity = 3



    # --------------
    # Time stepping:
    # --------------

    # if dt_variable==1: variable time steps used based on cfl_desired,
    # if dt_variable==0: fixed time steps dt = dt_initial will always be used.
    clawdata.dt_variable = True

    # Initial time step for variable dt.
    # If dt_variable==0 then dt=dt_initial for all steps:
    clawdata.dt_initial = 0.016

    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = 1e+99

    # Desired Courant number if variable dt used, and max to allow without
    # retaking step with a smaller dt:
    clawdata.cfl_desired = 0.75
    clawdata.cfl_max = 1.0

    # Maximum number of time steps to allow between output times:
    clawdata.steps_max = 50000




    # ------------------
    # Method to be used:
    # ------------------

    # Order of accuracy:  1 => Godunov,  2 => Lax-Wendroff plus limiters
    clawdata.order = 2
    
    # Use dimensional splitting? (not yet available for AMR)
    clawdata.dimensional_split = 'unsplit'
    
    # For unsplit method, transverse_waves can be 
    #  0 or 'none'      ==> donor cell (only normal solver used)
    #  1 or 'increment' ==> corner transport of waves
    #  2 or 'all'       ==> corner transport of 2nd order corrections too
    clawdata.transverse_waves = 2

    # Number of waves in the Riemann solution:
    clawdata.num_waves = 3
    
    # List of limiters to use for each wave family:  
    # Required:  len(limiter) == num_waves
    # Some options:
    #   0 or 'none'     ==> no limiter (Lax-Wendroff)
    #   1 or 'minmod'   ==> minmod
    #   2 or 'superbee' ==> superbee
    #   3 or 'mc'       ==> MC limiter
    #   4 or 'vanleer'  ==> van Leer
    clawdata.limiter = ['mc', 'mc', 'mc']

    clawdata.fwave = True    # True ==> use f-wave version of algorithms
    
    # Source terms splitting:
    #   src_split == 0 or 'none'    ==> no source term (src routine never called)
    #   src_split == 1 or 'godunov' ==> Godunov (1st order) splitting used, 
    #   src_split == 2 or 'strang'  ==> Strang (2nd order) splitting used,  not recommended.
    clawdata.source_split = 'godunov'


    # --------------------
    # Boundary conditions:
    # --------------------

    # Number of ghost cells (usually 2)
    clawdata.num_ghost = 2

    # Choice of BCs at xlower and xupper:
    #   0 => user specified (must modify bcN.f to use this option)
    #   1 => extrapolation (non-reflecting outflow)
    #   2 => periodic (must specify this at both boundaries)
    #   3 => solid wall for systems where q(2) is normal velocity

    clawdata.bc_lower[0] = 'extrap'
    clawdata.bc_upper[0] = 'extrap'

    clawdata.bc_lower[1] = 'extrap'
    clawdata.bc_upper[1] = 'extrap'


    # ---------------
    # AMR parameters:
    # ---------------


    # max number of refinement levels:
    clawdata.amr_levels_max = 4

    # List of refinement ratios at each level (length at least mxnest-1)
    #clawdata.inratx = [3,4,5,8]    #  was 4,30] for 6 levels
    clawdata.refinement_ratios_x = [3,4,5,8,45]  # 5: 15sec, 6: 1/3 sec
    clawdata.refinement_ratios_y = [3,4,5,8,45]
    clawdata.refinement_ratios_t = [1,1,1,1,1,1]


    # Specify type of each aux variable in clawdata.auxtype.
    # This must be a list of length maux, each element of which is one of:
    #   'center',  'capacity', 'xleft', or 'yleft'  (see documentation).

    clawdata.aux_type = ['center','capacity','yleft','center']


    # Flag using refinement routine flag2refine rather than richardson error
    clawdata.flag_richardson = False    # use Richardson?
    clawdata.flag2refine = True

    # steps to take on each level L between regriddings of level L+1:
    clawdata.regrid_interval = 3

    # width of buffer zone around flagged points:
    # (typically the same as regrid_interval so waves don't escape):
    clawdata.regrid_buffer_width  = 3

    # clustering alg. cutoff for (# flagged pts) / (total # of cells refined)
    # (closer to 1.0 => more small grids may be needed to cover flagged cells)
    clawdata.clustering_cutoff = 0.700000

    # print info about each regridding up to this level:
    clawdata.verbosity_regrid = 0  

    # Specify when checkpoint files should be created that can be
    # used to restart a computation.

    clawdata.checkpt_style = 0

    if clawdata.checkpt_style == 0:
        # Do not checkpoint at all
        pass

    elif clawdata.checkpt_style == 1:
        # Checkpoint only at tfinal.
        pass

    elif clawdata.checkpt_style == 2:
        # Specify a list of checkpoint times.  
        clawdata.checkpt_times = [0.1,0.15]

    elif clawdata.checkpt_style == 3:
        # Checkpoint every checkpt_interval timesteps (on Level 1)
        # and at the final time.
        clawdata.checkpt_interval = 5


    #  ----- For developers ----- 
    # Toggle debugging print statements:
    clawdata.dprint = False      # print domain flags
    clawdata.eprint = False      # print err est flags
    clawdata.edebug = False      # even more err est flags
    clawdata.gprint = False      # grid bisection/clustering
    clawdata.nprint = False      # proper nesting output
    clawdata.pprint = False      # proj. of tagged points
    clawdata.rprint = False      # print regridding summary
    clawdata.sprint = False      # space/memory output
    clawdata.tprint = False      # time step reporting each level
    clawdata.uprint = False      # update/upbnd reporting
    
    # More AMR parameters can be set -- see the defaults in pyclaw/data.py

    #------------------------------------------------------------------
    # GeoClaw specific parameters:
    #------------------------------------------------------------------

    rundata = setgeo(rundata)   # Defined below

    return rundata
    # end of function setrun
    # ----------------------


#-------------------
def setgeo(rundata):
#-------------------
    """
    Set GeoClaw specific runtime parameters.
    For documentation see ....
    """

    try:
	geodata = rundata.geodata
    except:
        print "*** Error, this rundata has no geodata attribute"
	raise AttributeError("Missing geodata attribute")

    geodata.variable_dt_refinement_ratios = True
       
    # == Physics ==
    geodata.gravity = 9.81
    geodata.coordinate_system = 2
    geodata.earth_radius = 6367.5e3

    # == Forcing Options
    geodata.coriolis_forcing = True

    # == Algorithm and Initial Conditions ==
    geodata.eta_init = 0.0
    geodata.dry_tolerance = 1.e-3
    geodata.wave_tolerance = 2.e-2  # 4.e-3 ###Bre changed this number
    geodata.deep_depth = 2.e2 
    geodata.max_level_deep = 4 
    geodata.friction_forcing = 1
    geodata.manning_coefficient =.025
    geodata.friction_depth = 500.

    # == settopo.data values ==
    geodata.topofiles = []
    # for topography, append lines of the form
    #    [topotype, minlevel, maxlevel, t1, t2, fname]

    topodir =   './HAI/'
    geodata.topofiles.append([3, 1, 1, 0., 1.e10, topodir+'1minlgr.asc'])
    geodata.topofiles.append([3, 1, 1, 0., 1.e10, topodir+'hawaii_6s.txt'])
    geodata.topofiles.append([3, 1, 1, 0., 1.e10, topodir+'hilo_3s.asc'])
    geodata.topofiles.append([3, 6, 6, 7.5*3600., 1.e10, topodir+'hilo_port_1_3s.asc'])

    # == setdtopo.data values ==
    geodata.dtopofiles = []
    # for moving topography, append lines of the form :   (<= 1 allowed for now!)
    #   [topotype, minlevel,maxlevel,fname]
    dtopodir = topodir
    fname = dtopodir + '/Fujii.txydz'  # name of earthquake deformation file
    geodata.dtopofiles.append([1,3,3,fname])


    # == setqinit.data values ==
    geodata.qinit_type = 0
    geodata.qinitfiles = []
    # for qinit perturbations, append lines of the form: (<= 1 allowed for now!)
    #   [minlev, maxlev, fname]

    # == setregions.data values ==
    geodata.regions = []
    # to specify regions of refinement append lines of the form
    #  [minlevel,maxlevel,t1,t2,x1,x2,y1,y2]
    geodata.regions.append([1, 2, 0., 1e9, 0, 360, -90, 90])
    geodata.regions.append([1, 3, 0., 5.*3600., 132., 220., 5., 40.])
    geodata.regions.append([1, 3, 5.*3600.,  8.*3600., 180., 220., 5., 40.])
    geodata.regions.append([4, 4, 7.*3600., 1e9, 204,205.5,19.4,20.4])
    geodata.regions.append([5, 5, 7.3*3600., 1e9, 204.85, 205, 19.68, 19.85])
    geodata.regions.append([6, 6, 7.5*3600., 1e9, 204.905,204.95,19.72,19.745])
    
    # == setgauges.data values ==
    rundata.gaugedata.gauges = []
    # for gauges append lines of the form  [gaugeno, x, y, t1, t2]
    # Velocity gauges:
    rundata.gaugedata.gauges.append([1125, 204.91802, 19.74517, 7.0*3600., 1.e9]) #Hilo
    rundata.gaugedata.gauges.append([1126, 204.93003, 19.74167, 7.0*3600., 1.e9]) #Hilo
    # geodata.gauges.append([11261, 204.93003, 19.739, 7.0*3600., 1.e9]) #Hilo
    # Tide gauge:
    rundata.gaugedata.gauges.append([7760, 204.9437, 19.7306,  7.0*3600., 1.e9]) #Hilo

    # == setfixedgrids.data values ==
    geodata.fixedgrids = []
    # for fixed grids append lines of the form
    # [t1,t2,noutput,x1,x2,y1,y2,xpoints,ypoints,\
    #  ioutarrivaltimes,ioutsurfacemax]
    
    # == Multilayer ==
    geodata.num_layers = 1
    geodata.rho = 1025.0
    geodata.richardson_tolerance = 0.95

    return rundata
    # end of function setgeo
    # ----------------------



if __name__ == '__main__':
    # Set up run-time parameters and write all data files.
    import sys
    rundata = setrun(*sys.argv[1:])
    rundata.write()

