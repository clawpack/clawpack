"""
Module to set up run time parameters for Clawpack.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.

"""

from __future__ import absolute_import
from __future__ import print_function
import numpy as numpy

import clawpack.geoclaw.data
import clawpack.geoclaw.topotools as tt


# Rotation transformations
def transform_c2p(x,y,x0,y0,theta):
    return ((x+x0)*numpy.cos(theta) - (y+y0)*numpy.sin(theta),
            (x+x0)*numpy.sin(theta) + (y+y0)*numpy.cos(theta))

def transform_p2c(x,y,x0,y0,theta):
    return ( x*numpy.cos(theta) + y*numpy.sin(theta) - x0,
            -x*numpy.sin(theta) + y*numpy.cos(theta) - y0)


# Class containing some setup for the qinit especially for multilayer tests 
class QinitMultilayerData(clawpack.geoclaw.data.QinitData):
    r"""
    Modified Qinit data object for multiple layers

    """

    def __init__(self):
        super(QinitMultilayerData, self).__init__()

        # Test qinit data > 4
        self.add_attribute("init_location", [0.0, 0.0])
        self.add_attribute("wave_family", 1)
        self.add_attribute("epsilon", 0.02)
        self.add_attribute("angle", 0.0)
        self.add_attribute("sigma", 0.02)

    def write(self, data_source='setrun.py'):

        # Initial perturbation
        self.open_data_file('qinit.data',data_source)
        self.data_write('qinit_type')

        # Perturbation requested
        if self.qinit_type == 0:
            pass
        elif 0 < self.qinit_type < 5:
            # Check to see if each qinit file is present and then write the data
            for tfile in self.qinitfiles:
                try:
                    fname = "'%s'" % os.path.abspath(tfile[-1])
                except:
                    raise Warning("File %s was not found." % fname)
                    # raise MissingFile("file not found")
                self._out_file.write("\n%s  \n" % fname)
                self._out_file.write("%3i %3i \n" % tuple(tfile[:-1]))
        elif self.qinit_type >= 5 and self.qinit_type <= 9:
            self.data_write('epsilon')
            self.data_write("init_location")
            self.data_write("wave_family")
            self.data_write("angle")
            self.data_write("sigma")
        else:
            raise ValueError("Invalid qinit_type parameter %s." % self.qinit_type)
        self.close_data_file()



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

    from clawpack.clawutil import data

    assert claw_pkg.lower() == 'geoclaw',  "Expected claw_pkg = 'geoclaw'"

    num_dim = 2
    rundata = data.ClawRunData(claw_pkg, num_dim)

    #------------------------------------------------------------------
    # GeoClaw specific parameters:
    #------------------------------------------------------------------
    rundata = setgeo(rundata)
    rundata = set_multilayer(rundata)

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
    clawdata.lower[0] = -1      # west longitude
    clawdata.upper[0] = 2.0       # east longitude

    clawdata.lower[1] = -1.0       # south latitude
    clawdata.upper[1] = 2.0         # north latitude



    # Number of grid cells: Coarsest grid
    clawdata.num_cells[0] = 150
    clawdata.num_cells[1] = 150

    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    clawdata.num_eqn = 6

    # Number of auxiliary variables in the aux array (initialized in setaux)
    clawdata.num_aux = 4 + rundata.multilayer_data.num_layers

    # Index of aux array corresponding to capacity function, if there is one:
    clawdata.capa_index = 0

    
    
    # -------------
    # Initial time:
    # -------------

    clawdata.t0 = 0.0


    # Restart from checkpoint file of a previous run?
    # If restarting, t0 above should be from original run, and the
    # restart_file 'fort.chkNNNNN' specified below should be in 
    # the OUTDIR indicated in Makefile.

    clawdata.restart = False               # True to restart from prior results
    clawdata.restart_file = 'fort.chk00036'  # File to use for restart data

    # -------------
    # Output times:
    #--------------

    # Specify at what times the results should be written to fort.q files.
    # Note that the time integration stops after the final output time.
    # The solution at initial time t0 is always written in addition.

    clawdata.output_style = 1

    if clawdata.output_style==1:
        # Output nout frames at equally spaced times up to tfinal:
        clawdata.num_output_times = 40
        clawdata.tfinal = 1.0
        clawdata.output_t0 = True  # output at initial (or restart) time?

    elif clawdata.output_style == 2:
        # Specify a list of output times.
        clawdata.output_times = [0.5, 1.0]

    elif clawdata.output_style == 3:
        # Output every iout timesteps with a total of ntot time steps:
        clawdata.output_step_interval = 1
        clawdata.total_steps = 10
        clawdata.output_t0 = True
        

    clawdata.output_format = 'ascii'      # 'ascii' or 'binary' 

    clawdata.output_q_components = 'all'   # need all
    clawdata.output_aux_components = 'all'  # eta=h+B is in q
    clawdata.output_aux_onlyonce = False    # output aux arrays each frame



    # ---------------------------------------------------
    # Verbosity of messages to screen during integration:
    # ---------------------------------------------------

    # The current t, dt, and cfl will be printed every time step
    # at AMR levels <= verbosity.  Set verbosity = 0 for no printing.
    #   (E.g. verbosity == 2 means print only on levels 1 and 2.)
    clawdata.verbosity = 1



    # --------------
    # Time stepping:
    # --------------

    # if dt_variable==1: variable time steps used based on cfl_desired,
    # if dt_variable==0: fixed time steps dt = dt_initial will always be used.
    clawdata.dt_variable = True

    # Initial time step for variable dt.
    # If dt_variable==0 then dt=dt_initial for all steps:
    clawdata.dt_initial = 0.00225

    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = 1e+99

    # Desired Courant number if variable dt used, and max to allow without
    # retaking step with a smaller dt:
    clawdata.cfl_desired = 0.75
    clawdata.cfl_max = 1.0
    # clawdata.cfl_desired = 0.45
    # clawdata.cfl_max = 0.5

    # Maximum number of time steps to allow between output times:
    clawdata.steps_max = 5000




    # ------------------
    # Method to be used:
    # ------------------

    # Order of accuracy:  1 => Godunov,  2 => Lax-Wendroff plus limiters
    clawdata.order = 2
    
    # Use dimensional splitting? (not yet available for AMR)
    #  0 or 'unsplit' or none'  ==> Unsplit
    #  1 or 'increment'         ==> corner transport of waves
    #  2 or 'all'               ==> corner transport of 2nd order corrections too
    clawdata.dimensional_split = "unsplit"
    
    # For unsplit method, transverse_waves can be 
    #  0 or 'none'      ==> donor cell (only normal solver used)
    #  1 or 'increment' ==> corner transport of waves
    #  2 or 'all'       ==> corner transport of 2nd order corrections too
    clawdata.transverse_waves = 2

    # Number of waves in the Riemann solution:
    clawdata.num_waves = 6
    
    # List of limiters to use for each wave family:  
    # Required:  len(limiter) == num_waves
    # Some options:
    #   0 or 'none'     ==> no limiter (Lax-Wendroff)
    #   1 or 'minmod'   ==> minmod
    #   2 or 'superbee' ==> superbee
    #   3 or 'mc'       ==> MC limiter
    #   4 or 'vanleer'  ==> van Leer
    clawdata.limiter = ['mc', 'mc', 'mc', 'mc', 'mc', 'mc']

    clawdata.use_fwaves = True    # True ==> use f-wave version of algorithms
    
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



    # --------------
    # Checkpointing:
    # --------------

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


    # ---------------
    # AMR parameters:
    # ---------------
    amrdata = rundata.amrdata

    # max number of refinement levels:
    amrdata.amr_levels_max = 1

    # List of refinement ratios at each level (length at least mxnest-1)
    amrdata.refinement_ratios_x = [2,6]
    amrdata.refinement_ratios_y = [2,6]
    amrdata.refinement_ratios_t = [2,6]


    # Specify type of each aux variable in amrdata.auxtype.
    # This must be a list of length maux, each element of which is one of:
    #   'center',  'capacity', 'xleft', or 'yleft'  (see documentation).

    amrdata.aux_type = ['center','center','yleft','center','center','center']


    # Flag using refinement routine flag2refine rather than richardson error
    amrdata.flag_richardson = False    # use Richardson?
    amrdata.flag2refine = True

    # steps to take on each level L between regriddings of level L+1:
    amrdata.regrid_interval = 3

    # width of buffer zone around flagged points:
    # (typically the same as regrid_interval so waves don't escape):
    amrdata.regrid_buffer_width  = 2

    # clustering alg. cutoff for (# flagged pts) / (total # of cells refined)
    # (closer to 1.0 => more small grids may be needed to cover flagged cells)
    amrdata.clustering_cutoff = 0.700000

    # print info about each regridding up to this level:
    amrdata.verbosity_regrid = 0  

    #  ----- For developers ----- 
    # Toggle debugging print statements:
    amrdata.dprint = False      # print domain flags
    amrdata.eprint = False      # print err est flags
    amrdata.edebug = False      # even more err est flags
    amrdata.gprint = False      # grid bisection/clustering
    amrdata.nprint = False      # proper nesting output
    amrdata.pprint = False      # proj. of tagged points
    amrdata.rprint = False      # print regridding summary
    amrdata.sprint = False      # space/memory output
    amrdata.tprint = True       # time step reporting each level
    amrdata.uprint = False      # update/upbnd reporting
    
    # More AMR parameters can be set -- see the defaults in pyclaw/data.py

    # ---------------
    # Regions:
    # ---------------
    rundata.regiondata.regions = []
    # to specify regions of refinement append lines of the form
    #  [minlevel,maxlevel,t1,t2,x1,x2,y1,y2]

    # ---------------
    # Gauges:
    # ---------------
    rundata.gaugedata.gauges = []
    # for gauges append lines of the form  [gaugeno, x, y, t1, t2]
    gauge_locations = [-0.1,0.0,0.1,0.2,0.3]
    for (i,x_c) in enumerate(gauge_locations):
        # y0 = (self.run_data.clawdata.yupper - self.run_data.clawdata.ylower) / 2.0
        # x_p,y_p = transform_c2p(x_c,0.0,location[0],location[1],angle)
        x_p = x_c * numpy.cos(0.0)
        y_p = x_c * numpy.sin(0.0)
        # print "+=====+"
        # print x_c,0.0
        # print x_p,y_p
        if (rundata.clawdata.lower[0] < x_p < rundata.clawdata.upper[0] and
                rundata.clawdata.lower[1] < y_p < rundata.clawdata.upper[1]):
            rundata.gaugedata.gauges.append([i, x_p, y_p, 0.0, 1e10])
            # print "Gauge %s: (%s,%s)" % (i,x_p,y_p)
    # print "+=====+"

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
        geo_data = rundata.geo_data
    except:
        print("*** Error, this rundata has no geo_data attribute")
        raise AttributeError("Missing geo_data attribute")
       
    # == Physics ==
    geo_data.gravity = 9.81
    geo_data.coordinate_system = 1
    geo_data.earth_radius = 6367.5e3

    # == Forcing Options
    geo_data.coriolis_forcing = False

    # == Algorithm and Initial Conditions ==
    geo_data.sea_level = 0.0
    geo_data.dry_tolerance = 1.e-3
    geo_data.friction_forcing = True
    geo_data.manning_coefficient = 0.025
    geo_data.friction_depth = 1e6

    # Refinement settings
    refinement_data = rundata.refinement_data
    refinement_data.variable_dt_refinement_ratios = True
    refinement_data.wave_tolerance = 1.e-1
    refinement_data.deep_depth = 1e2
    refinement_data.max_level_deep = 3

    # == settopo.data values ==
    topo_data = rundata.topo_data
    # for topography, append lines of the form
    #    [topotype, minlevel, maxlevel, t1, t2, fname]
    topo_data.topofiles.append([2, 1, 5, 0.0, 1e10, 'topo.tt2'])
    
    # == setdtopo.data values ==
    dtopo_data = rundata.dtopo_data
    # for moving topography, append lines of the form :   (<= 1 allowed for now!)
    #   [topotype, minlevel,maxlevel,fname]


    return rundata
    # end of function setgeo
    # ----------------------


def set_multilayer(rundata):

    data = rundata.multilayer_data

    # Physics parameters
    data.num_layers = 2
    data.rho = [0.9,1.0]
    data.eta = [0.0,-0.6]
    
    # Algorithm parameters
    data.eigen_method = 2
    data.inundation_method = 2
    data.richardson_tolerance = 0.95
    # data.wave_tolerance = [0.1,0.1]
    # data.dry_limit = True

    rundata.replace_data('qinit_data', QinitMultilayerData())
    rundata.qinit_data.qinit_type = 6
    rundata.qinit_data.epsilon = 0.02
    rundata.qinit_data.angle = 0.0
    rundata.qinit_data.sigma = 0.02
    rundata.qinit_data.wave_family = 4
    rundata.qinit_data.init_location = [-0.1,0.0]

    return rundata


def bathy_step(x, y, location=0.15, angle=0.0, left=-1.0, right=-0.2):
    x_c,y_c = transform_p2c(x, y, location, 0.0, angle)
    return ((x_c <= 0.0) * left 
          + (x_c >  0.0) * right)


def write_topo_file(run_data, out_file, **kwargs):

    # Make topography
    topo_func = lambda x, y: bathy_step(x, y, **kwargs)
    topo = tt.Topography(topo_func=topo_func)
    topo.x = numpy.linspace(run_data.clawdata.lower[0], 
                            run_data.clawdata.upper[0], 
                            run_data.clawdata.num_cells[0] + 8)
    topo.y = numpy.linspace(run_data.clawdata.lower[1], 
                            run_data.clawdata.upper[1], 
                            run_data.clawdata.num_cells[1] + 8)
    topo.write(out_file)

    # Write out simple bathy geometry file for communication to the plotting
    with open("./bathy_geometry.data", 'w') as bathy_geometry_file:
        if "location" in kwargs:
            location = kwargs['location']
        else:
            location = 0.15
        if "angle" in kwargs:
            angle = kwargs['angle']
        else:
            angle = 0.0
        bathy_geometry_file.write("%s\n%s" % (location, angle) )


if __name__ == '__main__':
    # Set up run-time parameters and write all data files.
    import sys
    if len(sys.argv) == 2:
        rundata = setrun(sys.argv[1])
    else:
        rundata = setrun()

    rundata.write()

    write_topo_file(rundata, 'topo.tt2')
