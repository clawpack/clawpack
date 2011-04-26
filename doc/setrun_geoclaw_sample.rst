

.. _setrun_geoclaw_sample:

*****************************************************************
Sample GeoClaw `setrun.py` module
*****************************************************************

This shows an example of a :mod:`setrun.py` module for a GeoClaw
application.  Note that *claw_pkg='geoclaw'* is expected.  The line::

        rundata = data.ClawRunData(claw_pkg, ndim)

creates an object of class :class:`ClawRunData` with a :attr:`geodata`
attribute.  

The new parameters needed by GeoClaw are set in the function
`setrun_setgeo`_ below.


::


    """ 
    Module to set up run time parameters for Clawpack.
    
    The values set in the function setrun are then written out to data files
    that will be read in by the Fortran code.
        
    """ 
    
    import os
    from pyclaw import data 
    
    
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
        
        assert claw_pkg.lower() == 'geoclaw',  "Expected claw_pkg = 'geoclaw'"
    
        ndim = 2
        rundata = data.ClawRunData(claw_pkg, ndim)
    
        #------------------------------------------------------------------
        # GeoClaw specific parameters:
        #------------------------------------------------------------------
    
        rundata = setgeo(rundata)   # Defined below
        
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
        clawdata.ndim = ndim
        
        # Lower and upper edge of computational domain:
        clawdata.xlower = -100.
        clawdata.xupper = 100.
        
        clawdata.ylower = -100.
        clawdata.yupper = 100.
            
    
        # Number of grid cells:
        clawdata.mx = 50
        clawdata.my = 50
            
    
        # ---------------
        # Size of system:
        # ---------------
    
        # Number of equations in the system:
        clawdata.meqn = 3
    
        # Number of auxiliary variables in the aux array (initialized in setaux)
        clawdata.maux = 3
        
        # Index of aux array corresponding to capacity function, if there is one:
        clawdata.mcapa = 0
        
        
        
        # -------------
        # Initial time:
        # -------------
    
        clawdata.t0 = 0.0
        
        
        # -------------
        # Output times:
        #--------------
    
        # Specify at what times the results should be written to fort.q files.
        # Note that the time integration stops after the final output time.
        # The solution at initial time t0 is always written in addition.
    
        clawdata.outstyle = 1
    
        if clawdata.outstyle==1:
            # Output nout frames at equally spaced times up to tfinal:
            clawdata.nout = 14
            clawdata.tfinal = 7.0
    
        elif clawdata.outstyle == 2:
            # Specify a list of output times.  
            clawdata.tout =  [0.5, 1.0]   # used if outstyle == 2
            clawdata.nout = len(clawdata.tout)
    
        elif clawdata.outstyle == 3:
            # Output every iout timesteps with a total of ntot time steps:
            iout = 1
            ntot = 5
            clawdata.iout = [iout, ntot]
        
    
    
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
        clawdata.dt_variable = 1
        
        # Initial time step for variable dt.  
        # If dt_variable==0 then dt=dt_initial for all steps:
        clawdata.dt_initial = 0.016
        
        # Max time step to be allowed if variable dt used:
        clawdata.dt_max = 1e+99
        
        # Desired Courant number if variable dt used, and max to allow without 
        # retaking step with a smaller dt:
        clawdata.cfl_desired = 0.5
        clawdata.cfl_max = 1.0
        
        # Maximum number of time steps to allow between output times:
        clawdata.max_steps = 5000
    
        
        
    
        # ------------------
        # Method to be used:
        # ------------------
    
        # Order of accuracy:  1 => Godunov,  2 => Lax-Wendroff plus limiters
        clawdata.order = 2
        
        # Transverse order for 2d or 3d (not used in 1d):
        clawdata.order_trans = 2
        
        # Number of waves in the Riemann solution:
        clawdata.mwaves = 3
        
        # List of limiters to use for each wave family:  
        # Required:  len(mthlim) == mwaves
        clawdata.mthlim = [3,3,3]
        
        # Source terms splitting:
        #   src_split == 0  => no source term (src routine never called)
        #   src_split == 1  => Godunov (1st order) splitting used, 
        #   src_split == 2  => Strang (2nd order) splitting used,  not recommended.
        clawdata.src_split = 1
        
        
        # --------------------
        # Boundary conditions:
        # --------------------
    
        # Number of ghost cells (usually 2)
        clawdata.mbc = 2
        
        # Choice of BCs at xlower and xupper:
        #   0 => user specified (must modify bcN.f to use this option)
        #   1 => extrapolation (non-reflecting outflow)
        #   2 => periodic (must specify this at both boundaries)
        #   3 => solid wall for systems where q(2) is normal velocity
        
        clawdata.mthbc_xlower = 1
        clawdata.mthbc_xupper = 1
        
        clawdata.mthbc_ylower = 1
        clawdata.mthbc_yupper = 1
        
    
        # ---------------
        # AMR parameters:
        # ---------------
    
    
        # max number of refinement levels:
        mxnest = 3
    
        clawdata.mxnest = -mxnest   # negative ==> anisotropic refinement in x,y,t
    
        # List of refinement ratios at each level (length at least mxnest-1)
        clawdata.inratx = [2,4,2]
        clawdata.inraty = [2,4,2]

        clawdata.inratt = [2,4,2]
        # Instead of setting these ratios, set the next parameter to True
        # to automatically choose refinement ratios in time based on estimate
        # of maximum wave speed on all grids at each level.
        variable_dt_refinement_ratios = False

    
    
        # Specify type of each aux variable in clawdata.auxtype.
        # This must be a list of length maux, each element of which is one of:
        #   'center',  'capacity', 'xleft', or 'yleft'  (see documentation).
    
        clawdata.auxtype = ['center','center','yleft']
    
    
        clawdata.tol = -1.0     # negative ==> don't use Richardson estimator
        clawdata.tolsp = 0.5    # used in default flag2refine subroutine
                                # (Not used in geoclaw!)
    
        clawdata.kcheck = 2     # how often to regrid (every kcheck steps)
        clawdata.ibuff  = 2     # width of buffer zone around flagged points
    
        # More AMR parameters can be set -- see the defaults in pyclaw/data.py
    
        return rundata
        # end of function setrun
        # ----------------------
    

.. comment
    _setrun_geoclaw_sample_parameters:

    GeoClaw parameters
    ------------------

.. _setrun_setgeo:

setgeo function
---------------
    
::

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
    
        # == setgeo.data values ==
        geodata.igravity = 1
        geodata.gravity = 9.81
        geodata.icoordsys = 1
    
        # == settsunami.data values ==
        geodata.sealevel = 0.
        geodata.drytolerance = 1.e-3
        geodata.wavetolerance = 1.e-2
        geodata.depthdeep = 1.e2
        geodata.maxleveldeep = 3
        geodata.ifriction = 0
        geodata.coeffmanning = 0.
        geodata.frictiondepth = 20.
    
        # == settopo.data values ==
        geodata.topofiles = []
        # for topography, append lines of the form
        #   [topotype, minlevel, maxlevel, t1, t2, fname]
        geodata.topofiles.append([2, 1, 1, 0., 1.e10, 'bowl.topotype2'])
    
        # == setdtopo.data values ==
        geodata.dtopofiles = []
        # for moving topography, append lines of the form:  (<= 1 allowed for now!)
        #   [minlevel,maxlevel,fname]
    
        # == setqinit.data values ==
        geodata.iqinit = 1
        geodata.qinitfiles = []  
        # for qinit perturbations, append lines of the form: (<= 1 allowed for now!)
        #   [minlev, maxlev, fname]
        geodata.qinitfiles.append([1, 2, 'hump.xyz'])
    
        # == setregions.data values ==
        geodata.regions = []
        # to specify regions of refinement append lines of the form
        #  [minlevel,maxlevel,t1,t2,x1,x2,y1,y2]
        geodata.regions.append([1, 1, 0.e0, 1.e10, -100.,100., -100.,100.])
        geodata.regions.append([1, 2, 0.e0, 1.e10,    0.,100.,    0.,100.])
        geodata.regions.append([2, 3, 3.e0, 1.e10,   52., 72.,   52., 72.])
        geodata.regions.append([2, 3, 3.e0, 1.e10,   75., 95.,   -10.,  10.])
    
        # == setgauges.data values ==
        geodata.gauges = []
        # for gauges append lines of the form  [gaugeno, x, y, t1, t2]
        from numpy import linspace
        # gauges along x-axis:
        gaugeno = 0
        for r in linspace(86., 93., 9):
            gaugeno = gaugeno+1
            x = r + .0005  # shift a bit away from cell corners
            y = .0005
            geodata.gauges.append([gaugeno, x, y, 0., 1e10])
    
    
        # == setfixedgrids.data values ==
        geodata.fixedgrids = []
        # for fixed grids append lines of the form
        # [t1,t2,noutput,x1,x2,y1,y2,xpoints,ypoints,\
        #  ioutarrivaltimes,ioutsurfacemax]
        geodata.fixedgrids.append([1., 2., 4, 0., 100., 0., 100., 11, 11, 0, 0])
    
        return rundata
        # end of function setgeo
        # ----------------------
    
    
    
    if __name__ == '__main__':
        # Set up run-time parameters and write all data files.
        import sys
        if len(sys.argv) == 2:
    	rundata = setrun(sys.argv[1])
        else:
    	rundata = setrun()
    
        rundata.write()
        
