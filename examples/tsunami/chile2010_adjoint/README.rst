
.. _geoclaw_examples_tsunami_chile2010_adjoint:

Tsunami arising offshore Maule, Chile, Feb. 27, 2010 
=====================================================

This is a modified version of `$CLAW/geoclaw/examples/tsunami/chile2010`.
See that directory and README for more information about the problem and data.

Adjoint flagging
----------------

The adjoint method is used to flag cells needing refinement, as described in
the papers:

- Adjoint Methods for Guiding Adaptive Mesh Refinement in Tsunami Modeling 
  by Brisa N. Davis and R. J. LeVeque, Pure Appl. Geophys. 173 (2016), pp.
  4055-4074. 
  `<http://faculty.washington.edu/rjl/pubs/adjoint2016>`_

- Analysis and Performance Evaluation of Adjoint-Guided Adaptive Mesh
  Refinement for Linear Hyperbolic PDEs Using Clawpack, by
  B. N. Davis and R. J. LeVeque, 2018.
  `<http://faculty.washington.edu/rjl/pubs/adjoint2018>`_




Folder Organization
--------------------

- **adjoint:**

  Contains code to solve the adjoint problem.

  The output times specified in this directory should to at least as
  far out in time as the forward solution is desired, with sufficiently
  dense outputs to properly capture the evolving adjoint solution.

Running the Code
--------------------

Go to the folder `adjoint` and run in a terminal::

    make topo  # downloads topo data and creates adjoint initial hump
    make new   # compile everything
    make .plots

The code will produce two new folders: _output and _plots. 
The first one contains all the output files, while the latter one contains
the plots and interactive visualization apps.

Then return to this directory and 

    make topo  # created dtopo file modeling earthquake
    make new
    make .plots

to run the forward solution and make plots that also show the inner product
of the forward and adjoint solution on levels 1 and 2 (not on level 3 since 
there is no need to flag further in this 3-level run).

Alternatively, to run first the adjoint and then the forward problem a
single script can be invoked.  

Running Variations
--------------------

In `setrun.py`, the following flags are currently set (in various places)::

    adjointdata.use_adjoint = True

    # Flag for refinement using routine flag2refine:
    amrdata.flag2refine = True
    rundata.amrdata.flag2refine_tol = 0.0005

    # time period of interest:
    adjointdata.t1 = rundata.clawdata.t0
    adjointdata.t2 = rundata.clawdata.tfinal

Setting `adjointdata.use_adjoint` to `False` will go back to using standard
flagging based on the magnitude of undivided differences or an estimate of
the one-step error.  Flagging based on Richardson estimate of the 1-step
error does not generally work well in GeoClaw (with or without the adjoint
inner product), so using `flag2refine` is recommended.

The refinement tolerance is set appropriate for the adjoint flagging. Note
that waves are no longer refined after passing the DART gauge.
The location of interest is specified in `adjoint/maketopo.py`, where
the functional used as initial data (at the final time) in the adjoint
problem is set to be a Gaussian hump centered at the DART gauge location.

The time period of interest can be changed to some subset of the full run
time.  Try changing this to see how the AMR adapts to only capture waves
reaching the gauge over a shorter specified time period, e.g. try `t1 =
3*3600.` and `t2 = 4.5*3600.` to capture only the leading wave at the gauge.


