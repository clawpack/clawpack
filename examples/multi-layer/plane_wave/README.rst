
.. _geoclaw_examples_multilayer_plane_wave:

Basic 2D Multi-Layer Shallow Water Examples
===========================================

These examples include simple setups for plane wave tests of the 2-layer 
shallow water equations with a plane wave initial condition and a jump in 
bathymetry.  The location and angles of the plane wave and bathymetry jump are 
controlled via the *setrun.py* file.  Additionally the ability to run with a 
similar topography but a gaussian hump as the initial condition has also been
included.  To run the basic plane-wave example use::

    make .plots

Running the *run_tests.py* script will execute a family of different setups for
both plane-waves incident at different angles to the grid and topography and a 
set of the ``bubble tests``.  The basics of the Riemann solvers involved in
these examples can be found in [1].

1.	Mandli, K. T. A Numerical Method for the Two Layer Shallow Water Equations with Dry States. Ocean Modelling 72, 80–91 (2013).
