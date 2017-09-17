
.. _geoclaw_examples_storm_surge_hurricane_ike:

Basic Example of Storm Surge from Hurricane Ike
===============================================

This example contains the data and setup for running a storm surge forecast for
Hurricane Ike.  The example can be run via::

    make all

which should download the necessary topography and storm data, start the 
simulation and plot the results.  

For the example to run in a reasonable amount of time the 
resolution has been limited to only the first two levels or refinement.  If you
would like to run at the full resolution matching that which was run in
[1] change line 280 in *setrun.py* to be::

    amrdata.amr_levels_max = 5

There is also additional topography data needed for the region around Galveston
Bay which can be obtained from NOAA or by contacting Kyle Mandli.  Note also 
that the observed tide gauge data plotted in [1] was obtained from 
Andrew Kennedy of the University of Notre Dame.

1.	Mandli, K. T. & Dawson, C. N. Adaptive Mesh Refinement for Storm Surge. Ocean Modelling 75, 36–50 (2014).
