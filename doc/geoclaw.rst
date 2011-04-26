

.. _geoclaw:

***************
GeoClaw
***************

The `$CLAW/geoclaw` directory contains a specialized version of some Clawpack
and AMRClaw routines that have been modified to work well for certain
geophysical flow problems.  

.. warning:: As with all of Clawpack, this code is provided as a research
   and teaching tool with no guarantee of suitability for any particular
   purpose, and no liability on the part of the authors.  See the
   :ref:`license` for more details.

Currently the focus is on 2d depth-averaged
shallow water equations for flow over varying topography.  The term
*bathymetry* is often used for underwater topography (sea floor or lake
bottom), but in this documentation and in the code the term *topography* is
often used to refer to either.

A primary concern with such flows is handling the margins of the flow where
the depth goes to 0, for example at the shore line.  In GeoClaw this is
handled by letting the depth variable *h* in the shallow water equations be
0 in some cells.  Robust Riemann solvers are used that allow for dry cells
adjacent to wet cells and that allow wetting and drying, for example as a
tsunami inundates dry land.

Some sample calculations can be viewed in the :ref:`apps`.



.. _geoclaw_run:

Running a GeoClaw code
----------------------

Setting up, running, and plotting a GeoClaw application follows the same pattern
as other AMRClaw applications, which in turn use many of the same
conventions as the classic single grid Clawpack code, in particular:

 * Setting parameters is done in `setrun.py`, as for other versions
   of Clawpack, as described in :ref:`setrun`.  However, there are several
   new parameters that may or must be set for GeoClaw.  See
   :ref:`setrun_geoclaw` for more details on these.

 * The program can be compiled and run using *make* and *make .output* as
   for other versions, see :ref:`fortran`.

 * Plots of results can be created either as a set of webpages via
   *make .plots* or interactively using *Iplotclaw*.  See
   :ref:`plotting` for more details.  Some additional Python plotting tools 
   that are useful for GeoClaw output (e.g. plotting land and water with
   different colormaps) are described in the section
   :ref:`plotting_geoclaw`.


.. _topo_intro:

Topography
----------

To simulate  flow over topography it is of course necessary to specify 
the topography.  This is usually done by providing one or more files of
surface elevation (relative to some reference, e.g. sea level) at a set of
points on a rectangular grid (with x-y locations in Cartesian units or in
latitude-longitude, depending on the application).

Several file formats are recognized by GeoClaw.  See :ref:`topo` for more
information on how to specify topography and some on-line resources for
obtaining topography.

.. _geoclaw_plotting:

Plotting GeoClaw results
------------------------

GeoClaw results can be plotted with the usual Python plotting tools (see
:ref:`plotting`).  

Some special tools and colormaps are available, see :ref:`geoplot`.

Setting up a new example
------------------------

 * Hints to appear.

