
.. _geoclaw_examples_tsunami_chile2010:

Tsunami arising offshore Maule, Chile, Feb. 27, 2010 
=====================================================

This uses a very simple source model with a single fault plane and constant
slip, from an early inversion by the USGS.  The input parameters are
specified in `usgs100227.cfg` and converted into the dtopo file
`usgs100227.txt` by the script `maketopo.py`.  This can be run via::

    make topo

which also downloads a topo file for the ocean bathymetry.
This bathymetry originally came from the NOAA National Geophysical Data
Center (NGDC)
using `Design-a-grid <http://www.ngdc.noaa.gov/mgg/gdas/gd_designagrid.html>`_ .

A single gauge captures the sea surface elevation at the location of 
`DART buoy 32412
<http://www.ndbc.noaa.gov/station_page.php?station=32412>`_.

setplot_speeds.py is a version of setplot that also plots velocities.
