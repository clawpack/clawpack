
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

`setplot_speeds.py` is a version of setplot that also plots velocities.


**Creating kml files to view on Google Earth**

By default `make data` will create kml files showing the computational
domain, the extents of the topofiles and dtopofiles used and the
AMR regions specified, and the gauge location.  These can be viewed
using Google Earth or other platforms supporting `kml`.  These are
created by the line in `setrun.py` invoking
`kmltools.make_input_data_kmls(rundata)`.

Note that clicking on a rectangle or gauge in Google Earth displays
information about it.  You may need to unselect some of the kml files under
Places in order to click on ones underneath.

`setplot_kml.py` is a version of the setplot file that produces `kml` files
from the GeoClaw output for viewing the tsunami propagation.

To use::

    make .output
    make plots SETPLOT_FILE=setplot_kml.py

and then open `_plots/Chile_2010.kmz` from Google Earth.  Note that this is
a zipped directory containing a number of `kml` files and opening this way
should open them all.  You can then select which ones to view (e.g. the
simultion, the grid patches, regions specified in `setrun.py`, and/or gauges).

You can `unzip` the `kmz` file if you want to extract individual `kml`
files.

See `<http://www.clawpack.org/googleearth_plotting.html>`_ for
documentation.

