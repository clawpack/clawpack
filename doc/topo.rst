

.. _topo:

*****************************************************************
Topography data
*****************************************************************

The :ref:`geoclaw` software for flow over topography requires at least one
topo file to be input, see :ref:`setrun_geoclaw`.

Currently topo files are restricted to three possible formats as ASCII files.
A future project is to allow other formats including NetCDF.

In the descriptions below it is assumed that the topo file gives the
elevation of the topography (relative to some reference level) as a value of
z at each (x,y) point on a rectangular grid.  Only uniformly spaced
rectangular topo grids are currently recognized.  

More than one topo file can be specified (see :ref:`setrun_topo`) that might
cover overlapping regions at different resolutions.  The union of all the
topo files should cover the full computational domain specified (and may
extend outside it).  Internally in :ref:`geoclaw` a single
piecewise-bilinear function is constructed from the union of the topo files,
using the best information available in regions of overlap.  This function
is then integrated over computational grid cells to obtain the single topo value
in each grid cell needed when solving depth averaged equations such as the
shallow water equations with these finite volume methods.  Note that this
has the feature that if a grid cell is refined at some stage in the
computation, the topo used in the fine cells have an average value that is
equal to the coarse cell value.  This is crucial in maintaining the
ocean-at-rest steady state, for example.

The recognized topotypes are:

  **topotype = 1**

    x,y,z values on each line, progressing from upper left (NW) corner across
    rows (moving east), then down in standard GIS form.  
    The size of the grid and spacing
    between the grid points is deduced from the data.  

    *Example:* if you want a flat bottom at B = -1000.
    over a domain  0. <= x <= 10. and  20. <= y <= 30.
    then the topo file could be simply::

        0.  30.  -1000.
        10. 30.  -1000.
        0.  20.  -1000.
        10. 20.  -1000.

    These files are larger than necessary since they store the x,y values at
    each point even though the points are required to be equally spaced.
    Many data sets come this way, but note that you can convert a file of
    this type to one of the more compact types below using
    `pyclaw.geotools.topotools.converttopotype(inputfile, outputfile,
    topotypein=1, topotypeout=2, nodata_value=None)`.



  **topotype = 2**

    The file starts with a header consisting of 6 lines containing::

      mx
      my
      xllcorner
      yllcorner
      cellsize
      nodataval

    and is followed by mx*my lines containing the z values at each x,y,
    again progressing from upper left (NW) corner across
    rows (moving east), then down in standard GIS form.  
    The lower left corner of the grid
    is *(xllcorner, yllcorner)* and the distance between grid points in both
    x and y is *cellsize*.  The value *nodataval* indicates what value of z
    is specified for missing data points (often something like 9999. in data
    sets with missing values).

    *Example:*  For the same example as above, the topo file with
    topotype==2 would be::

      2         mx
      2         my
      0.        xllcorner
      20.       yllcorner
      10.       cellsize
      9999.     nodatavalue
      -1000.
      -1000.
      -1000.
      -1000.


  **topotype = 3**

    The file starts with a header consisting of 6 lines as for *topotype=2*,
    followed by *my* lines, each containing *mx* values for one row of data
    (ordered as before, so the first line of data is the northernmost line
    of data, going from west to east).

    *Example:*  For the same example as above, the topo file with
    topotype==3 would be::

      2         mx
      2         my
      0.        xllcorner
      20.       yllcorner
      10.       cellsize
      9999.     nodatavalue
      -1000.  -1000.
      -1000.  -1000.


It is also possible to specify values -1, -2, or -3 for *topotype*, in which
case the *z* values will be negated as they are read in (since some data
sets use different convensions for positive and negative values relative to
sea level). 

For :ref:`geoclaw` applications in the ocean or lakes (such as tsunami
modeling), it is generally assumed that *sealevel = 0* has been set in
:ref:`setrun_tsunami` and that *z<0* corresponds to subsurface bathymetry
and *z>0* to topograpy above sea level.

.. _topo_sources:

Downloading topography files
----------------------------

The example
`$CLAW/apps/tsunami/chile2010
<claw/apps/tsunami/chile2010/README.html>`_
is set up to automatically download topo files via::

	$ make topo

See the `maketopo.py <claw/apps/tsunami/chile2010/maketopo.py.html>`_
file in that directory.

Other such examples will appear in the future.  

Several on-line databases are available for topograpy, e.g.

 * NOAA National Geophysical Data Center (NGDC)
   `Design-a-grid <http://www.ngdc.noaa.gov/mgg/gdas/gd_designagrid.html>`_


.. _topo_dtopo:

Topography displacement files
-----------------------------

For tsunami generation a file *dtopo* is generally used to specify the
displacement of the topography relative to that specified in the topo files.

Currently only one format is allowed for this file: it must be similar to
topo files with *topotype=1* as described above, except that each line
starts with a *t* value for the time, so each line contains t,x,y,dz

The x,y,dz values give the displacement dz at x,y at time t.  It is assumed
that the grid is uniform and that the file contains mx*my*mt lines if mt
different times are specified for an mx*my grid.  

**To do:** it would be better to have a header as for *topotype=2,3* that also
lists *dt* and then just list the *dz* values.

The Okada model can be used to generate *dtopo* files from fault parameters.
See `$CLAW/apps/tsunami/chile2010/maketopo.py
<claw/apps/tsunami/chile2010/maketopo.py.html>`_ for an example.


.. _qinit_file:

qinit data file
---------------

Instead of (or in addition to) specifying a displacement of the topography
it is possible to specify a perturbation to the depth, momentum, or surface
elevation of the initial data.  This is generally useful only for tsunami
modeling where the initial data specified in the default *qinit_geo.f* function
is the stationary water with surface elevation equal to *sealevel* as set in
:ref:`setrun_tsunami`.  

Of course it is possible to copy the *qinit_geo.f* function to your
directory and modify it, but for some applications the initial elevation may
be given on grid of the same type as described above.  In this case file can
be provided as described at :ref:`setrun_qinit` containing this
perturbation.

The file format is similar to what is described above for *topotype=1*, but
now each line contains *x,y,dq* where *dq* is a perturbation to one of the 
components of *q* as specified by the value of *iqinit* specified (see
:ref:`setrun_qinit`).  If *iqinit = 4*, the value *dq* is instead the
surface elevation desired for the initial data and the depth *h* (first
component of *q*) is set accordingly.

