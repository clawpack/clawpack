r"""
kmltools module: $CLAW/geoclaw/src/python/geoclaw/kmltools.py

Tools to make kml files to overlay on Google Earth.
Note that color is in KML format, BGR with 2 hex digits for each, e.g.

  FF0000 is blue, 00FF00 is green,  0000FF is red, 00FF00 is yellow.

Actually it's an 8 hex digit number, where the first two digits are
transparency, but in this module these default to 'FF' (but you can specify
the full 8 digits if you want it transparent).

:Functions:
 - deg2dms - convert decimal degrees to (degrees, minutes, seconds)
 - regions2kml - create a kml outline for each regions specified in setrun
 - box2kml - create a kml outline from a rectangular box
 - quad2kml - create a kml outline for an arbitrary quadrilateral
 - poly2kml - create a kml outline for an arbitrary polygon
 - line2kml - create a kml line connecting 2 points
 - gauges2kml - create a kml marker for each gauge specified in setrun
 - kml_header - used internally
 - kml_footer - used internally
 - kml_region - used internally
 - kml_gauge - used internally

 - strip_archive_extensions - strip off things like .tar or .gz
"""


from __future__ import absolute_import
from __future__ import print_function
from six.moves import range

def f2s(x, num_digits=6):
    r"""
    Convert float to string in fixed point notation with at most
    *num_digits* digits of precision and trailing zeros removed, 
    for printing nicely in kml description boxes.
    """
    format = '%' + '.%sf' % num_digits
    s = (format % x).rstrip('0')
    return s
    
def deg2dms(dy):
    r"""
    Convert decimal degrees to tuple (degrees, minutes, seconds)
    """

    from numpy import floor
    dy_deg = floor(dy)
    dy_min = floor((dy-dy_deg)*60.)
    dy_sec = (dy-dy_deg-dy_min/60.)*3600.
    return dy_deg,dy_min,dy_sec


def regions2kml(rundata=None,fname='regions.kml',verbose=True,combined=True):

    """
    Create a KML box for each AMR region specified for a GeoClaw run.

    :Inputs:

      - *rundata* - an object of class *ClawRunData* or None

        If *rundata==None*, try to create based on executing function *setrun*
        from the `setrun.py` file in the current directory.

      - *fname* (str) - resulting kml file.

      - *verbose* (bool) - If *True*, print out info about each region found

      - *combined* (bool) - If *True*, combine into single kml file with
        name given by *fname*.  This is the default. 
        If False, *fname* is ignored and individual files are created for
        each region with names are Domain.kml, Region00.kml, etc.
        These will show up separately in GoogleEarth so they can be turned
        on or off individually.

    First create a box for the entire domain (in red) and then a box
    for each region (in white).

    :Example:

        >>> from clawpack.geoclaw import kmltools
        >>> kmltools.regions2kml()

    is equivalent to:

        >>> from clawpack.geoclaw import kmltools
        >>> from setrun import setrun
        >>> rundata = setrun()
        >>> kmltools.regions2kml(rundata)

    By default this creates a file named *regions.kml* that can be opened in
    Google Earth.

    """

    from numpy import cos,pi,floor

    if rundata is None:
        try:
            import setrun
            reload(setrun)
            rundata = setrun.setrun()
        except:
            raise IOError("*** cannot execute setrun file")

    clawdata = rundata.clawdata
    x1,y1 = clawdata.lower[0:]
    x2,y2 = clawdata.upper[0:]
    description = "  x1 = %s, x2 = %s\n" % (f2s(x1),f2s(x2)) \
            + "  y1 = %s, y2 = %s\n" % (f2s(y1),f2s(y2))

    mx,my = clawdata.num_cells[0:]
    dx = (x2-x1)/float(mx)
    dx_meters = dx*111e3*cos(pi*0.5*(y1+y2)/180.)
    dy = (y2-y1)/float(my)
    dy_meters = dy*111e3
    if verbose:
        print("Domain:   %10.6f  %10.6f  %10.6f  %10.6f" % (x1,x2,y1,y2))
    dx_deg,dx_min,dx_sec = deg2dms(dx)
    dy_deg,dy_min,dy_sec = deg2dms(dy)
    #print "Level 1 resolution:  dx = %g deg, %g min, %g sec = %g meters" \
    #   % (dx_deg,dx_min,dx_sec,dx_meters)
    levtext = "Level 1 resolution:  dy = %g deg, %g min, %g sec = %g meters\n" \
        % (dy_deg,dy_min,dy_sec,dy_meters)
    if verbose:
        print(levtext)
    description = description + levtext

    amr_levels_max = rundata.amrdata.amr_levels_max
    refinement_ratios_y = rundata.amrdata.refinement_ratios_y
    num_ref_ratios = len(refinement_ratios_y)
    if amr_levels_max > num_ref_ratios+1:
        raise IOError("*** Too few refinement ratios specified for " \
            + "amr_levels_max = %i" % amr_levels_max)
    dy_levels = (num_ref_ratios+1) * [dy]
    for k,r in enumerate(refinement_ratios_y):
        level = k+2
        dy = dy_levels[k] / r
        dy_levels[k+1] = dy
        dy_meters = dy*111e3
        dy_deg,dy_min,dy_sec = deg2dms(dy)
        levtext = "Level %s resolution:  dy = %g deg, %g min, %g sec = %g meters  (refined by %i)\n" \
                % (level,dy_deg,dy_min,dy_sec,dy_meters,r)
        if verbose:
            print(levtext)
        description = description + levtext

    if verbose:
        print("Allowing maximum of %i levels" % amr_levels_max)

    elev = 0.
    if not combined:
        fname = 'Domain.kml'

    kml_text = kml_header(fname)

    mapping = {}
    mapping['x1'] = x1
    mapping['x2'] = x2
    mapping['y1'] = y1
    mapping['y2'] = y2
    mapping['elev'] = elev
    mapping['name'] = 'Computational Domain'
    mapping['desc'] = description
    mapping['color'] = "0000FF"  # red
    mapping['width'] = 2

    region_text = kml_region(mapping)
    kml_text = kml_text + region_text

    if not combined:
        kml_text = kml_text + kml_footer()
        kml_file = open(fname,'w')
        kml_file.write(kml_text)
        kml_file.close()
        if verbose:
            print("Created ",fname)

            

    regions = rundata.regiondata.regions
    if len(regions)==0 and verbose:
        print("No regions found in setrun.py")


    for rnum,region in enumerate(regions):
        if not combined:
            fname = 'Region_%s.kml' % str(rnum).zfill(2)
            kml_text = kml_header(fname)

        minlevel,maxlevel = region[0:2]
        t1,t2 = region[2:4]
        x1,x2,y1,y2 = region[4:]

        if verbose:
            print("Region %i: %10.6f  %10.6f  %10.6f  %10.6f" \
                    % (rnum,x1,x2,y1,y2))
            print("           minlevel = %i,  maxlevel = %i" \
                    % (minlevel,maxlevel) \
                    + "  t1 = %s,  t2 = %s" % (f2s(t1),f2s(t2)))
        mapping = {}
        mapping['minlevel'] = minlevel
        mapping['maxlevel'] = maxlevel
        mapping['t1'] = t1
        mapping['t2'] = t2
        mapping['x1'] = x1
        mapping['x2'] = x2
        mapping['y1'] = y1
        mapping['y2'] = y2
        mapping['elev'] = elev
        mapping['name'] = 'Region %i' % rnum
        description = "minlevel = %i, maxlevel = %i\n" % (minlevel,maxlevel) \
            + "  t1 = %s, t2 = %s\n" % (f2s(t1),f2s(t2)) \
            + "  x1 = %s, x2 = %s\n" % (f2s(x1),f2s(x2)) \
            + "  y1 = %s, y2 = %s\n\n" % (f2s(y1),f2s(y2))
        if len(dy_levels) >= minlevel:
            dy = dy_levels[minlevel-1]
            dy_deg,dy_min,dy_sec = deg2dms(dy)
            dy_meters = dy*111e3
            levtext = "Level %s resolution:  \ndy = %g deg, %g min, %g sec \n= %g meters\n" \
                    % (minlevel,dy_deg,dy_min,dy_sec,dy_meters)
            description = description + levtext
        if (maxlevel > minlevel) and (len(dy_levels) >= maxlevel):
            dy = dy_levels[maxlevel-1]
            dy_deg,dy_min,dy_sec = deg2dms(dy)
            dy_meters = dy*111e3
            levtext = "\nLevel %s resolution:  \ndy = %g deg, %g min, %g sec \n= %g meters\n" \
                    % (maxlevel,dy_deg,dy_min,dy_sec,dy_meters)
            description = description + levtext
        mapping['desc'] = description
        mapping['color'] = "FFFFFF"  # white
        mapping['width'] = 3

        region_text = kml_region(mapping)
        kml_text = kml_text + region_text
        if not combined:
            kml_text = kml_text + kml_footer()
            kml_file = open(fname,'w')
            kml_file.write(kml_text)
            kml_file.close()
            if verbose:
                print("Created ",fname)

    if combined:
        kml_text = kml_text + kml_footer()
        kml_file = open(fname,'w')
        kml_file.write(kml_text)
        kml_file.close()
        if verbose:
            print("Created ",fname)


def line2kml(xy,fname='line.kml',name='line',color='00FFFF',width=3,
             verbose=True):
    """
    Make a KML line with default color yellow.

    :Inputs:

     - *xy* a tuple ((x1,x2),(y1,y2)) (preferred) 
            or (x1,x2,y1,y2) (for backward compatibility)
     - *fname* (str) name of resulting kml file
     - *name* (str) name to appear on line on Google Earth
     - *color* (str) Color in format aabbggrr
     - *width* (str) line width
     - *verbose* (bool) - If *True*, print out info

    """
     
    if type(xy[0]) is tuple:
        x1,x2 = xy[0]
        y1,y2 = xy[1]
    else:
        x1,x2,y1,y2 = xy[0:]

    if verbose:
        print("Line:   %10.6f  %10.6f  %10.6f  %10.6f" % (x1,x2,y1,y2))

    elev = 0.
    kml_text = kml_header(fname)

    mapping = {}
    mapping['x1'] = x1
    mapping['x2'] = x2
    mapping['y1'] = y1
    mapping['y2'] = y2
    mapping['elev'] = elev
    mapping['name'] = name
    mapping['desc'] = "  x1 = %s, x2 = %s\n" % (f2s(x1),f2s(x2)) \
            + "  y1 = %s, y2 = %s" % (f2s(y1),f2s(y2))
    mapping['color'] = color
    mapping['width'] = width

    region_text = kml_line(mapping)

    kml_text = kml_text + region_text + kml_footer()
    kml_file = open(fname,'w')
    kml_file.write(kml_text)
    kml_file.close()
    if verbose:
        print("Created ",fname)


def box2kml(xy,fname=None,name='box',color='FF0000',width=3,verbose=True):
    """
    Make a KML box with default color blue.

    :Inputs:

     - *xy* a tuple ((x1,x2),(y1,y2)) (preferred) 
            or (x1,x2,y1,y2) (for backward compatibility)
     - *fname* (str) name of resulting kml file
     - *name* (str) name to appear in box on Google Earth
     - *color* (str) Color in format aabbggrr
     - *width* (str) line width
     - *verbose* (bool) - If *True*, print out info

    """

    if fname is None:
        fname = name + '.kml'

    if type(xy[0]) is tuple:
        x1,x2 = xy[0]
        y1,y2 = xy[1]
    else:
        x1,x2,y1,y2 = xy[0:]

    if verbose:
        print("Box:   %10.6f  %10.6f  %10.6f  %10.6f" % (x1,x2,y1,y2))

    elev = 0.
    kml_text = kml_header(fname)

    mapping = {}
    mapping['x1'] = x1
    mapping['x2'] = x2
    mapping['y1'] = y1
    mapping['y2'] = y2
    mapping['elev'] = elev
    mapping['name'] = name
    mapping['desc'] = "  x1 = %s, x2 = %s\n" % (f2s(x1),f2s(x2)) \
            + "  y1 = %s, y2 = %s" % (f2s(y1),f2s(y2))
    mapping['color'] = color
    mapping['width'] = width

    region_text = kml_region(mapping)

    kml_text = kml_text + region_text + kml_footer()
    kml_file = open(fname,'w')
    kml_file.write(kml_text)
    kml_file.close()
    if verbose:
        print("Created ",fname)


def quad2kml(xy,fname=None,name='quad',color='FF0000',width=3,verbose=True):
    """
    Make a KML quadrilateral with default color blue.

    :Inputs:

     - *xy* a tuple ((x1,x2,x3,x4),(y1,y2,y3,y4)) (preferred) 
            or (x1,x2,y1,y2,x3,y3,x4,y4) (for backward compatibility)
     - *fname* (str) name of resulting kml file
     - *name* (str) name to appear in box on Google Earth
     - *color* (str) Color in format aabbggrr
     - *width* (str) line width
     - *verbose* (bool) - If *True*, print out info

    """

    if fname is None:
        fname = name + '.kml'

    if type(xy[0]) is tuple:
        x1,x2,x3,x4 = xy[0]
        y1,y2,y3,y4 = xy[1]
    else:
        x1,y1,x2,y2,x3,y3,x4,y4 = xy[0:]

    if verbose:
        print("Quadrilateral:   %10.6f  %10.6f" % (x1,y1))
        print("                 %10.6f  %10.6f" % (x2,y2))
        print("                 %10.6f  %10.6f" % (x3,y3))
        print("                 %10.6f  %10.6f" % (x4,y4))

    elev = 0.
    kml_text = kml_header(fname)

    mapping = {}
    mapping['x1'] = x1
    mapping['x2'] = x2
    mapping['x3'] = x3
    mapping['x4'] = x4
    mapping['y1'] = y1
    mapping['y2'] = y2
    mapping['y3'] = y3
    mapping['y4'] = y4
    mapping['elev'] = elev
    mapping['name'] = name
    mapping['desc'] = "  x1 = %s, y1 = %s\n" % (f2s(x1),f2s(y1)) \
            + "  x2 = %s, y2 = %s" % (f2s(x2),f2s(y2)) \
            + "  x3 = %s, y3 = %s" % (f2s(x3),f2s(y3)) \
            + "  x4 = %s, y4 = %s" % (f2s(x4),f2s(y4))
    mapping['color'] = color
    mapping['width'] = 3

    region_text = kml_region(mapping)

    kml_text = kml_text + region_text + kml_footer()
    kml_file = open(fname,'w')
    kml_file.write(kml_text)
    kml_file.close()
    if verbose:
        print("Created ",fname)


def poly2kml(xy,fname=None,name='poly',color='00FF00', width=3,
             verbose=True):
    """
    Make a KML polygon with default color blue.

    :Inputs:

     - *xy* a tuple (x,y) where x and y are lists of vertices
     - *fname* (str) name of resulting kml file
     - *name* (str) name to appear in box on Google Earth
     - *color* (str) Color in format aabbggrr
     - *width* (str) line width
     - *verbose* (bool) - If *True*, print out info

    """

    if fname is None:
        fname = name + '.kml'

    x,y = xy

    if verbose:
        print("Polygon:     %10.6f  %10.6f" % (x[0],y[0]))
        for j in range(1,len(x)):
            print("             %10.6f  %10.6f" % (x[j],y[j]))

    elev = 0.
    kml_text = kml_header(fname)

    mapping = {}
    mapping['x'] = x
    mapping['y'] = y
    mapping['elev'] = elev
    mapping['name'] = name
    d = "  x[0] = %s, y[0] = %s\n" % (x[0],y[0]) 
    for j in range(1,len(x)):
        d = d + "  x[%i] = %s, y[%i] = %s" % (j,f2s(x[j]),j,f2s(y[j]))
    mapping['desc'] = d
    mapping['color'] = color
    mapping['width'] = width

    v = "\n"
    for j in range(len(x)):
        v = v + "%s,%s,%s\n" % (f2s(x[j]),f2s(y[j]),f2s(elev))
    v = v + "%s,%s,%s\n" % (f2s(x[0]),f2s(y[0]),f2s(elev))
    v.replace(' ','')
    
    region_text = kml_region(mapping, v)
    for j in range(1,len(x)):
        d = d + "  x[%i] = %s, y[%i] = %s" % (j,f2s(x[j]),j,f2s(y[j]))

    kml_text = kml_text + region_text + kml_footer()
    kml_file = open(fname,'w')
    kml_file.write(kml_text)
    kml_file.close()
    if verbose:
        print("Created ",fname)


def gauges2kml(rundata=None, fname='gauges.kml', verbose=True):

    """

    Create a KML marker for each gauge specified for a GeoClaw run.

    :Inputs:

      - *rundata* - an object of class *ClawRunData* or None

        If *rundata==None*, try to create based on executing function *setrun*
        from the `setrun.py` file in the current directory.

      - *fname* (str) - resulting kml file.

      - *verbose* (bool) - If *True*, print out info about each region found


    :Example:

        >>> from clawpack.geoclaw import kmltools
        >>> kmltools.gauges2kml()

    is equivalent to:

        >>> from clawpack.geoclaw import kmltools
        >>> from setrun import setrun
        >>> rundata = setrun()
        >>> kmltools.gauges2kml(rundata)

    By default this creates a file named *gauges.kml* that can be opened in
    Google Earth.

    """


    if rundata is None:
        try:
            import setrun
            reload(setrun)
            rundata = setrun.setrun()
        except:
            raise IOError("*** cannot execute setrun file")

    elev = 0.
    kml_text = kml_header(fname)


    gauges = rundata.gaugedata.gauges
    if len(gauges)==0 and verbose:
        print("No gauges found in setrun.py")


    for rnum,gauge in enumerate(gauges):
        t1,t2 = gauge[3:5]
        x1,y1 = gauge[1:3]
        gaugeno = gauge[0]
        if verbose:
            print("Gauge %i: %s, %s  \n" % (gaugeno,f2s(x1),f2s(y1)) \
                    + "  t1 = %s,  t2 = %s" % (f2s(t1),f2s(t2)))
        mapping = {}
        mapping['gaugeno'] = gaugeno
        mapping['t1'] = t1
        mapping['t2'] = t2
        mapping['x1'] = x1
        mapping['y1'] = y1
        mapping['elev'] = elev
        mapping['name'] = 'Gauge %i' % rnum
        description = "  t1 = %s, t2 = %s\n" % (f2s(t1),f2s(t2)) \
            + "  x1 = %s, y1 = %s\n" % (f2s(x1),f2s(y1))
        mapping['desc'] = description

        gauge_text = kml_gauge(mapping)
        kml_text = kml_text + gauge_text
    kml_text = kml_text + kml_footer()
    kml_file = open(fname,'w')
    kml_file.write(kml_text)
    kml_file.close()
    if verbose:
        print("Created ",fname)



def kml_header(name='GeoClaw kml file'):
    header = """<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2">
<Document><name>%s</name>
""" % name
    return header

def kml_footer():
    footer = """
</Document>
</kml>
"""
    return footer


def kml_region(mapping, vertex_text=None):

    if vertex_text is None:
        if 'x3' in mapping:
            # quadrilateral with 4 corners specified
            vertex_text = """
{x1:.9f},{y1:.9f},{elev:.9f}
{x2:.9f},{y2:.9f},{elev:.9f}
{x3:.9f},{y3:.9f},{elev:.9f}
{x4:.9f},{y4:.9f},{elev:.9f}
{x1:.9f},{y1:.9f},{elev:.9f}
""".format(**mapping).replace(' ','')

        else:
            # rectangle with 2 corners specified
            vertex_text = """
{x1:.9f},{y1:.9f},{elev:.9f}
{x2:.9f},{y1:.9f},{elev:.9f}
{x2:.9f},{y2:.9f},{elev:.9f}
{x1:.9f},{y2:.9f},{elev:.9f}
{x1:.9f},{y1:.9f},{elev:.9f}
""".format(**mapping).replace(' ','')

    mapping['vertices'] = vertex_text
    if len(mapping['color'])==6:
        mapping['color'] = 'FF' + mapping['color']

    kml_text = """
<Style id="Path">
<LineStyle><color>{color:s}</color><width>{width:d}</width></LineStyle>
<PolyStyle><color>00000000</color></PolyStyle>
</Style>
<Placemark><name>{name:s}</name>
<description>{desc:s}</description>
<styleUrl>#Path</styleUrl>
<Polygon>
<tessellate>1</tessellate>
<altitudeMode>clampToGround</altitudeMode>
<outerBoundaryIs><LinearRing><coordinates>
{vertices:s}
</coordinates></LinearRing></outerBoundaryIs>
</Polygon>
</Placemark>
""".format(**mapping)

    return kml_text

def kml_line(mapping):

    if len(mapping['color'])==6:
        mapping['color'] = 'FF' + mapping['color']

        line_text = """
{x1:.9f},{y1:.9f},{elev:.9f}
{x2:.9f},{y2:.9f},{elev:.9f}
""".format(**mapping).replace(' ','')

    mapping['line'] = line_text
    kml_text = """
<Style id="Path">
<LineStyle><color>{color:s}</color><width>{width:d}</width></LineStyle>
<PolyStyle><color>00000000</color></PolyStyle>
</Style>
<Placemark><name>{name:s}</name>
<description>{desc:s}</description>
<styleUrl>#Path</styleUrl>
<LineString>
<tessellate>1</tessellate>
<altitudeMode>clampToGround</altitudeMode>
<coordinates>
{line:s}
</coordinates>
</LineString>
</Placemark>
""".format(**mapping)

    return kml_text

def kml_gauge(mapping):
    gauge_text = "{x1:.9f},{y1:.9f},{elev:.9f}".format(**mapping).replace(' ','')

    mapping['gauge'] = gauge_text

    kml_text = """
<Placemark><name>Gauge {gaugeno:d}</name>
<description>{desc:s}</description>
<styleUrl>#markerstyle</styleUrl>
<Point>
<coordinates>
{gauge:s}
</coordinates>
</Point>
</Placemark>
""".format(**mapping)

    return kml_text



def kml_timespan(t1,t2,event_time=None,tz=None,tscale=1):

    r"""
    Create time strings necessary for sliders in Google Earth.  The time
    span will cover time [t1,t2], with the start of the event given by
    event_time.

    [t1,t2]    : time span,

    event_time : Start of event in UTC :  [Y,M,D,H,M,S], e.g. [2010,2,27,3,34,0]
    tz         : time zone offset to UTC.  e.g. +3 for Chile; -9 for Japan.

    Time span element looks like ::

        <TimeSpan>
          <begin>2010-02-27T06:34:00+03:00</begin>
          <end>2010-02-27T07:04:00+03:00</end>
        </TimeSpan>

    As for how well this handles  Daylight  Savings time, here is what the documentation
    on the Python 'time' module has to say :

    "DST is Daylight Saving Time, an adjustment of the timezone by (usually) one hour
    during part of the year. DST rules are magic (determined by local law) and can
    change from year to year. The C library has a table containing the local rules
    (often it is read from a system file for flexibility) and is the only source of
    True Wisdom in this respect."

    """

    t1 = t1*tscale   # Time converted to seconds
    t2 = t2*tscale

    import time
    # to adjust time from UTC to time in event locale.
    if event_time == None:
        # Use local time.
        starttime = time.mktime(time.localtime())  # seconds UTC
        tz_offset = time.timezone/3600.0   # in seconds
    else:
        ev = tuple(event_time) + (0,0,0)   # Extend to 9 tuple; no DST
        # mktime returns time in seconds + timezone offset, i.e. seconds UTC
        # Subtract out the timezone offset here, since it will get added back
        # in when we do gmtime(starttime + ...) below.
        starttime = time.mktime(ev) - time.timezone
        if tz is None:
            print("===> Time zone offset not defined;  assuming zero offset. " \
                "Set plotdata.kml_tz_offset to define an offset (in hours) from "\
                "UTC (positive west of UTC; negative east of UTC)")
            tz = 0

        tz_offset = tz

    if (tz_offset == None):
        tzstr = "Z"  # no offset; could also just set to "+00:00"
    else:
        # Google Earth will show time slider time in local time, where
        # local + offset = UTC.
        tz_offset = tz_offset*3600.    # Offset in seconds
        tz = time.gmtime(abs(tz_offset))
        if (tz_offset > 0):
            tzstr = time.strftime("+%H:%M",tz)  # Time to UTC
        else:
            tzstr = time.strftime("-%H:%M",tz)

    # Get time strings for start and end of time span
    gbegin = time.gmtime(starttime + t1)
    timestrbegin = "%s%s" % (time.strftime("%Y-%m-%dT%H:%M:%S", gbegin),tzstr)

    gend = time.gmtime(starttime + t2)
    timestrend = "%s%s" % (time.strftime("%Y-%m-%dT%H:%M:%S", gend),tzstr)

    return timestrbegin,timestrend

def topo2kml(topo_file_name, topo_type, color='00FF00'):       
    """
    Create a kml file putting a box around the region covered by a topofile.
    Color is green by default.
    """

    import os
    from clawpack.geoclaw import topotools
    topo = topotools.Topography(topo_file_name, topo_type=topo_type)
    topo.read_header()
    xy = topo.extent
    name = os.path.splitext(os.path.split(topo_file_name)[-1])[0]
    file_name = '%s.kml' % name
    box2kml(xy, file_name, name, color)

def dtopo2kml(dtopo_file_name, dtopo_type, color='8888FF'):       
    """
    Create a kml file putting a box around the region covered by a dtopofile.
    Color is pink by default.
    """

    import os
    from clawpack.geoclaw import dtopotools
    dtopo = dtopotools.DTopography()
    dtopo.read(dtopo_file_name, dtopo_type)
    x1 = dtopo.x.min()
    x2 = dtopo.x.max()
    y1 = dtopo.y.min()
    y2 = dtopo.y.max()
    xy = (x1,x2,y1,y2)
    name = os.path.splitext(os.path.split(dtopo_file_name)[-1])[0]
    file_name = '%s.kml' % name
    box2kml(xy, file_name, name, color)
        

def make_input_data_kmls(rundata):
    """
    Produce kml files for the computational domain, all gauges and regions 
    specified, and all topo and dtopo files specified in rundata.
    This can be used, e.g. by adding the lines 

        from clawpack.geoclaw import kmltools
        kmltools.make_input_data_kmls(rundata)

    to the end of a `setrun.py` file so that `make data` will generate all
    kml files in addition to the `*.data` files.
    """
    
    import os
    from . import topotools, dtopotools

    regions2kml(rundata, combined=False)
    gauges2kml(rundata)

    topofiles = rundata.topo_data.topofiles
    for f in topofiles:
        topo_file_name = f[-1]
        topo_type = f[0]
        topo2kml(topo_file_name, topo_type)
        
    dtopofiles = rundata.dtopo_data.dtopofiles
    for f in dtopofiles:
        dtopo_file_name = f[-1]
        dtopo_type = f[0]
        dtopo2kml(dtopo_file_name, dtopo_type)
        
        
