"""
Tools to make kml files to overlay on Google Earth.
Note that color is in KML format, GBR with 2 hex digits for each, e.g.
  FF0000 is blue, 00FF00 is green,  0000FF is red, 00FF00 is yellow.
Actually it's an 8 hex digit number, where the first two digits are
transparency, but in this module these default to 'FF' (but you can specify
the full 8 digits if you want it transparent).

"""


def deg2dms(dy):
    from numpy import floor
    dy_deg = floor(dy)
    dy_min = floor((dy-dy_deg)*60.)
    dy_sec = (dy-dy_deg-dy_min/60.)*3600.
    return dy_deg,dy_min,dy_sec


def regions2kml(rundata=None,fname='regions.kml'):

    """
    Read in the AMR regions from setrun.py and create a kml box for each.
    First create a box for the entire domain (in red) and then a box
    for each region (in white).
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
    description = "  x1 = %g, x2 = %g\n" % (x1,x2) \
            + "  y1 = %g, y2 = %g\n" % (y1,y2) 

    mx,my = clawdata.num_cells[0:]
    dx = (x2-x1)/float(mx)
    dx_meters = dx*111e3*cos(pi*0.5*(y1+y2)/180.)
    dy = (y2-y1)/float(my)
    dy_meters = dy*111e3
    print "Domain:   %10.6f  %10.6f  %10.6f  %10.6f" % (x1,x2,y1,y2)
    dx_deg,dx_min,dx_sec = deg2dms(dx)
    dy_deg,dy_min,dy_sec = deg2dms(dy)
    #print "Level 1 resolution:  dx = %g deg, %g min, %g sec = %g meters" \
    #   % (dx_deg,dx_min,dx_sec,dx_meters)
    levtext = "Level 1 resolution:  dy = %g deg, %g min, %g sec = %g meters\n" \
        % (dy_deg,dy_min,dy_sec,dy_meters)
    print levtext
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
        print levtext
        description = description + levtext

    print "Allowing maximum of %i levels" % amr_levels_max

    elev = 0.
    kml_text = kml_header()
    
    mapping = {}
    mapping['x1'] = x1
    mapping['x2'] = x2
    mapping['y1'] = y1
    mapping['y2'] = y2
    mapping['elev'] = elev
    mapping['name'] = 'Computational Domain'
    mapping['desc'] = description
    mapping['color'] = "0000FF"  # red

    region_text = kml_region(mapping)
    kml_text = kml_text + region_text

    regions = rundata.regiondata.regions
    if len(regions)==0:
        print "No regions found in setrun.py"


    for rnum,region in enumerate(regions):
        minlevel,maxlevel = region[0:2]
        t1,t2 = region[2:4]
        x1,x2,y1,y2 = region[4:]
        print "Region %i: %10.6f  %10.6f  %10.6f  %10.6f" % (rnum,x1,x2,y1,y2)
        print "           minlevel = %i,  maxlevel = %i" % (minlevel,maxlevel) \
                + "  t1 = %10.1f,  t2 = %10.1f" % (t1,t2)
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
            + "  t1 = %g, t2 = %g\n" % (t1,t2) \
            + "  x1 = %g, x2 = %g\n" % (x1,x2) \
            + "  y1 = %g, y2 = %g\n\n" % (y1,y2) 
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

        region_text = kml_region(mapping)
        kml_text = kml_text + region_text
    kml_text = kml_text + kml_footer()
    kml_file = open(fname,'w')
    kml_file.write(kml_text)
    kml_file.close()
    print "Created ",fname
        


def box2kml(xy,fname='box.kml',name='box',color='FF0000'):
    """
    Make a box with default color blue.
    xy should be a tuple (x1,x2,y1,y2).
    """

    x1,x2,y1,y2 = xy[0:]
    print "Box:   %10.6f  %10.6f  %10.6f  %10.6f" % (x1,x2,y1,y2)

    elev = 0.
    kml_text = kml_header() 
    
    mapping = {}
    mapping['x1'] = x1
    mapping['x2'] = x2
    mapping['y1'] = y1
    mapping['y2'] = y2
    mapping['elev'] = elev
    mapping['name'] = name
    mapping['desc'] = "  x1 = %g, x2 = %g\n" % (x1,x2) \
            + "  y1 = %g, y2 = %g" % (y1,y2) 
    mapping['color'] = color

    region_text = kml_region(mapping)

    kml_text = kml_text + region_text + kml_footer()
    kml_file = open(fname,'w')
    kml_file.write(kml_text)
    kml_file.close()
    print "Created ",fname
        

def gauges2kml(rundata=None, fname='gauges.kml'):

    """
    Read in the gauge locations from setrun.py and create a kml point for each.
    """


    if rundata is None:
        try:
            import setrun
            reload(setrun)
            rundata = setrun.setrun()
        except:
            raise IOError("*** cannot execute setrun file")

    elev = 0.
    kml_text = kml_header()
    

    gauges = rundata.gaugedata.gauges
    if len(gauges)==0:
        print "No gauges found in setrun.py"


    for rnum,gauge in enumerate(gauges):
        t1,t2 = gauge[3:5]
        x1,y1 = gauge[1:3]
        gaugeno = gauge[0]
        print "Gauge %i: %10.6f  %10.6f  \n" % (gaugeno,x1,y1) \
                + "  t1 = %10.1f,  t2 = %10.1f" % (t1,t2)
        mapping = {}
        mapping['gaugeno'] = gaugeno
        mapping['t1'] = t1
        mapping['t2'] = t2
        mapping['x1'] = x1
        mapping['y1'] = y1
        mapping['elev'] = elev
        mapping['name'] = 'Gauge %i' % rnum
        description = "  t1 = %g, t2 = %g\n" % (t1,t2) \
            + "  x1 = %g, y1 = %g\n" % (x1,y1)
        mapping['desc'] = description

        gauge_text = kml_gauge(mapping)
        kml_text = kml_text + gauge_text
    kml_text = kml_text + kml_footer()
    kml_file = open(fname,'w')
    kml_file.write(kml_text)
    kml_file.close()
    print "Created ",fname
        


def kml_header():
    """
    Color is a BGR hex string used to set color of lines.
    default is red.
    """
    header = """<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2">
<Document><name>My document</name>
<description>Content</description>
""" 
    return header

def kml_footer():
    footer = """
</Document>
</kml>
"""
    return footer

def kml_region(mapping):
    region_text = """
{x1:10.4f},{y1:10.4f},{elev:10.4f}
{x2:10.4f},{y1:10.4f},{elev:10.4f}
{x2:10.4f},{y2:10.4f},{elev:10.4f}
{x1:10.4f},{y2:10.4f},{elev:10.4f}
{x1:10.4f},{y1:10.4f},{elev:10.4f}
""".format(**mapping).replace(' ','')
    
    mapping['region'] = region_text
    if len(mapping['color'])==6: 
        mapping['color'] = 'FF' + mapping['color']

    kml_text = """
<Style id="Path">
<LineStyle><color>{color:s}</color><width>3</width></LineStyle>
<PolyStyle><color>00000000</color></PolyStyle>
</Style>
<Placemark><name>{name:s}</name>
<description>{desc:s}</description>
<styleUrl>#Path</styleUrl>
<Polygon>
<tessellate>1</tessellate>
<altitudeMode>clampToGround</altitudeMode>
<outerBoundaryIs><LinearRing><coordinates>
{region:s}
</coordinates></LinearRing></outerBoundaryIs>
</Polygon>
</Placemark>
""".format(**mapping)

    return kml_text

def kml_gauge(mapping):
    gauge_text = "{x1:10.4f},{y1:10.4f},{elev:10.4f}".format(**mapping).replace(' ','')
    
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
