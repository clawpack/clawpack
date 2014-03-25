"""
Tools to make kml files to overlay on Google Earth.
Note that color is in KML format, GBR with 2 hex digits for each, e.g.
  FF0000 is blue, 00FF00 is green,  0000FF is red, 00FF00 is yellow.
Actually it's an 8 hex digit number, where the first two digits are
transparency, but in this module these default to 'FF' (but you can specify
the full 8 digits if you want it transparent).

"""

def regions2kml(fname='regions.kml'):

    """
    Read in the AMR regions from setrun.py and create a kml box for each.
    First create a box for the entire domain (in red) and then a box
    for each region (in white).
    """

    from pylab import text
    try:
        import setrun
        reload(setrun)
        rundata = setrun.setrun()
    except:
        raise IOError("*** cannot execute setrun file")

    clawdata = rundata.clawdata
    x1,y1 = clawdata.lower[0:]
    x2,y2 = clawdata.upper[0:]
    print "Domain:   %10.6f  %10.6f  %10.6f  %10.6f" % (x1,x2,y1,y2)


    elev = 0.
    kml_text = kml_header()
    
    mapping = {}
    mapping['x1'] = x1
    mapping['x2'] = x2
    mapping['y1'] = y1
    mapping['y2'] = y2
    mapping['elev'] = elev
    mapping['name'] = 'Computational Domain'
    mapping['desc'] = "  x1 = %g, x2 = %g\n" % (x1,x2) \
            + "  y1 = %g, y2 = %g" % (y1,y2) 
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
        mapping['desc'] = "minlevel = %i, maxlevel = %i\n" % (minlevel,maxlevel) \
            + "  t1 = %g, t2 = %g\n" % (t1,t2) \
            + "  x1 = %g, x2 = %g\n" % (x1,x2) \
            + "  y1 = %g, y2 = %g" % (y1,y2) 
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
