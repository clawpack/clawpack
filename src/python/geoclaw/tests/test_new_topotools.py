
import numpy
from clawpack.geoclaw import topotools 

def topo_bowl(x,y):
    """Sample topo"""
    z = 1000.*(x**2 + y**2 - 1.)
    return z


def test_read_write_topo_bowl():
    """
    Test writing and reading topo files with small number of points
    Note that ordering should go from NW corner.
    """

    nxpoints = 5
    nypoints = 4
    xlower = -1.
    xupper = 3.
    ylower = 0.
    yupper = 3.
    topo = topotools.Topography(topo_func=topo_bowl)
    topo.x = numpy.linspace(xlower,xupper,nxpoints)
    topo.y = numpy.linspace(ylower,yupper,nypoints)

    for ttype in [1,2,3]:
        fname = 'bowl.tt%s' % ttype
        if ttype==1:
            topotools.topo1writer(fname,topo_bowl,xlower,xupper,ylower,yupper,\
                              nxpoints,nypoints)
        elif ttype==2:
            topotools.topo2writer(fname,topo_bowl,xlower,xupper,ylower,yupper,\
                              nxpoints,nypoints)
        elif ttype==3:
            topotools.topo3writer(fname,topo_bowl,xlower,xupper,ylower,yupper,\
                              nxpoints,nypoints)
        print "Created ",fname
        topo_in = topotools.Topography(fname)
        print "Read back in and max difference in z is ",abs(topo.Z - topo_in.Z).max()
        
        topo_in = topotools.Topography()
        topo_in.read(path=fname)  # should figure out topo_type properly
        print "Read back in second way and max difference in z is ",abs(topo.Z - topo_in.Z).max()


def test_against_old():
    """
    Test against the old topotools from 5.1.0.
    Compare bowl.tt1 to bowl_old.tt1
    """

    try:
        from clawpack.geoclaw import topotools_5_1_0 as old_topotools
    except:
        print "*** can't compare to old: need topotools_5_1_0.py"
        return

    nxpoints = 5
    nypoints = 4
    xlower = -1.
    xupper = 3.
    ylower = 0.
    yupper = 3.
    topo = topotools.Topography(topo_func=topo_bowl)
    topo.x = numpy.linspace(xlower,xupper,nxpoints)
    topo.y = numpy.linspace(ylower,yupper,nypoints)

    fname = 'bowl_old.tt1'
    old_topotools.topo1writer(fname,topo_bowl,xlower,xupper,ylower,yupper,\
                      nxpoints,nypoints)
    print "Created ",fname
    X,Y,Z = old_topotools.topofile2griddata(fname,topotype=1)
    # flip because new topotools uses different convention:
    Y = numpy.flipud(Y)
    Z = numpy.flipud(Z)
    print "Read back in and max difference in x is ",abs(topo.X - X).max()
    print "Read back in and max difference in y is ",abs(topo.Y - Y).max()
    print "Read back in and max difference in z is ",abs(topo.Z - Z).max()


def topo_bowl_hill(x,y):
    """
    Sample topography
    """
    # Parabolic bowl
    z = 1000.*(x**2 + y**2 - 1.)
    # Add a Gaussian hill
    z = z + 1000.*numpy.exp(-100*((x-0.7)**2 + (y-0.8)**2))
    return z

def make_topo_bowl_hill():
    """
    Create topo file as a sample to read and plot.
    """

    nxpoints = 101
    nypoints = 76
    xlower = -1.5
    xupper = 2.5
    ylower = -1.
    yupper = 2.

    fname = 'bowl_hill.tt2'
    topotools.topo2writer(fname,topo_bowl_hill,xlower,xupper,ylower,yupper,\
                          nxpoints,nypoints)
    print "Created ",fname


def test_read_write_topo_bowl_hill():
    """
    Test writing and reading topo files.
    """

    nxpoints = 101
    nypoints = 76
    xlower = -1.5
    xupper = 2.5
    ylower = -1.
    yupper = 2.
    topo = topotools.Topography(topo_func=topo_bowl_hill)
    topo.x = numpy.linspace(xlower,xupper,nxpoints)
    topo.y = numpy.linspace(ylower,yupper,nypoints)

    for ttype in [1,2,3]:
        fname = 'bowl_hill.tt%s' % ttype
        topo.write(fname, topo_type=ttype)
        print "Created ",fname
        topo_in = topotools.Topography(fname,topo_type=ttype)
        print "Read back in and max difference in z is ",abs(topo.Z - topo_in.Z).max()


def test_plot_topo_bowl_hill():

    """
    Create topo and write out, then read in again and plot.
    Note that center of bowl should be at (0,0).
    """

    import matplotlib.pyplot as plt
    make_topo_bowl_hill()
    fname = 'bowl_hill.tt2'
    topo = topotools.Topography(fname,topo_type=2)

    topo.plot()

    topo2 = topo.crop([0.5, 1.5, 0., 2.])
    topo2.plot()
    plt.title("Cropped topography")

