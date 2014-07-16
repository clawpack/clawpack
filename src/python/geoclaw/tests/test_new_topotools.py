
import os
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

    assert numpy.allclose(topo.x, numpy.array([-1.,  0.,  1.,  2., 3.])), \
        "*** topo.x does not match"
    assert numpy.allclose(topo.X, \
        numpy.array([[-1.,  0.,  1.,  2.,  3.],
                   [-1.,  0.,  1.,  2.,  3.],
                   [-1.,  0.,  1.,  2.,  3.],
                   [-1.,  0.,  1.,  2.,  3.]])), \
        "*** topo.X does not match"
    assert numpy.allclose(topo.y, numpy.array([ 0.,  1.,  2.,  3.])), \
        "*** topo.y does not match"
    assert numpy.allclose(topo.Y, \
        numpy.array([[ 0.,  0.,  0.,  0.,  0.],
                   [ 1.,  1.,  1.,  1.,  1.],
                   [ 2.,  2.,  2.,  2.,  2.],
                   [ 3.,  3.,  3.,  3.,  3.]])), \
        "*** topo.Y does not match"
    assert numpy.allclose(topo.Z, \
        numpy.array([[     0.,  -1000.,      0.,   3000., 8000.],
                   [  1000.,      0.,   1000.,   4000.,   9000.],
                   [  4000.,   3000.,   4000.,   7000.,  12000.],
                   [  9000.,   8000.,   9000.,  12000.,  17000.]])), \
        "*** topo.Z does not match"

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
    return topo

def test_crop_topo_bowl():
    """
    Test cropping a topo file.
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

    # topo.Z should be created automatically when referenced below:
    assert numpy.allclose(topo.Z, \
        numpy.array([[     0.,  -1000.,      0.,   3000., 8000.],
                   [  1000.,      0.,   1000.,   4000.,   9000.],
                   [  4000.,   3000.,   4000.,   7000.,  12000.],
                   [  9000.,   8000.,   9000.,  12000.,  17000.]])), \
        "*** topo.Z does not match"

    topo2 = topo.crop([0,1,0,2])

    assert numpy.allclose(topo2.x, numpy.array([0.,  1.])), \
        "*** topo.x does not match"
    assert numpy.allclose(topo2.y, numpy.array([ 0.,  1.,  2.])), \
        "*** topo.y does not match"
    assert numpy.allclose(topo2.Z, \
        numpy.array([[-1000.,     0.],
                   [    0.,  1000.],
                   [ 3000.,  4000.]])), \
        "*** topo.Z does not match"

    return topo2


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


def test_plot_kahului():
    r"""
    Example illustrating reading in a topo file and plotting.
    Uses the test data kahului_sample_1s.tt2, created by cropping 
    the data file obtained from the NGDC site
        http://www.ngdc.noaa.gov/dem/squareCellGrid/download/604
    In addition to using the Topography.plot function, also 
    illustrate how to do a contour data of the data directly.
    """
    import matplotlib.pyplot as plt

    K = topotools.Topography('kahului_sample_1s.tt2',topo_type=2)
    K.plot()
    plt.title("Kahului Harbor at 1 second resolution")

    plt.title("Kahului Harbor at 1 second resolution")

    assert K.Z.shape == (46, 65), "*** K.Z is wrong shape"
    assert numpy.allclose(K.Z[:3,:3], \
                          numpy.array([[ 11.339,  11.339, 11.339],
                                       [ 13.339,  11.339,  11.339],
                                       [ 13.339,  11.339, 10.339]])), \
                "*** Topography K does not match"

    # Make a contour plot of topography / bathymetry:
    plt.figure()
    ax = plt.axes()
    plt.contour(K.X, K.Y, K.Z, numpy.linspace(-20,-2,10), colors='b', \
                linestyles='-')
    plt.contour(K.X, K.Y, K.Z, numpy.linspace(2,20,10), colors='g')
    plt.contour(K.X, K.Y, K.Z, [0.], colors='r')  # mean high water
    ax.ticklabel_format(format="plain", useOffset=False)
    plt.title("2-meter contours of topo (green) and bathymetry (blue)",\
              fontsize=12)

def test_fetch_topo_url():

    """
    Fetch topography file from the web.
    """

    K = topotools.Topography('kahului_sample_1s.tt2',topo_type=2)
    K.read()

    url = 'https://raw.githubusercontent.com/rjleveque/geoclaw/5f675256c043e59e5065f9f3b5bdd41c2901702c/src/python/geoclaw/tests/kahului_sample_1s.tt2'
    fname = 'kahului_web.tt2'
    topotools.fetch_topo_url(url, local_fname=fname, force=True)
    K2 = topotools.Topography(fname)
    K2.read()
    
    assert numpy.allclose(K.Z, K2.Z), "*** K2 != K after fetching from URL"
    os.system("rm %s" % fname)
    os.system("rm %s.txt" % fname)

