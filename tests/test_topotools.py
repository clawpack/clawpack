#!/usr/bin/env python

from __future__ import absolute_import
from __future__ import print_function
import os
import sys
import tempfile
import shutil

try:
    # For Python 3.0 and later
    from urllib.error import URLError
except ImportError:
    # Fall back to Python 2's urllib2
    from urllib2 import URLError

import numpy

import nose

import clawpack.geoclaw.topotools as topotools
import clawpack.clawutil.data
from six.moves import range

# Set local test directory to get local files
testdir = os.path.dirname(__file__)
if len(testdir) == 0:
     testdir = "./"


def topo_bowl(x,y):
    """Sample topo"""
    z = 1000.*(x**2 + y**2 - 1.)
    return z


def test_read_write_topo_bowl():
    """
    Test writing and reading topo files with small number of points
    Note that ordering should go from NW corner.
    """

    # Base topography
    topo = topotools.Topography(topo_func=topo_bowl)
    topo.x = numpy.linspace(-1.0, 3.0, 5)
    topo.y = numpy.linspace( 0.0, 3.0, 4)

    assert numpy.allclose(topo.x, numpy.array([-1.,  0.,  1.,  2., 3.])), \
           "Topography x values are incorrect."
    assert numpy.allclose(topo.X,
                          numpy.array([[-1.,  0.,  1.,  2.,  3.],
                                       [-1.,  0.,  1.,  2.,  3.],
                                       [-1.,  0.,  1.,  2.,  3.],
                                       [-1.,  0.,  1.,  2.,  3.]])), \
           "Topography X values are incorrect."
    assert numpy.allclose(topo.y, numpy.array([ 0.,  1.,  2.,  3.])), \
           "Topography y values are incorrect."
    assert numpy.allclose(topo.Y,
                          numpy.array([[ 0.,  0.,  0.,  0.,  0.],
                                       [ 1.,  1.,  1.,  1.,  1.],
                                       [ 2.,  2.,  2.,  2.,  2.],
                                       [ 3.,  3.,  3.,  3.,  3.]])), \
           "Topography Y values are incorrect."
    assert numpy.allclose(topo.Z,
                numpy.array([[     0.,  -1000.,      0.,   3000., 8000.],
                             [  1000.,      0.,   1000.,   4000.,   9000.],
                             [  4000.,   3000.,   4000.,   7000.,  12000.],
                             [  9000.,   8000.,   9000.,  12000.,  17000.]])), \
           "Topography Z values are incorrect."

    temp_path = tempfile.mkdtemp()
    try:
        for topo_type in range(1, 4):
            path = os.path.join(temp_path, 'bowl.tt%s' % topo_type)
            topo.write(path, topo_type=topo_type,Z_format="%22.15e")

            topo_in = topotools.Topography(path)
            assert numpy.allclose(topo.Z, topo_in.Z), \
                   "Differnece in written and read topography found."
    except AssertionError as e:
        # If the assertion failed then copy the contents of the directory
        shutil.copytree(temp_path, os.path.join(os.getcwd(), 
                                                "test_read_write_topo_bowl"))
        raise e
    finally:
        shutil.rmtree(temp_path)


def test_crop_topo_bowl():
    """
    Test cropping a topo file.
    """

    topo = topotools.Topography(topo_func=topo_bowl)
    topo.x = numpy.linspace(-1.0, 3.0, 5)
    topo.y = numpy.linspace( 0.0, 3.0, 4)

    # topo.Z should be created automatically when referenced below:
    assert numpy.allclose(topo.Z,
           numpy.array([[     0.,  -1000.,      0.,   3000., 8000.],
                        [  1000.,      0.,   1000.,   4000.,   9000.],
                        [  4000.,   3000.,   4000.,   7000.,  12000.],
                        [  9000.,   8000.,   9000.,  12000.,  17000.]])), \
           "Basic topography does not match test data."

    cropped_topo = topo.crop([0, 1, 0, 2])
    assert numpy.allclose(cropped_topo.x, numpy.array([0.0, 1.0])), \
           "Cropped topography y values do not match"
    assert numpy.allclose(cropped_topo.y, numpy.array([ 0.,  1.,  2.])), \
           "Cropped topography y values do not match."
    assert numpy.allclose(cropped_topo.Z,
                          numpy.array([[-1000.,     0.],
                                       [    0.,  1000.],
                                       [ 3000.,  4000.]])), \
           "Cropped topography Z values do not match."



def test_against_old():
    """
    Test against the old topotools from 5.1.0.
    Compare bowl.tt1 to bowl_old.tt1
    """
    
    from . import old_topotools

    nxpoints = 5
    nypoints = 4
    xlower = -1.0
    xupper = 3.0
    ylower = 0.0
    yupper = 3.0
    topo = topotools.Topography(topo_func=topo_bowl)
    topo.x = numpy.linspace(xlower, xupper, nxpoints)
    topo.y = numpy.linspace(ylower, yupper, nypoints)

    temp_path = tempfile.mkdtemp()
    try:
        file_path = os.path.join(temp_path, "bowl_old.tt1")
        old_topotools.topo1writer(file_path, topo_bowl, xlower, xupper, ylower, 
                                  yupper, nxpoints, nypoints)
        X, Y, Z = old_topotools.topofile2griddata(file_path, topotype=1)
        Y = numpy.flipud(Y)
        Z = numpy.flipud(Z)

        assert numpy.allclose(topo.X, X), "Difference in X grid."
        assert numpy.allclose(topo.Y, Y), "Difference in Y grid."
        assert numpy.allclose(topo.Z, Z), "Difference in Z grid."

    except AssertionError as e:
        # If the assertion failed then copy the contents of the directory
        shutil.copytree(temp_path, os.path.join(os.getcwd(), 
                                                "test_against_old"))
        raise e
    finally:
        shutil.rmtree(temp_path)



def topo_bowl_hill(x,y):
    """
    Sample topography
    """
    # Parabolic bowl
    z = 1000.*(x**2 + y**2 - 1.)
    # Add a Gaussian hill
    z = z + 1000.*numpy.exp(-100*((x-0.7)**2 + (y-0.8)**2))
    return z


def test_read_write_topo_bowl_hill():
    """
    Test writing and reading topo files.
    """
    temp_path = tempfile.mkdtemp()

    try:
        topo = topotools.Topography(topo_func=topo_bowl_hill)
        topo.x = numpy.linspace(-1.5, 2.5, 101)
        topo.y = numpy.linspace(-1.0, 2.0, 76)

        for topo_type in range(1,4):
            file_path = os.path.join(temp_path, 'bowl_hill.tt%s' % topo_type)
            topo.write(file_path, topo_type=topo_type,Z_format="%22.15e")
            topo_in = topotools.Topography(path=file_path, topo_type=topo_type)
            assert numpy.allclose(topo.Z, topo_in.Z), \
                   "Written file of topo_type=%s does not equal read in" + \
                   " file." % topo_type

    except AssertionError as e:
        # If the assertion failed then copy the contents of the directory
        shutil.copytree(temp_path, os.path.join(os.getcwd(), 
                                              "test_read_write_topo_bowl_hill"))
        raise e
    finally:
        shutil.rmtree(temp_path)


def test_netcdf():
    r"""Test Python NetCDF formatted topography reading"""

    temp_path = tempfile.mkdtemp()

    try:
        # Fetch comparison data
        url = "".join(('https://raw.githubusercontent.com/rjleveque/geoclaw/',
                       '5f675256c043e59e5065f9f3b5bdd41c2901702c/src/python/',
                       'geoclaw/tests/kahului_sample_1s.tt2'))
        clawpack.clawutil.data.get_remote_file(url, output_dir=temp_path,
                                                    force=True)
        
        # Paths
        local_path = os.path.join(temp_path, os.path.basename(url))
        nc_path = os.path.join(temp_path, "test.nc")

        # Write out NetCDF version of file
        ascii_topo = topotools.Topography(path=local_path)
        ascii_topo.read()
        ascii_topo.write(nc_path, topo_type=4,Z_format="%22.15e")

        # Read back in NetCDF file
        nc_topo = topotools.Topography(path=nc_path)
        nc_topo.read()

        # Compare arrays - use tolerance based on 30 arcsecond accuracy
        assert numpy.allclose(ascii_topo.x, nc_topo.x), \
                    "Flat x-arrays did not match."
        assert numpy.allclose(ascii_topo.y, nc_topo.y), \
                    "Flat y-arrays did not match."
        assert numpy.allclose(ascii_topo.Z, nc_topo.Z), \
                    "Flat y-arrays did not match."

    except AssertionError as e:
        shutil.copytree(temp_path, os.path.join(os.getcwd()),
            'test_read_netcdf')
        raise e

    except ImportError as e:
        raise nose.SkipTest("Skipping test since NetCDF support not found.")

    except RuntimeError as e:
        raise nose.SkipTest("NetCDF topography test skipped due to " +
                            "runtime failure.")
    except URLError:
        raise nose.SkipTest("Could not fetch remote file, skipping test.")
    
    finally:
        shutil.rmtree(temp_path)


def test_get_remote_file():
    """Test the ability to fetch a remote file from the web."""
    
    temp_path = tempfile.mkdtemp()
    try:

        url = "".join(('https://raw.githubusercontent.com/rjleveque/geoclaw/',
                       '5f675256c043e59e5065f9f3b5bdd41c2901702c/src/python/',
                       'geoclaw/tests/kahului_sample_1s.tt2'))
        clawpack.clawutil.data.get_remote_file(url, output_dir=temp_path,
            force=True)

        local_path = os.path.join(temp_path, os.path.basename(url))
        download_topo = topotools.Topography(path=local_path)

        test_path = os.path.join(testdir, "data", os.path.basename(url))
        test_topo = topotools.Topography(path=test_path)

        assert numpy.allclose(download_topo.Z, test_topo.Z), \
               "Downloaded file does not match %s" % test_path
    except AssertionError as e:
        shutil.copy(local_path, os.path.join(os.getcwd(), "remote_file.tt2"))
        raise e

    except URLError:
        raise nose.SkipTest("Could not fetch remote file, skipping test.")

    finally:
        shutil.rmtree(temp_path)


def test_unstructured_topo(save=False, plot=False):

    try:
        import scipy
    except:
        raise nose.SkipTest("Skipping test since scipy not found")

    # Create random test data
    def func(x, y):
        return x * (1 - x) * numpy.cos(4 * numpy.pi * x) * numpy.sin(4 * numpy.pi * y**2)**2

    fill_topo = topotools.Topography()
    fill_topo.x = numpy.linspace(0, 1, 100)
    fill_topo.y = numpy.linspace(0, 1, 200)
    fill_topo.Z = func(fill_topo.X, fill_topo.Y)

    points = numpy.loadtxt(os.path.join(testdir, "data", 
                                                 "unstructured_points.txt"))
    values = func(points[:,0], points[:,1])

    # Create topography object
    topo = topotools.Topography(unstructured=True)
    topo.x = points[:,0]
    topo.y = points[:,1]
    topo.z = values

    if plot:
        import matplotlib.pyplot as plt

        fig = plt.figure(figsize=(16,6))
        axes = fig.add_subplot(1, 3, 1)
        fill_topo.plot(axes=axes)
        axes.set_title("True Field")
        axes = fig.add_subplot(1, 3, 2)
        topo.plot(axes=axes, region_extent=[0, 1, 0, 1])
        axes.set_title("Unstructured Field")

    topo.interp_unstructured(fill_topo, extent=[0, 1, 0, 1], delta=(1e-2,1e-2))
    assert not topo.unstructured

    # Load (and save) test data and make the comparison
    test_data_path = os.path.join(testdir, "data", "unstructured_test_data.tt3")
    if save:
        topo.write(test_data_path,Z_format="%22.15e")

    compare_data = topotools.Topography(path=test_data_path)

    assert numpy.allclose(compare_data.Z, topo.Z)

    if plot:
        axes = fig.add_subplot(1, 3, 3)
        topo.plot(axes=axes)
        axes.set_title("Interpolated Field")

        plt.show()


def plot_topo_bowl_hill():

    """
    Create topo and write out, then read in again and plot.
    Note that center of bowl should be at (0,0).
    """

    try:
        import matplotlib
    except ImportError:
        raise nose.SkipTest("Skipping test since matplotlib not found.")

    matplotlib.use("Agg")  # use image backend -- needed for Travis tests
    import matplotlib.pyplot as plt

    topo = topotools.Topography(topo_func=topo_bowl_hill)
    topo.x = numpy.linspace(-1.5, 2.5, 101)
    topo.y = numpy.linspace(-1.0, 2.0, 76)
    fname = 'bowl_hill.tt2'
    topo = topotools.Topography(fname,topo_type=2)

    topo.plot()
    fname = "bowl_hill.png"
    plt.savefig(fname)
    print("Created ",fname)

    topo2 = topo.crop([0.5, 1.5, 0., 2.])
    topo2.plot()
    plt.title("Cropped topography")
    fname = "bowl_hill_crop.png"
    plt.savefig(fname)
    print("Created ",fname)


def plot_kahului():
    r"""
    Example illustrating reading in a topo file and plotting.
    Uses the test data kahului_sample_1s.tt2, created by cropping 
    the data file obtained from the NGDC site
        http://www.ngdc.noaa.gov/dem/squareCellGrid/download/604
    In addition to using the Topography.plot function, also 
    illustrate how to do a contour data of the data directly.
    """

    try:
        import matplotlib
    except ImportError:
        raise nose.SkipTest("Skipping test since matplotlib not found.")

    matplotlib.use("Agg")  # use image backend -- needed for Travis tests
    import matplotlib.pyplot as plt

    path = os.path.join(testdir,'kahului_sample_1s.tt2')
    K = topotools.Topography(path,topo_type=2)
    K.plot()
    plt.title("Kahului Harbor at 1 second resolution")

    plt.title("Kahului Harbor at 1 second resolution")
    fname = "kahului_imshow.png"
    plt.savefig(fname)
    print("Created ",fname)

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

    # fix aspect ratio based on latitude:
    mean_lat = 0.5 * (K.y.max() + K.y.min())
    ax.set_aspect(1.0 / numpy.cos(numpy.pi / 180.0 * mean_lat))

    # fix tick marks so readable:
    ax.ticklabel_format(format="plain", useOffset=False)
    plt.xticks(rotation=20)

    plt.title("2-meter contours of topo (green) and bathymetry (blue)",\
              fontsize=12)
    fname = "kahului_contour.png"
    plt.savefig(fname)
    print("Created ",fname)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        if "plot" in sys.argv[1].lower():
            plot_kahului()
            plot_topo_bowl_hill()
            test_unstructured_topo(save=False, plot=True)
        elif bool(sys.argv[1]):
            test_unstructured_topo(save=True)
    else:
        # Run tests one at a time
        test_read_write_topo_bowl()
        test_crop_topo_bowl()
        test_against_old()
        test_read_write_topo_bowl_hill()
        test_get_remote_file()
        test_unstructured_topo()
        test_netcdf()

        print("All tests passed.")
