#!/usr/bin/env python

import os
import sys
import shutil
import tempfile
import inspect

import numpy
import nose

import clawpack.geoclaw.topotools as topotools
import clawpack.geoclaw.dtopotools as dtopotools

# Set local test directory to get local files
testdir = os.path.dirname(__file__)
if len(testdir) == 0:
     testdir = "./"

def test_read_csv_make_dtopo(save=False):
    r"""Test reading and making of a CSV subfault speficied dtopo."""
    
    subfault_path = os.path.join(testdir, 'data', 'alaska1964.csv')
    input_units = {"length":"km", "width":"km", "depth":"km", "slip":"m", 
             "mu":"dyne/cm^2"}
    fault = dtopotools.CSVFault()
    fault.read(subfault_path, input_units=input_units, 
                              coordinate_specification="noaa sift")

    assert abs(fault.Mw() - 8.53336) < 1e-4, "*** Mw is wrong: %g" % fault.Mw()

    xlower = 203
    xupper = 214.  # approximate - adjusted below
    ylower = 55
    yupper = 60.  # approximate - adjusted below

    # dtopo parameters:
    points_per_degree = 4  # 15 minute resolution
    dx = 1. / points_per_degree
    mx = int((xupper - xlower) / dx + 1)
    xupper = xlower + (mx - 1) * dx
    my = int((yupper - ylower) / dx + 1)
    yupper = ylower + (my - 1) * dx

    x = numpy.linspace(xlower, xupper, mx)
    y = numpy.linspace(ylower, yupper, my)

    dtopo = fault.create_dtopography(x, y, times=[1.])
    # print "max dz = ",dtopo.dz_list[-1].max()
    
    # temp_path = tempfile.mkdtemp()
    # try:
    test_data_path = os.path.join(testdir, "data", "alaska1964_test_data.tt3")
    if save:
        dtopo.write(test_data_path, dtopo_type=3)
    compare_data = dtopotools.DTopography(path=test_data_path)
    compare_data.read(path=test_data_path, dtopo_type=3)

    assert dtopo.DZ.shape == compare_data.DZ.shape, \
        "dtopo.DZ.shape is %s, should be %s" \
        % (dtopo.DZ.shape, compare_data.DZ.shape)

    assert numpy.allclose(compare_data.DZ, dtopo.DZ)

    # except AssertionError as e:
    #     test_name = inspect.stack()[1][-2][0][:-3]
    #     test_dump_path = os.path.join(os.getcwd(), test_name)
    #     shutil.mkdir(test_dump_path)
    #     shutil.copy(temp_path, test_dump_path)
    #     raise e
    # finally:
    #     shutil.rmtree(temp_path)


def test_read_ucsb_make_dtopo(save=False):
    r"""Test reading and making of a UCSB subfault speficied dtopo."""

    subfault_path = os.path.join(testdir, 'data', 'tohoku_ucsb.txt')
    fault = dtopotools.UCSBFault()
    fault.read(subfault_path)

    assert abs(fault.Mw() - 9.13957973) < 1e-4, "*** Mw is wrong: %g" % fault.Mw()

    xlower = 140.
    xupper = 146.
    ylower = 35.
    yupper = 41.

    # dtopo parameters:
    points_per_degree = 4  # 15 minute resolution
    dx = 1. / points_per_degree
    mx = int((xupper - xlower) / dx + 1)
    xupper = xlower + (mx - 1) * dx
    my = int((yupper - ylower) / dx + 1)
    yupper = ylower + (my - 1) * dx

    x = numpy.linspace(xlower, xupper, mx)
    y = numpy.linspace(ylower, yupper, my)

    tmax = 0.
    for s in fault.subfaults:
        tmax = max(tmax, s.rupture_time + s.rise_time + s.rise_time_ending)

    fault.rupture_type = 'dynamic'
    times = numpy.linspace(0, tmax, 10)
    dtopo = fault.create_dtopography(x, y, times)
    # print "max dz = ",dtopo.dz_list[-1].max()
    
    # temp_path = tempfile.mkdtemp()
    # try:
    # Load (and save) test data and make the comparison
    test_data_path = os.path.join(testdir, "data", "tohoku_test_data.tt3")
    if save:
        dtopo.write(test_data_path, dtopo_type=3)
    compare_data = dtopotools.DTopography(path=test_data_path)
    compare_data.read(path=test_data_path, dtopo_type=3)

    assert dtopo.DZ.shape == compare_data.DZ.shape, \
        "dtopo.DZ.shape is %s, should be %s" \
        % (dtopo.DZ.shape, compare_data.DZ.shape)

    assert numpy.allclose(compare_data.DZ, dtopo.DZ)

    # except AssertionError as e:
    #     shutil.copy(test_data_path, os.path.split(test_data_path)[1])
    #     raise e
    # finally:
    #     shutil.rmtree(temp_path)


def test_read_sift_make_dtopo(save=False):
    r"""Test reading and making of a SIFT subfault speficied dtopo"""

    sift_slip = {'acsza1':2, 'acszb1':3}
    fault = dtopotools.SiftFault(sift_slip)

    assert abs(fault.Mw() - 7.3) < 1e-4, "*** Mw is wrong: %g" % fault.Mw()

    xlower = 162.
    xupper = 168.
    ylower = 53.
    yupper = 59.

    # dtopo parameters:
    points_per_degree = 4  # 15 minute resolution
    dx = 1./points_per_degree
    mx = int((xupper - xlower)/dx + 1)
    xupper = xlower + (mx-1)*dx
    my = int((yupper - ylower)/dx + 1)
    yupper = ylower + (my-1)*dx

    x = numpy.linspace(xlower,xupper,mx)
    y = numpy.linspace(ylower,yupper,my)

    times = [1.]
    dtopo = fault.create_dtopography(x,y,times)
    # print "max dz = ",dtopo.dz_list[-1].smax()
    
    # temp_path = tempfile.mkdtemp()
    # try:
    # Load (and save) test data and make the comparison
    test_data_path = os.path.join(testdir, "data", "sift_test_data.tt3")
    if save:
        dtopo.write(test_data_path, dtopo_type=3)
    compare_data = dtopotools.DTopography(path=test_data_path)
    compare_data.read(path=test_data_path, dtopo_type=3)

    assert dtopo.DZ.shape == compare_data.DZ.shape, \
        "dtopo.DZ.shape is %s, should be %s" \
        % (dtopo.DZ.shape, compare_data.DZ.shape)

    assert numpy.allclose(compare_data.DZ, dtopo.DZ)

    # except AssertionError as e:
    #     test_name = inspect.stack()[1][-2][0][:-3]
    #     test_dump_path = os.path.join(os.getcwd(), test_name)
    #     shutil.mkdir(test_dump_path)
    #     shutil.copy(temp_path, test_dump_path)
    #     raise e
    # finally:
    #     shutil.rmtree(temp_path)


def test_SubdividedPlaneFault_make_dtopo(save=False):
    r""""""

    # get a unit source fault plane as starting point:
    sift_slip = {'acsza1':1.}
    fault = dtopotools.SiftFault(sift_slip)
    fault_plane = fault.subfaults[0]
    # Mo = fault_plane.Mo()
    # print "original Mo = ",Mo

    fault2 = dtopotools.SubdividedPlaneFault(fault_plane, nstrike=5, ndip=3)
    # print "new Mo = ",fault2.Mo()
    #fault2.plot_subfaults(slip_color=True)

    assert abs(fault2.Mw() - 6.83402) < 1e-4, \
           "*** Mw is wrong: %g" % fault.Mw()

    xlower = 162.
    xupper = 168.
    ylower = 53.
    yupper = 59.

    # dtopo parameters:
    points_per_degree = 4  # 15 minute resolution
    dx = 1./points_per_degree
    mx = int((xupper - xlower)/dx + 1)
    xupper = xlower + (mx-1)*dx
    my = int((yupper - ylower)/dx + 1)
    yupper = ylower + (my-1)*dx

    x = numpy.linspace(xlower,xupper,mx)
    y = numpy.linspace(ylower,yupper,my)

    times = [1.]
    dtopo = fault2.create_dtopography(x,y,times)
    # print "max dz = ", dtopo.dz_list[-1].max()
    
    # temp_path = tempfile.mkdtemp()
    # try:
    # Load (and save) test data and make the comparison
    test_data_path = os.path.join(testdir, "data", 
            "SubdividedFaultPlane_test_data.tt3")
    if save:
        dtopo.write(test_data_path, dtopo_type=3)
    compare_data = dtopotools.DTopography(path=test_data_path)
    compare_data.read(path=test_data_path, dtopo_type=3)

    assert dtopo.DZ.shape == compare_data.DZ.shape, \
        "dtopo.DZ.shape is %s, should be %s" \
        % (dtopo.DZ.shape, compare_data.DZ.shape)

    assert numpy.allclose(compare_data.DZ, dtopo.DZ)

    # except AssertionError as e:
    #     shutil.copy(test_data_path, os.path.split(test_data_path)[1])
    #     raise e
    # finally:
    #     shutil.rmtree(temp_path)


def test_dtopo_io():
    r"""Test IO of dtopography class"""

    test_data_path = os.path.join(testdir, "data", "alaska1964_test_data.tt3")
    test_dtopo = dtopotools.DTopography(path=test_data_path)

    temp_path = tempfile.mkdtemp()
    try:
        dtopo_paths = [os.path.join(temp_path, 'alaska1964.tt1'),
                       os.path.join(temp_path, 'alaska1964.tt3')]
                       # os.path.join(temp_path, 'alaska1964.tt2'),

        for path in dtopo_paths:
            test_dtopo.write(path)
            dtopo = dtopotools.DTopography(path=path)

            assert test_dtopo.DZ.shape == dtopo.DZ.shape, \
                   "Shape of DZ not equal for topo_type = %s." % k

            assert numpy.allclose(test_dtopo.DZ, dtopo.DZ), \
                "DZ not equal for %s" % path



    except AssertionError as e:
        test_dump_path = os.path.join(os.getcwd(), "test_dtopo_io")
        shutil.mkdir(test_dump_path)
        shutil.copy(temp_path, test_dump_path)
        raise e
    finally:
        shutil.rmtree(temp_path)


def test_geometry():
    pass


def test_vs_old_dtopo():
    r"""Test new dtopotools with old version from 5.2"""

    raise nose.SkipTest("Skipping comparison with old tools.")

    import old_dtopotools

    temp_path = tempfile.mkdtemp()
    try:
        subfault_path = os.path.join(testdir, 'data', 'alaska1964.csv')
        old_dtopo_file = os.path.join(temp_path, 'old_alaska.tt1')
        old_dtopotools.builddynamicdeffile(subfault_path, subfault_path, 
                                           old_dtopo_file)
    
        input_units = {"length":"km", "width":"km", "depth":"km", "slip":"m", 
                       "mu":"dyne/cm^2"}
        fault = dtopotools.CSVFault()
        fault.read(subfault_path, input_units=input_units, 
                                  coordinate_specification="noaa sift")
        new_dtopo = fault.create_dtopography()

        X, Y, dZ = old_dtopotools.read_dtopo_old(old_dtopo_file, 
                                                 deftype='dynamic', 
                                                 only_last=False)

        assert numpy.allclose(X, new_dtopo.X), \
               "X values do not agree between old and new dtopotools."

        assert numpy.allclose(Y, new_dtopo.Y), \
               "Y values do not agree between old and new dtopotools."

        assert numpy.allclose(dZ, new_dtopo.dZ), \
               "dZ values do not agree between old and new dtopotools."

    except AssertionError as e:
        shutil.copy(temp_path, )
        raise e
    finally:
        shutil.rmtree(temp_path)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        if "plot" in sys.argv[1].lower():
            pass
        elif bool(sys.argv[1]):
            test_read_csv_make_dtopo(save=True)
    else:
        test_read_csv_make_dtopo()
        test_read_ucsb_make_dtopo()
        test_read_sift_make_dtopo()
        test_SubdividedPlaneFault_make_dtopo()
        test_dtopo_io()
        # test_vs_old_dtopo()
        # test_geometry()
