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

    assert abs(fault.Mw() - 9.2) < 1e-4, "*** Mw is wrong: %g" % fault.Mw()

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

    test_data_path = os.path.join(testdir, "data", "alaska1964_test_data.tt3")
    if save:
        dtopo.write(test_data_path, dtopo_type=3)
    compare_data = dtopotools.DTopography(path=test_data_path)
    compare_data.read(path=test_data_path, dtopo_type=3)

    assert len(dtopo.dz_list) == len(compare_data.dz_list), \
        "len(dtopo.dz_list) is %s, should be %s" \
        % (len(dtopo.dz_list), len(compare_data.dz_list))

    for k in range(len(dtopo.dz_list)):
        assert numpy.allclose(compare_data.dz_list[k], dtopo.dz_list[k])


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

    test_data_path = os.path.join(testdir, "data", "tohoku_test_data.tt3")
    if save:
        dtopo.write(test_data_path, dtopo_type=3)
    compare_data = dtopotools.DTopography(path=test_data_path)
    compare_data.read(path=test_data_path, dtopo_type=3)

    assert len(dtopo.dz_list) == len(compare_data.dz_list), \
        "len(dtopo.dz_list) is %s, should be %s" \
        % (len(dtopo.dz_list), len(compare_data.dz_list))

    for k in range(len(dtopo.dz_list)):
        assert numpy.allclose(compare_data.dz_list[k], dtopo.dz_list[k])


def test_read_sift_make_dtopo(save=False):
    r"""Test reading and making of a SIFT subfault speficied dtopo"""

    sift_slip = {'acsza1':2, 'acszb1':3}
    fault = dtopotools.SiftFault(sift_slip)

    assert abs(fault.Mw() - 7.966666666) < 1e-4, "*** Mw is wrong: %g" % fault.Mw()

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

    assert len(dtopo.dz_list) == len(compare_data.dz_list), \
        "len(dtopo.dz_list) is %s, should be %s" \
        % (len(dtopo.dz_list), len(compare_data.dz_list))

    for k in range(len(dtopo.dz_list)):
        assert numpy.allclose(compare_data.dz_list[k], dtopo.dz_list[k])

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

    assert abs(fault2.Mw() - 7.50068666) < 1e-4, \
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

    test_data_path = os.path.join(testdir, "data", 
            "SubdividedFaultPlane_test_data.tt3")
    if save:
        dtopo.write(test_data_path, dtopo_type=3)
    compare_data = dtopotools.DTopography(path=test_data_path)
    compare_data.read(path=test_data_path, dtopo_type=3)

    assert len(dtopo.dz_list) == len(compare_data.dz_list), \
        "len(dtopo.dz_list) is %s, should be %s" \
        % (len(dtopo.dz_list), len(compare_data.dz_list))

    for k in range(len(dtopo.dz_list)):
        assert numpy.allclose(compare_data.dz_list[k], dtopo.dz_list[k])


def test_dtopo_io():
    r"""Test IO of dtopography class"""

    test_data_path = os.path.join(testdir, "data", "alaska1964_test_data.tt3")
    test_dtopo = dtopotools.DTopography(path=test_data_path)

    temp_path = tempfile.mkdtemp()
    try:
        dtopo_paths = [os.path.join(temp_path, 'alaska1964.tt1'),
                       os.path.join(temp_path, 'alaska1964.tt3')]
                       # os.path.join(temp_path, 'alaska1964.tt2'),

        dtopo = []
        for (k, path) in enumerate(dtopo_paths):
            test_dtopo.write(path)
            dtopo = dtopotools.DTopography(path=path)

            assert len(test_dtopo.dz_list) == len(dtopo.dz_list), \
                   "Number of dz fields not equal for topo_type = %s." % k

            for (n, dz) in enumerate(test_dtopo.dz_list):
                assert numpy.allclose(dz, dtopo.dz_list[n]), \
                       "dz field not equal for topo_type = %s " + \
                       "at time = %s." % (k, dtopo.times[n])


    except AssertionError as e:
        test_dump_path = os.path.join(os.getcwd(), "test_dtopo_io")
        shutil.mkdir(test_dump_path)
        shutil.copy(temp_path, test_dump_path)
        raise e
    finally:
        shutil.rmtree(temp_path)


def test_geometry():
    r"""Test subfault geometry calculation."""

    import old_dtopotools

    subfault_path = os.path.join(testdir, 'data', 'alaska1964.csv')
    input_units = {"length":"km", "width":"km", "depth":"km", "slip":"m", 
         "mu":"dyne/cm^2"}
    specifications = ['top center', 'centroid', 'bottom center', 'noaa sift']

    for specification in specifications:
        fault = dtopotools.CSVFault()
        fault.read(subfault_path, input_units=input_units, 
                                  coordinate_specification=specification)

        # Subfault 10 is chosen at random, maybe do all?
        subfault = fault.subfaults[10]
        geometry = old_dtopotools.set_geometry(subfault)

        coord_tests = {"top center":{'test':[geometry['x_top'], 
                                             geometry['y_top'], 
                                             geometry['depth_top']], 
                                 'computed':subfault.fault_plane_centers[0]},
                       "centroid":{'test':[geometry['x_centroid'],
                                           geometry['y_centroid']],
                               'computed':subfault.fault_plane_centers[1][:2]},
                       "bottom center":{"test":[geometry['x_bottom'],
                                                geometry['y_bottom'],
                                                geometry['depth_bottom']],
                                    "computed":subfault.fault_plane_centers[2]},
                       "Corner A":{"test":[geometry["x_corners"][2],
                                           geometry["y_corners"][2]],
                               "computed":subfault.fault_plane_corners[0][:2]},
                       "Corner B":{"test":[geometry["x_corners"][3],
                                           geometry["y_corners"][3]],
                               "computed":subfault.fault_plane_corners[1][:2]},
                       "Corner C":{"test":[geometry["x_corners"][0],
                                           geometry["y_corners"][0]],
                               "computed":subfault.fault_plane_corners[2][:2]},
                       "Corner D":{"test":[geometry["x_corners"][1],
                                           geometry["y_corners"][1]],
                               "computed":subfault.fault_plane_corners[3][:2]}

                    }

        for (values, coord_test) in coord_tests.iteritems():
            assert numpy.allclose(coord_test['test'], coord_test['computed']), \
                   "Specification = %s, coords= %s:\n%s !=\n%s" % (
                                                         specification, 
                                                         values, 
                                                         coord_test['test'], 
                                                         coord_test['computed'])


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
        try:
            test_read_csv_make_dtopo()
            test_read_ucsb_make_dtopo()
            test_read_sift_make_dtopo()
            test_SubdividedPlaneFault_make_dtopo()
            test_dtopo_io()
            test_geometry()
            test_vs_old_dtopo()
        except nose.SkipTest as e:
            print e.message
