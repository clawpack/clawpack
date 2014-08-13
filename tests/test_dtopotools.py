#!/usr/bin/env python

import os
import sys
import shutil
import tempfile

import numpy

import clawpack.geoclaw.topotools as topotools
import clawpack.geoclaw.dtopotools as dtopotools

# Set local test directory to get local files
testdir = os.path.dirname(__file__)
if len(testdir) == 0:
     testdir = "./"

def test_read_csv_make_dtopo(save=False):
    r"""Test reading and making of a CSV subfault speficied dtopo"""
    
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
    dx = 1./points_per_degree
    mx = int((xupper - xlower)/dx + 1)
    xupper = xlower + (mx-1)*dx
    my = int((yupper - ylower)/dx + 1)
    yupper = ylower + (my-1)*dx

    x = numpy.linspace(xlower,xupper,mx)
    y = numpy.linspace(ylower,yupper,my)

    dtopo = fault.create_dtopography(x,y,times=[1.])
    print "max dz = ",dtopo.dz_list[-1].max()
    
    temp_path = tempfile.mkdtemp()
    try:
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


    except AssertionError as e:
        shutil.copy(test_data_path, os.path.split(test_data_path)[1])
        raise e
    finally:
        shutil.rmtree(temp_path)

def test_read_ucsb_make_dtopo(save=False):

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
    dx = 1./points_per_degree
    mx = int((xupper - xlower)/dx + 1)
    xupper = xlower + (mx-1)*dx
    my = int((yupper - ylower)/dx + 1)
    yupper = ylower + (my-1)*dx

    x = numpy.linspace(xlower,xupper,mx)
    y = numpy.linspace(ylower,yupper,my)

    tmax = 0.
    for s in fault.subfaults:
        tmax = max(tmax, s.rupture_time + s.rise_time + s.rise_time_ending)

    fault.rupture_type = 'dynamic'
    times = numpy.linspace(0,tmax,10)
    dtopo = fault.create_dtopography(x,y,times)
    print "max dz = ",dtopo.dz_list[-1].max()
    
    temp_path = tempfile.mkdtemp()
    try:
        # Load (and save) test data and make the comparison
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


    except AssertionError as e:
        shutil.copy(test_data_path, os.path.split(test_data_path)[1])
        raise e
    finally:
        shutil.rmtree(temp_path)


def test_read_sift_make_dtopo(save=False):
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
    print "max dz = ",dtopo.dz_list[-1].max()
    
    temp_path = tempfile.mkdtemp()
    try:
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

    except AssertionError as e:
        shutil.copy(test_data_path, os.path.split(test_data_path)[1])
        raise e
    finally:
        shutil.rmtree(temp_path)


def test_SubdividedPlaneFault_make_dtopo(save=False):

    # get a unit source fault plane as starting point:
    sift_slip = {'acsza1':1.}
    fault = dtopotools.SiftFault(sift_slip)
    fault_plane = fault.subfaults[0]
    Mo = fault_plane.Mo()
    print "original Mo = ",Mo

    fault2 = dtopotools.SubdividedPlaneFault(fault_plane, nstrike=5, ndip=3)
    print "new Mo = ",fault2.Mo()
    #fault2.plot_subfaults(slip_color=True)

    assert abs(fault2.Mw() - 7.50068666) < 1e-4, "*** Mw is wrong: %g" % fault.Mw()

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
    print "max dz = ",dtopo.dz_list[-1].max()
    
    temp_path = tempfile.mkdtemp()
    try:
        # Load (and save) test data and make the comparison
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

    except AssertionError as e:
        shutil.copy(test_data_path, os.path.split(test_data_path)[1])
        raise e
    finally:
        shutil.rmtree(temp_path)



def test_geometry():
    pass


def test_read_write_dtopo():
    r"""Test new dtopo read functionality vs. old version"""
    
    import old_dtopotools

    # Write out supported dtopo files that are supported
    old_dtopotools

    old_dtopotools.read_dtopo()



def test_vs_old_dtopo():
    r"""Test new dtopotools with old version from 5.2"""

    import old_dtopotools

    subfault_path = os.path.join(testdir, 'data', 'alaska1964.csv')
    old_subfaults = old_dtopotools.read_subfault_model_csv(subfault_path)
    dtopo_params = {}
    old_dtopo = old_dtopotools.make_dtopo_from_subfaults(old_subfaults, 
                                                         dtopo_params)
    


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
        
        test_geometry()
