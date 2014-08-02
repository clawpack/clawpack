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

    subfault_path = os.path.join(testdir, 'data', 'alaska1964.csv')
    units = {"length":"km", "width":"km", "depth":"km", "slip":"m", 
             "mu":"dyne/cm^2"}
    fault = dtopotools.CSVFault()
    fault.read(subfault_path, units=units, coordinate_specification="noaa sift")

    assert abs(fault.Mw() - 9.2) < 1e-4, "*** Mw is wrong: %g" % fault.Mw()

    xlower = 203.5
    xupper = 214.  # approximate - adjusted below
    ylower = 54.5
    yupper = 60.  # approximate - adjusted below

    # dtopo parameters:
    points_per_degree = 60  # 1 minute resolution
    dx = 1./points_per_degree
    mx = int((xupper - xlower)/dx + 1)
    xupper = xlower + (mx-1)*dx
    my = int((yupper - ylower)/dx + 1)
    yupper = ylower + (my-1)*dx

    x = numpy.linspace(xlower,xupper,mx)
    y = numpy.linspace(ylower,yupper,my)

    dtopo = fault.create_deformation_array(x,y,times=[1.])
    
    temp_path = tempfile.mkdtemp()
    try:
        dtopo_fname = os.path.join(temp_path, 'alaska1964.tt3')
        dtopo.write(dtopo_fname, dtopo_type=3)
        dtopo.dz_list[-1].min()
        dtopo2 = dtopotools.DTopography()
        dtopo2.read(dtopo_fname, dtopo_type=3)

        # Check that this data looks right:
        assert len(dtopo2.x) == 631,  "*** length of x is wrong"
        assert len(dtopo2.y) == 331,  "*** length of y is wrong"
        dz_max = dtopo2.dz_list[-1].max()
        assert abs(dz_max - 15.368266081250006) < 1e-5, \
                 "*** dz_max is wrong: %g" % dz_max

    except AssertionError as e:
        shutil.copy(dtopo_fname, os.path.join(os.getcwd(), "alaska1964.tt3"))
        raise e
    finally:
        shutil.rmtree(temp_path)


def test_NOAA_Sift_make_dtopo(save=False):
    pass


def test_geometry():
    pass


if __name__ == "__main__":
    if len(sys.argv) > 1:
        if "plot" in sys.argv[1].lower():
            pass
        elif bool(sys.argv[1]):
            test_read_csv_make_dtopo(save=True)
    else:
        test_read_csv_make_dtopo()
        test_NOAA_Sift_make_dtopo()
        test_geometry()