
from __future__ import absolute_import
import os
import sys
import tempfile
import shutil

from clawpack.geoclaw import etopotools

def test_fetch_etopo():
    """
    Test fetching etopo1 data from the NCEI (NGDC) website.
    """
    
    xlimits = (-120,-60)
    ylimits = (-60,0)
    resolution = 10./60.   # in degrees

    try:
        temp_path = tempfile.mkdtemp()

        topo = etopotools.etopo1_download(xlimits,ylimits, dx=resolution, \
                output_dir=temp_path, return_topo=True)

        assert topo.Z.shape==(361,361), "*** topo file has wrong shape"
        assert topo.Z[0,0]==-4516, "*** topo.Z[0,0]=%g has unexpected value"\
                                % topo.Z[0,0]

    except AssertionError as e:
        # If the assertion failed then copy the contents of the directory
        shutil.copytree(temp_path, os.path.join(os.getcwd(), 
                                                "test_fetch_etopo"))
        raise e
    finally:
        shutil.rmtree(temp_path)
