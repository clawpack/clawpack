#!/usr/bin/env python

from __future__ import print_function
import os
import numpy
import nose
import os
from clawpack.geoclaw import topotools


# Set local test directory to get local files
testdir = os.path.dirname(__file__)
if len(testdir) == 0:
     testdir = "./"

extent = [-125,-124, 48, 48.5]

def test_etopo1_topo(make_plot=False, save=False):
    
    try:
        import netCDF4
    except:
        raise nose.SkipTest("netCDF4 not installed, skipping test")
        
    topo1 = topotools.read_netcdf('etopo1', extent=extent, verbose=True)

    topo10 = topotools.read_netcdf('etopo1', extent=extent, 
                                   coarsen=10, verbose=True)

    testdata_path = os.path.join(os.path.dirname(__file__), 'data', 'etopo1_10min.asc')
    if save:
        topo10.write(testdata_path, topo_type=3, Z_format='%.0f')
        print('Created %s' % testdata_path)

    topo10input = topotools.Topography()
    topo10input.read(testdata_path, topo_type=3)
    
    assert numpy.allclose(topo10.Z, topo10input.Z), \
           "topo10.Z does not agree with archived data"
    
    if make_plot:
        import matplotlib.pyplot as plt
        plt.figure(figsize=(12,5))
        ax1 = plt.subplot(1,2,1)
        topo1.plot(axes=ax1)
        plt.title('1 minute etopo1 data')
        ax10 = plt.subplot(1,2,2)
        topo10.plot(axes=ax10)
        plt.title('10 minute etopo1 data')
        pname = 'etopo1_test_plot.png'
        plt.savefig(pname)
        print('Created %s' % pname)
    
def test_etopo1_xarray():

    try:
        import xarray
    except:
        raise nose.SkipTest("xarray not installed, skipping test")
        
    topo10,topo10_xarray = topotools.read_netcdf('etopo1', extent=extent, 
                                                 return_xarray=True,
                                                 coarsen=10, verbose=True)

    testdata_path = os.path.join(testdir, 'data', 'etopo1_10min.asc')
    topo10input = topotools.Topography()
    topo10input.read(testdata_path, topo_type=3)
    
    assert numpy.allclose(topo10_xarray['z'], topo10input.Z), \
           "topo10_xarray['z'] does not agree with archived data"
    

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        if "plot" in sys.argv[1].lower():
            test_etopo1_topo(make_plot=True)
        elif bool(sys.argv[1]):
            test_etopo1_topo(save=True)
    else:
        # Run tests
        test_etopo1_topo()
        test_etopo1_xarray()
        print("All tests passed.")

