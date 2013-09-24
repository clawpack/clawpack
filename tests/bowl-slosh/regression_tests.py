
"""
Regression tests.  Execute via:
    python regression_tests.py
to test, or
    python regression_tests.py True
to create new regression data for archiving.
"""

from clawpack.visclaw import data
import os,sys
import numpy as np
import subprocess

# Create plotdata object for reading in gauges in tests below:
plotdata = data.ClawPlotData()
plotdata.outdir = '_output'

def test1():
    """
    Compile and run the code
    """
    job = subprocess.Popen(['make', 'clean'])
    return_code = job.wait()
    assert return_code == 0, "Problem with 'make clean'"
    job = subprocess.Popen(['make', 'topo'])
    return_code = job.wait()
    assert return_code == 0, "Problem with 'make topo'"
    job = subprocess.Popen(['make', '.output'])
    return_code = job.wait()
    assert return_code == 0, "Problem with 'make .output'"


def test2(save_new_regression_data=False):
    """
    Check Gauge 1, only test sum of t and sum of q values.
    """
    
    # unique for this test:
    gaugeno = 1
    fname_data = 'regression_data_test2.txt'

    g = plotdata.getgauge(gaugeno)
    tsum = g.t.sum()
    qsum = g.q[0,:].sum()

    new_data = np.array([tsum,qsum])
    
    if save_new_regression_data:
        np.savetxt(fname_data, new_data)
        print "*** Created new regression_data file ", fname_data

    # Read in archived data for comparison:
    regression_data = np.loadtxt(fname_data)

    tol = 1e-14
    assert np.allclose(new_data,regression_data,tol), \
        "\n  new_data: %s, \n  expected: %s"  % (new_data, regression_data)
    print "Gauge %i ok" % gaugeno


if __name__=="__main__":
    save_new_regression_data = (len(sys.argv) > 1) and (sys.argv[1]=='True')
    test1()
    test2(save_new_regression_data)
    


