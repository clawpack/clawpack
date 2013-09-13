
from clawpack.visclaw import data
import os
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

    g1 = plotdata.getgauge(1)
    tsum = g1.t.sum()
    qsum = g1.q[0,:].sum()

    if save_new_regression_data:
        np.savetxt('regression_data_test2.txt',np.array([tsum,qsum]))

    # Read in archived data for comparison:
    regression_data = np.loadtxt('regression_data_test2.txt')
    tsum_expected = regression_data[0]
    qsum_expected = regression_data[1]

    tol = 1e-14
    assert np.allclose(tsum,tsum_expected,tol), \
        "gauge 1: tsum = %s, expected: %s"  % (tsum, tsum_expected)
    assert np.allclose(qsum,qsum_expected,tol), \
        "gauge 1: qsum = %s, expected: %s"  % (qsum, qsum_expected)
    print "Gauge 1 OK"
    
    
if __name__=="__main__":
    save_new_regression_data = False  # Set to True to archive new results
    test1()
    test2(save_new_regression_data)
    if save_new_regression_data:
        print "*** Created new regression_data files ***"
    


