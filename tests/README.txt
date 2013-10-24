This directory is for regression tests.

To run all tests:
    python run_tests.py
or: 
    make tests
which does the same thing.  This runs nosetests but also saves all output to
a file in case you want to check it.

To just run nosetests in all subdirectories:
    nosetests */

To clean up afterwards, removing all executables, output, and test results:
    make clobber

Each test does a short run with regions and gauges set to exercise the code,
The test passes if the code runs and if the sum of t values and of q values
agree with archived results for each gauge.

Developers: To create new archived results for a test case, go into the
directory and type:
    python regression_tests.py True
Then 'git add' and issue a pull request if you believe the new results are
more correct.

