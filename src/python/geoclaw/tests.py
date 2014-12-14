r"""
Execute nosetests in all subdirectories, to run a series of quick
regression tests.

Sends output and result/errors to separate files to simplify checking
results and looking for errors.
"""

import os
import sys
import tempfile
import subprocess
import unittest
import shutil
import inspect
import time
import glob

import numpy

import clawpack.geoclaw.util

# Support for WIP decorator
from functools import wraps
from nose.plugins.attrib import attr
from nose.plugins.skip import SkipTest

# Clean library files whenever this module is used
if os.environ.has_key("CLAW"):
    CLAW = os.environ["CLAW"]
else:
    raise ValueError("Need to set CLAW environment variable.")

for lib_path in [os.path.join(CLAW,"amrclaw","src","2d"),
                 os.path.join(CLAW,"geoclaw","src","2d","shallow"),
                 os.path.join(CLAW,"geoclaw","src","2d","shallow","multilayer"),
                 os.path.join(CLAW,"geoclaw","src","2d","shallow","surge")]:
    for path in glob.glob(os.path.join(lib_path,"*.o")):
        os.remove(path)
    for path in glob.glob(os.path.join(lib_path,"*.mod")):
        os.remove(path)


# Work in progress test decorator
def fail(message):
    raise SkipTest(message)
 
def wip(f):
    @wraps(f)
    def run_test(*args, **kwargs):
        try:
            # Set to success so we don't save out the output when we know things
            # are awry
            args[0].success = True
            f(*args, **kwargs)
        except Exception as e:
            raise SkipTest("WIP test failed: " + str(e))
        fail("test passed but marked as work in progress")
    
    return attr('wip')(run_test)


# TODO: Maybe rename this to `GeoClawRegressionTest`
class GeoClawTest(unittest.TestCase):

    r"""Base GeoClaw regression test setup

    All regression tests for GeoClaw are derived from this base class.  The
    class conforms to the *unittest* modules standards including *setUp*, 
    *runTest* and *tearDown* methods which implement various parts of the 
    regression test.  Generally these methods are responsible for the following:

     - *setUp*: Creates the temprorary directory that will house test output and
       data.   Also instantiates catching of output to both *stdout* and 
       *stderr*.  Finally, this method also calls *build_executable* which will
       call the local Makefile's build process via *make .exe*.
     - *runTest*: Actually runs the test calling the following functions by
       default:
        - *load_rundata(): Creates the *rundata* objects via the local 
          *setrun.py* file.  The resulting *rundata* object is stored in the 
          class attribute *rundata*.
        - *write_rundata_objects*: Writes out all the data objects found in the
          *rundata* class attribute.
        - *run_code*: Runs the simulation based on the *runclaw.py* script 
          (what *make output* usually calls) in a subprocess to preserve 
          parallelism.  Note that all output is redirected to the class 
          attributes *stderr* and *stdout*.
        - *check_gauges*: Default test check provided by this class.  Simply 
          checks the gauges recorded from the simulation that they are equal in
          sum to the default test data.
     - *tearDown*:  Closes output redirection and tests for success via the
       class attribute of a similar name.  If the tests were not successful the
       temporary directory contents is copied to the local directory.  The 
       temporary directory is removed at this point.


    """

    def __init__(self, methodName="runTest"):

        super(GeoClawTest, self).__init__(methodName=methodName)

        self.stdout = None
        self.stderr = None
        self.success = False
        self.remote_files = []
        self.temp_path = None
        self.test_path = os.path.dirname(inspect.getfile(self.__class__))
        if len(self.test_path) == 0:
             self.test_path = "./"
        self.rundata = None

    def get_remote_file(self, url, **kwargs):
        r"""Fetch file located at *url* and store in object's *temp_path*.

        Will check downloaded file's suffix to see if the file needs to be
        un-archived.

        """

        if self.temp_path is None:
            raise ValueError("Temporary data directory has not been ",
                             "set yet.  Try calling super before attempting ",
                             "to get remote files.")

        output_path = clawpack.geoclaw.util.get_remote_file(url, 
                                                      output_dir=self.temp_path, 
                                                      **kwargs)
        self.remote_files.append(output_path)
        return output_path


    def setUp(self):
        r"""Create temp dir for data and setup log files.

        """

        self.temp_path = tempfile.mkdtemp()

        self.stdout = open(os.path.join(self.temp_path, "run_output.txt"), "w")
        self.stdout.write("Output from Test %s\n" % self.__class__.__name__)
        # TODO - Should change this to use the time module's formatting 
        # apparatus
        tm = time.localtime()
        year = str(tm[0]).zfill(4)
        month = str(tm[1]).zfill(2)
        day = str(tm[2]).zfill(2)
        hour = str(tm[3]).zfill(2)
        minute = str(tm[4]).zfill(2)
        second = str(tm[5]).zfill(2)
        date = 'Started %s/%s/%s-%s:%s.%s\n' % (year,month,day,hour,minute,second)
        self.stdout.write(date)
        self.stdout.write(("="*80 + "\n"))

        self.stderr = open(os.path.join(self.temp_path, "error_output.txt"), "w")
        self.stderr.write("Errors from Test %s\n" % self.__class__.__name__)
        self.stderr.write(date)
        self.stderr.write(("="*80 + "\n"))

        self.stdout.flush()
        self.stderr.flush()

        self.stdout.write("Paths:")
        self.stdout.write("  %s" % self.temp_path)
        self.stdout.write("  %s" % self.test_path)
        self.stdout.flush()
        self.build_executable()


    def build_executable(self, executable_name="xgeoclaw"):
        r"""Build executable by running `make .exe` in test directory.

        Moves the resulting executable to the temporary directory.


        """

        try:
            self.stdout.write("Test path and class info:\n")
            self.stdout.write("  class: %s\n" % str(self.__class__))
            self.stdout.write("  class file: %s\n" % str(inspect.getfile(self.__class__)))
            self.stdout.write("  test path: %s\n" % str(self.test_path))
            self.stdout.write("  temp path: %s\n" % str(self.temp_path))
            subprocess.check_call("cd %s ; make .exe" % self.test_path, 
                                                        stdout=self.stdout,
                                                        stderr=self.stderr,
                                                        shell=True)
        except subprocess.CalledProcessError as e:
            self.tearDown()
            raise e
        shutil.move(os.path.join(self.test_path, executable_name),  
                    self.temp_path)


    def load_rundata(self):
        r"""(Re)load setrun module and create *rundata* object


        """

        if sys.modules.has_key('setrun'):
            del(sys.modules['setrun'])
        sys.path.insert(0, self.test_path)
        import setrun
        self.rundata = setrun.setrun()
        sys.path.pop(0)


    def write_rundata_objects(self, path=None):
        r"""Write out data in the *rundata* object to *path*

        Defaults to the temporary directory path *temp_path*.

        """

        if path is None:
            path = self.temp_path

        orig_path = os.getcwd()
        os.chdir(path)
        self.rundata.write()
        os.chdir(orig_path)


    def run_code(self):
        r"""Run test code given an already compiled executable"""

        runclaw_cmd = " ".join((
                            "cd %s ;" % self.temp_path,
                            "python",
                            "$CLAW/clawutil/src/python/clawutil/runclaw.py",
                            "xgeoclaw",
                            self.temp_path,
                            "True",
                            "False",
                            self.temp_path))
        subprocess.check_call(runclaw_cmd, stdout=self.stdout, 
                                           stderr=self.stderr,
                                           shell=True)
        self.stdout.flush()
        self.stderr.flush()


    def runTest(self, save=False, indices=(2, 3)):
        r"""Basic run test functionality

        Note that this stub really only runs the code and performs no tests.

        :Input:
         - *save* (bool) - If *True* will save the output from this test to 
           the file *regresion_data.txt*.  Passed to *check_gauges*.  Default is
           *False*.
         - *indices* (tuple) - Contains indices to compare in the gague 
           comparison and passed to *check_gauges*.  Defaults to *(2, 3)*.

        """

        # Write out data files
        self.load_rundata()
        self.write_rundata_objects()

        # Run code
        self.run_code()

        # Perform tests
        self.check_gauges(save=save, indices=(2, 3))

        # If we have gotten here then we do not need to copy the run results
        self.success = True

    def check_gauges(self, save=False, indices=(2, 3)):
        r"""Basic test to assert gauge equality

        :Input:
         - *save* (bool) - If *True* will save the output from this test to 
           the file *regresion_data.txt*.  Default is *False*.
         - *indices* (tuple) - Contains indices to compare in the gague 
           comparison.  Defaults to *(2, 3)*.
        """

        # Get gauge data
        data = numpy.loadtxt(os.path.join(self.temp_path, 'fort.gauge'))
        data_sum = []
        for index in indices:
            data_sum.append(data[:, index].sum())

        # Get (and save) regression comparison data
        regression_data_file = os.path.join(self.test_path, "regression_data.txt")
        if save:
            numpy.savetxt(regression_data_file, data)
        regression_data = numpy.loadtxt(regression_data_file)
        regression_sum = []
        for index in indices:
            regression_sum.append(regression_data[:, index].sum())
        # regression_sum = regression_data

        # Compare data
        tolerance = 1e-14
        assert numpy.allclose(data_sum, regression_sum, tolerance), \
                "\n data: %s, \n expected: %s" % (data_sum, regression_sum)
        assert numpy.allclose(data, regression_data, tolerance), \
                "Full gauge match failed."



    def check_fgmax(self, save=False):
        r"""Basic test to assert fgmax equality
        Currently just records sum of fg.h and of fg.s.

        :Input:
         - *save* (bool) - If *True* will save the output from this test to 
           the file *regresion_data.txt*.  Default is *False*.
        """

        from clawpack.geoclaw import fgmax_tools

        fg = fgmax_tools.FGmaxGrid()
        fname = os.path.join(self.temp_path, 'fgmax1.txt')
        fg.read_input_data(fname)
        fg.read_output(outdir=self.temp_path)

        data_sum = numpy.array([fg.h.sum(), fg.s.sum()])

        # Get (and save) regression comparison data
        regression_data_file = os.path.join(self.test_path, \
                "regression_data_fgmax.txt")
        if save:
            numpy.savetxt(regression_data_file, data_sum)
        regression_sum = numpy.loadtxt(regression_data_file)

        # Compare data
        tolerance = 1e-14
        assert numpy.allclose(data_sum, regression_sum, tolerance), \
                "\n data: %s, \n expected: %s" % (data_sum, regression_sum)


    def tearDown(self):
        r"""Tear down test infrastructure.

        Closes *stdout* and *stderr*, removes the temporary directoy and if 
        *success* is *False*, copies the contents of the tempory directory to
        the current working directory.

        """
        self.stdout.close()
        self.stderr.close()

        if not self.success:
            output_dir = os.path.join(os.getcwd(),
                                         "%s_output" % self.__class__.__name__)

            # Remove directory if it exists already
            if os.path.exists(output_dir):
                shutil.rmtree(output_dir)

            # Copy out output files to local directory
            shutil.copytree(self.temp_path, output_dir)

        shutil.rmtree(self.temp_path)
        self.temp_path = None
