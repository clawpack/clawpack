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
import urllib
import tarfile
import time
import glob

import numpy

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


def strip_archive_extensions(path, extensions=["tar", "tgz", "bz2", "gz"]):
    r"""
    Strip off archive extensions defined in *extensions* list.

    Return stripped path calling this function recursively until all splitext
    does not provide an extension in the *extensions* list.

    """

    if os.path.splitext(path)[-1][1:] in extensions:
        return strip_archive_extensions(os.path.splitext(path)[0])
    else:
        return path


# TODO: Maybe rename this to `GeoClawRegressionTest`
class GeoClawTest(unittest.TestCase):

    r"""Base GeoClaw regression test setup


    :TODO:
     - Should we keep track of remote files?  Nothing is done with them now.
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

    def get_remote_file(self, url, force=False, verbose=False):
        r"""Fetch file located at *url* and store in object's *temp_path*.

        Will check downloaded file's suffix to see if the file needs to be
        un-archived.

        """

        if self.temp_path is None:
            raise ValueError("Temporary data directory has not been ",
                             "set yet.  Try calling super before attempting ",
                             "to get remote files.")

        file_name = os.path.basename(url)
        output_path = os.path.join(self.temp_path, file_name)
        unarchived_output_path = strip_archive_extensions(output_path)

        if not os.path.exists(unarchived_output_path) or force:
            if not os.path.exists(output_path):
                if verbose:
                    print "Downloading %s to %s..." % (url, output_path)
                urllib.urlretrieve(url, output_path)
                if verbose:
                    print "Done downloading."

            if tarfile.is_tarfile(output_path):
                if verbose:
                    print "Un-archiving %s to %s..." % (output_path, unarchived_output_path)
                with tarfile.open(output_path, mode="r:*") as tar_file:
                    tar_file.extractall(path=self.temp_path)
                if verbose:
                    print "Done un-archiving."
            # TODO: Should check here if a file is a bare compressed file (no tar)
        else:
            if verbose:
                print "Skipping %s because it already exists locally." % url

        self.remote_files.append(output_path)


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

        """

        # Write out data files
        self.load_rundata()
        self.write_rundata_objects()

        # Run code
        self.run_code()

        # Perform tests
        self.check_gauges(save=save, indices=(2, 3))


    def check_gauges(self, save=False, indices=(2, 3)):
        r"""Basic test to assert gauge equality

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

        # If we have gotten here then we do not need to copy the run results
        self.success = True


    def tearDown(self):
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

