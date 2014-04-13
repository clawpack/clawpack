#!/usr/bin/env python

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

class GeoClawTest(unittest.TestCase):

    r"""Base GeoClaw regression test setup


    :TODO:
    """

    def __init__(self):

        super(GeoClawTest, self).__init__()

        self.temp_path = None
        self.stdout = None
        self.stderr = None

        self.remote_files = []
        self.success = False

        self.test_path = os.path.dirname(inspect.getfile(self.__class__))

    def get_remote_file(self, url, force=False, verbose=False):
        r"""Fetch file located at *url* and store in object's *temp_path*

        """

        if self.temp_path is None:
            raise ValueError("Temporary data directory has not been ",
                             "set yet.  Try calling super before attempting ",
                             "to get remote files.")

        file_name = os.path.basename(url)
        output_path = os.path.join(self.temp_path, file_name)
        if not os.path.exists(output_path) or force:
            if verbose:
                print "Downloading %s to %s..." % (url, output_path)
            urllib.urlretrieve(url, output_path)
            if verbose:
                print "Finished downloading."
        else:
            if verbose:
                print "Skipping %s, file already exists." % file_name

        if tarfile.is_tarfile(output_path):
            with tarfile.open(output_path, mode="r:*") as tar_file:
                tar_file.extractall(path=self.temp_file)

        self.remote_files.append(output_path)


    def setUp(self):
        r"""Create temp dir for data and download files from web

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

        self.build_executable()


    def runTest(self):

        # Write out data files
        orig_path = os.getcwd()
        os.chdir(self.test_path)
        sys.path.append("./")
        import setrun
        rundata = setrun.setrun()
        os.chdir(self.temp_path)
        rundata.write()
        os.chdir(orig_path)

        # Run code
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
        subprocess.check_call("cd %s ; make .exe" % self.test_path, 
                                                        stdout=self.stdout, 
                                                        stderr=self.stderr, 
                                                        shell=True)
        shutil.move(os.path.join(self.test_path, executable_name),  
                    self.temp_path)

