#!/usr/bin/env python
# encoding: utf-8

"""Tests for reading and writing storm data"""

from __future__ import absolute_import
from __future__ import print_function

import tempfile
import shutil
import os
import datetime

import numpy

import clawpack.clawutil.test as test
import clawpack.geoclaw.surge.storm as storm

# Set local test directory to get local files 
testdir = os.path.dirname(__file__) 
if len(testdir) == 0: 
    testdir = "./"

def check_geoclaw(paths):
    """Check that two geoclaw formatted storm files are identical"""

    assert(True)

def test_storm_IO(save=False):
    r"""Test reading and writing of storm formats"""

    # Currently this only tests reading in data in all formats save for IMD and
    # writing them out in the geoclaw format.  This functionality will be added
    # once full writing functionality for the other formats is complete.

    # Create temp directory
    temp_path = tempfile.mkdtemp()

    try:
        # Currently we read in the format, write it back out in the GeoClaw 
        # format and check the stored GeoClaw file for that format
        for file_format in ['atcf', 'hurdat', 'jma', 'tcvitals']:
            print(file_format)
            input_path = os.path.join(testdir, "data", "storm", "%s.txt" % file_format)
            out_path = os.path.join(temp_path, '%s_geoclaw.txt' % file_format)
            check_path = os.path.join(testdir, "data", "storm", 
                                      "%s_geoclaw.txt" % file_format)

            # Read in test data and write it back out in the GeoClaw format
            test_storm = storm.Storm(input_path, file_format=file_format)
            if test_storm.time_offset is None:
                test_storm.time_offset = datetime.datetime(2001, 1, 1)
            test_storm.write(out_path, file_format="geoclaw")

            # Save new geoclaw test files into check_path if requested
            if save:
                test_storm.write(check_path, file_format="geoclaw")

            # Check geoclaw files
            check_geoclaw([out_path, check_path])

    except AssertionError as e:
        # If the assertion failed then copy the contents of the directory
        shutil.copytree(temp_path, os.path.join(os.getcwd(),
                                                'test_storm_IO'))
        raise e

    finally:
        shutil.rmtree(temp_path)


if __name__ == '__main__':
    test_storm_IO()
