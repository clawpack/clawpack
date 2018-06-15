#!/usr/bin/env python
# encoding: utf-8

"""Tests for reading and writing storm data"""

from __future__ import absolute_import
from __future__ import print_function

import tempfile
import shutil
import os 

import numpy

import clawpack.clawutil.test as test
import clawpack.geoclaw.surge.storm as storm

# Set local test directory to get local files 
testdir = os.path.dirname(__file__) 
if len(testdir) == 0: 
    testdir = "./"

def test_storm_IO():
    r"""Test reading and writing of storm formats"""

    # Currently this only tests reading in data in all formats save for IMD and
    # writing them out in the geoclaw format.  This functionality will be added
    # once full writing functionality for the other formats is complete.

    # Location of test data - Uses ATCF Ike data as basic test
    current_storm = storm.Storm(os.path.join(testdir, "data", "bal092008.dat"), 
                                file_format='atcf')

    # Create temp directory
    temp_path = tempfile.mkdtemp()

    # We only test a subset right now as some are not implemented
    test_list = ['atcf', 'geoclaw', 'atcf']
    try:
        for format_name in test_list:
            new_path = os.path.join(temp_path, "%s.txt" % format_name)
            current_storm.write(new_path, file_format=format_name)
            new_storm = storm.Storm(new_path, file_format=format_name)

            asseert(current_storm == new_storm)

            current_storm = new_storm.copy()

    except AssertionError as e:
        # If the assertion failed then copy the contents of the directory
        shutil.copytree(temp_path, os.path.join(os.getcwd(),
                                                'test_storm_IO'))
        raise e

    finally:
        shutil.rmtree(temp_path)


if __name__ == '__main__':
    test_storm_IO()
