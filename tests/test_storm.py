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
    base_storm = storm.Storm(os.path.join(testdir, "data", "bal092008.dat"), 
                             file_format='atcf')

    # Create temp directory
    temp_path = tempfile.mkdtemp()

    try:
        for file_format in ['geoclaw']:
            # Write out in chosen format
            base_storm.write(os.path.join(temp_path, '%s.txt' % file_format), 
                             file_format=file_format)

            # Read in format and compare
            new_storm = storm.Storm(os.path.join(temp_path, 
                                                 '%s.txt' % file_format),
                                    file_format=file_format)

            assert(base_storm == new_storm)
            
    except AssertionError as e:
        # If the assertion failed then copy the contents of the directory
        shutil.copytree(temp_path, os.path.join(os.getcwd(),
                                                'test_storm_IO'))
        raise e

    finally:
        shutil.rmtree(temp_path)


if __name__ == '__main__':
    test_storm_IO()
