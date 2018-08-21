#!/usr/bin/env python
# encoding: utf-8

"""Tests for reading and writing storm data"""

from __future__ import absolute_import
from __future__ import print_function

import tempfile
import shutil
import os
import sys
import datetime

import numpy

import clawpack.clawutil.test as test
import clawpack.geoclaw.surge.storm as storm

# Set local test directory to get local files 
testdir = os.path.dirname(__file__) 
if len(testdir) == 0: 
    testdir = "./"

# Current tests
file_format_tests = ['atcf', 'hurdat', 'jma', 'tcvitals', 'ibtracs']

def check_geoclaw(paths, check_header=False):
    """Check that two geoclaw formatted storm files are identical

    This does not use object equivalence due to round off errors that can occur
    due to format constraints.  If the *check_header* is True the routine also
    checks that the number of lines is equivalent (is implicitly checked 
    anyway) but more importantly that the *time_offset*s are equivalent.

    :Input:
     - *check_header* (bool) Check that the headers in the file are equivalent.
       Defaults to `False`.

    :Raises:
     - *AssertionError* - If the files do not agree to the precision required.

    """

    if check_header:
        with open(paths[0], 'r') as data_file[0], open(paths[1], 'r') as data_file[1]:
            # Check for number of lines
            assert(int(data_file[0].readline()) == int(data_file[1].readline()))

            # Check for time offset
            assert(data_file[0].readline() == data_file[1].readline())

    # Check rest of data
    data = []
    for path in paths:
        data.append(numpy.loadtxt(path, skiprows=3))
    numpy.testing.assert_almost_equal(data[0], data[1])

# TODO - turn this into a test generator for each file IO format rather than
#        a single loop
def test_storm_IO(save=False):
    r"""Test reading and writing of storm formats


    Currently this only tests reading in data in all formats save for IMD and
    writing them out in the geoclaw format.  This functionality will be added
    once full writing functionality for the other formats is complete.

    :Input:
     - *save* (list or bool) whether to save the data produced by the test as
       new test data.  This can either be a single `bool` that will be applied
       to all formats or a dictionary that should have keys for each format.  If
       a format is not included in the dict than it is assumed `False`.

    """

    save_dict = {}
    if isinstance(save, bool):
        for key in file_format_tests:
            save_dict[key] = save
    elif isinstance(save, dict):
        for key in file_format_tests:
            save_dict[key] = save.get(key, default=False)
    else:
        raise ValueError("Type %s is not valid for save argument." % type(save))

    # Create temp directory
    temp_path = tempfile.mkdtemp()

    try:
        # Currently we read in the format, write it back out in the GeoClaw 
        # format and check the stored GeoClaw file for that format
        for file_format in file_format_tests:
            if file_format=='ibtracs':
                file_suffix = 'nc'
                # Check here to see if we have xarray and pandas
                try:
                    import pandas
                    import xarray
                except ImportError as e:
                    print("Skipping IBTrACS IO test, missing pandas and xarray.")
                    continue

            else:
                file_suffix = 'txt'
            input_path = os.path.join(testdir, "data", "storm", "%s.%s" % (file_format,file_suffix))
            out_path = os.path.join(temp_path, '%s_geoclaw.txt' % file_format)
            check_path = os.path.join(testdir, "data", "storm", 
                                      "%s_geoclaw.txt" % file_format)

            # Read in test data and write it back out in the GeoClaw format
            # for IBTrACS input, need storm/year info
            if file_format=='ibtracs':
                kwargs = {'storm_name':'IKE',
                         'year':2008}
            else:
                kwargs = {}
            test_storm = storm.Storm(input_path, file_format=file_format, **kwargs)
            
            # Temporary testing thing to get around missing data in formats that
            # do not provide the proper radii
            if file_format in ['hurdat', 'jma']:
                test_storm.max_wind_radius[:] = 0.0
                test_storm.storm_radius[:] = 0.0

            test_storm.write(out_path, file_format="geoclaw")

            # Save new geoclaw test files into check_path if requested
            if save:
                test_storm.write(check_path, file_format="geoclaw")

            # Check geoclaw files
            check_geoclaw([out_path, check_path])

    except Exception as e:
        # If the assertion failed then copy the contents of the directory
        test_dump_path = os.path.join(os.getcwd(), 'test_storm_IO')
        if os.path.exists(test_dump_path):
            shutil.rmtree(test_dump_path)
        shutil.copytree(temp_path, os.path.join(os.getcwd(),
                                                'test_storm_IO'))
        print("Format test %s -> geoclaw  failed." % file_format)
        raise e

    finally:
        shutil.rmtree(temp_path)


if __name__ == '__main__':
    # Currently does not support only saving one of the format's data
    save = False
    if len(sys.argv) > 1:
        save = bool(sys.argv[1])
    test_storm_IO(save)
