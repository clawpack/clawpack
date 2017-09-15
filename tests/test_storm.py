#!/usr/bin/env python
# encoding: utf-8

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

#@test.wip
def test_storm_IO():
    r"""Test reading and writing of storm formats"""

    # Location of test data
    test_data_path = os.path.join(testdir, "data_storm", "geoclaw.storm")

    # Create temp directory
    temp_path = tempfile.mkdtemp()
    try:
        # Load original data
        base_storm = storm.Storm(path=test_data_path, file_format="geoclaw")
        # Run through formats, write out last in new format, read it back in,
        # and then compare to original
        # for file_format in enumerate(storm._supported_formats[1:]):
        for file_format in storm._supported_formats[1:]:
            new_storm_path = os.path.join(temp_path, "test.storm")
            # Write
            base_storm.write(path=new_storm_path, file_format=file_format)
            # Read
            new_storm = storm.Storm(path=new_storm_path,
                                    file_format=file_format)
            print('base_storm', base_storm) 
            print('new_storm', new_storm) 
            # Compare - TODO: may need to do a set of assert_allclose here
            assert base_storm == new_storm, \
                   "File format %s failed to be read in." % file_format

    except AssertionError as e:
        # If the assertion failed then copy the contents of the directory
        print('temp_path', temp_path) 
        raise e

    finally:
        shutil.rmtree(temp_path)


@test.wip
def test_storm_models():
    r"""Test storm model fields"""

    x = numpy.linspace(0.0, 1e5, 10)
    t = 1.0
    test_storm = storm.Storm(path=os.path.join(testdir, "data",
                                               "geoclaw.storm"),
                             file_format="geoclaw")
    for model_type in storm._supported_models:
        # Compute model fields
        model_data = storm.construct_fields(test_storm, x, t, model=model_type)

        # Load comparison data
        test_data_path = os.path.join(testdir, "data", "%s.txt" % model_type)
        test_data = numpy.loadtxt(test_data_path)

        # Assert
        numpy.testing.assert_allclose(model_data, test_data,
                                      "Storm model %s failed." % model_type)


if __name__ == '__main__':
    test_storm_IO()
