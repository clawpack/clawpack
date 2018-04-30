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

def test_storm_IO():
    r"""Test reading and writing of storm formats"""

    # Location of test data - Uses ATCF Ike data as basic test
    current_storm = storm.Storm(os.path.join(testdir, "data", "bal092008.dat"), 
                                file_format='atcf')

    # Create temp directory
    temp_path = tempfile.mkdtemp()

    # We only test a subset right now as some are not implemented
    test_list = ['atcf', 'geoclaw', 'hurdat', 'tcvitals', 'atcf']
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


def test_storm_models():
    r"""Test storm model fields"""

    pass

    # x = numpy.linspace(0.0, 1e5, 10)
    # t = 1.0
    # test_storm = storm.Storm(path=os.path.join(testdir, "data",
    #                                            "geoclaw.storm"),
    #                          file_format="geoclaw")
    # for model_type in storm._supported_models:
    #     # Compute model fields
    #     model_data = storm.construct_fields(test_storm, x, t, model=model_type)

    #     # Load comparison data
    #     test_data_path = os.path.join(testdir, "data", "%s.txt" % model_type)
    #     test_data = numpy.loadtxt(test_data_path)

    #     # Assert
    #     numpy.testing.assert_allclose(model_data, test_data,
    #                                   "Storm model %s failed." % model_type)


if __name__ == '__main__':
    test_storm_IO()
