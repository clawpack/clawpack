#!/usr/bin/env python
# encoding: utf-8

from __future__ import absolute_import
from __future__ import print_function
import six
from six.moves import range

import tempfile
import shutil

import numpy

import clawpack.clawutil.test as test
import clawpack.geoclaw.surge.storm as storm

test_data_path = os.path.join(testdir, "data", "geoclaw.storm")


@test.wip
def test_storm_IO():
    r"""Test reading and writing of storm formats"""

    # Create temp directory
    temp_path = tempfile.mkdtemp()
    try:
        # Load original data
        base_storm = storm.Storm(path=test_data_path, file_format="geoclaw")

        # Run through formats, write out last in new format, read it back in,
        # and then compare to original
        for file_format in enumerate(storm._supported_formats[1:]):
            new_storm_path = os.path.join(temp_path, "test.storm")
            # Write
            base_storm.write(path=new_storm_path, file_format=file_format)
            # Read
            new_storm = storm.Storm(path=new_storm_path,
                                    file_format=file_format)
            # Compare
            assert(base_storm == new_storm,
                   "File format %s failed to be read in." % file_format)

        # for topo_type in range(1, 4):
        #     path = os.path.join(temp_path, 'bowl.tt%s' % topo_type)
        #     topo.write(path, topo_type=topo_type,Z_format="%22.15e")

        #     topo_in = topotools.Topography(path)
        #     assert numpy.allclose(topo.Z, topo_in.Z), \
        #            "Differnece in written and read topography found."

    except AssertionError as e:
        # If the assertion failed then copy the contents of the directory
        shutil.copytree(temp_path, os.path.join(os.getcwd(), "test_storm"))
        raise e

    finally:
        shutil.rmtree(temp_path)


@test.wip
def test_storm_models():
    r"""Test storm model fields"""
    raise NotImplementedError("Storm model tests not implemented.")


if __name__ == '__main__':
    test_storm_IO()
    # test_storm_models()
