#!/usr/bin/env python

"""Simple implementation of a file fetcher"""

import sys
import os
import urllib
import tarfile
import subprocess

def get_bathy(url, destination=os.getcwd(), force=False):
    r"""Get bathymetry file located at `url`

    Will check downloaded file's suffix to see if the file needs to be extracted
    """

    file_name = os.path.basename(url)
    output_path = os.path.join(destination, file_name)
    if not os.path.exists(output_path) or force:
        print "Downloading %s to %s..." % (url, output_path)
        urllib.urlretrieve(url, output_path)
        print "Finished downloading."
    else:
        print "Skipping %s, file already exists." % file_name

    if tarfile.is_tarfile(output_path):
        with tarfile.open(output_path, mode="r:*") as tar_file:
            tar_file.extractall(path=destination)


if __name__ == "__main__":
    # Default URLs
    base_url = "http://users.ices.utexas.edu/~kyle/bathy/"

    # Override base_url
    if len(sys.argv) > 1:
        base_url = sys.argv[1]

    urls = [os.path.join(base_url, 'gulf_caribbean.tt3.tar.bz2')]

    for url in urls:
        get_bathy(url)