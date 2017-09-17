#!/usr/bin/env python

"""Simple implementation of a file fetcher"""

from __future__ import absolute_import
import sys
import os
import clawpack.clawutil.data

if __name__ == "__main__":
    # Default URLs
    base_url = "http://www.columbia.edu/~ktm2132/bathy/"

    # Override base_url
    if len(sys.argv) > 1:
        base_url = sys.argv[1]

    urls = [os.path.join(base_url, 'gulf_caribbean.tt3.tar.bz2')]

    for url in urls:
        clawpack.clawutil.data.get_remote_file(url)