#!/usr/bin/env python

"""Simple implementation of a file fetcher"""

import sys
import os
import clawpack.geoclaw.util as util

if __name__ == "__main__":
    # Default URLs
    base_url = "http://www.columbia.edu/~ktm2132/bathy/"

    # Override base_url
    if len(sys.argv) > 1:
        base_url = sys.argv[1]

    urls = [os.path.join(base_url, 'gulf_caribbean.tt3.tar.bz2')]

    for url in urls:
        util.get_remote_file(url)