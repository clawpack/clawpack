#!/usr/bin/env python3
import re
import os
import sys


def get_name(settings=None):
    return "mpi4py"


def get_version(settings=None):
    topdir = os.path.dirname(os.path.abspath(__file__))
    source = os.path.join(topdir, "clawpack", "__init__.py")
    with open(source, encoding="utf-8") as f:
        m = re.search(r"__version__\s*=\s*'(.*)'", f.read())
    version = m.groups()[0]
    local_version = os.environ.get("CLAWPACK_LOCAL_VERSION")
    if local_version:
        version = "{version}+{local_version}".format(**vars())
    return version



def dynamic_metadata(field, settings=None):
    getter = globals().get("get_" + field)
    if getter:
        return getter(settings)
    return globals()[field.replace(".", "_")]


if __name__ == "__main__":
    print(dynamic_metadata(sys.argv[1]))
