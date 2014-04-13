#!/usr/bin/env python

import os
import glob

# Test sub-packages
import bowl_slosh

# Determine directory:
try:
    CLAW = os.environ['CLAW']
except:
    raise ValueError("Need to set CLAW environment variable.")

# Clean library files
for lib_path in [os.path.join(CLAW,"amrclaw","src","2d"),
                 os.path.join(CLAW,"geoclaw","src","2d","shallow"),
                 os.path.join(CLAW,"geoclaw","src","2d","shallow","multilayer"),
                 os.path.join(CLAW,"geoclaw","src","2d","shallow","surge")]:
    for path in glob.glob(os.path.join(lib_path,"*.o")):
        os.remove(path)
    for path in glob.glob(os.path.join(lib_path,"*.mod")):
        os.remove(path)
