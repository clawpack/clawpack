"""
Run all examples and then compare to gallery versions, if available.
"""

from __future__ import absolute_import
from __future__ import print_function
from clawpack.clawutil import regression_tests, make_all
import os

env = os.environ
env['GIT_STATUS'] = 'True'
env['FFLAGS'] = '-O2 -fopenmp'
env['OMP_NUM_THREADS'] = '3'

make_all.make_all(make_clean_first=True, env=env)


print("\n-----------------------------------------------------------\n")

all_ok = regression_tests.test_subdirs()
if all_ok:
    print("===> All tests pass")
else:
    print("===> Some test(s) failed")
