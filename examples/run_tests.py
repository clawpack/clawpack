"""
Run all examples and then compare to gallery versions, if available.
"""

from clawpack.clawutil import regression_tests, make_all

make_all.make_all(make_clean_first=True)

print "\n-----------------------------------------------------------\n"

all_ok = regression_tests.test_subdirs()
if all_ok:
    print "===> All tests pass"
else:
    print "===> Some test(s) failed"
