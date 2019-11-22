#!/usr/bin/env bash

set -e

print_errors () {
  for failed_test_path in *_output
  do
    cat $failed_test_path/run_output.txt
    cat $failed_test_path/error_output.txt
  done
  exit 1
}

cd $CLAW

cat /proc/cpuinfo
gfortran -v
echo "DEBUGGING --------------------"
echo "CLAW="$CLAW
echo "PYTHONPATH="$PYTHONPATH
echo "------------------------------"

echo
echo "Beginning tests"
echo "------------------------------"

echo "testing pyclaw"
cd $CLAW/pyclaw/examples
nosetests -sv || print_errors

echo "testing classic"
cd $CLAW/classic
nosetests -sv || print_errors

echo "testing amrclaw"
cd $CLAW/amrclaw
nosetests -sv || print_errors

echo "testing geoclaw"
cd $CLAW/geoclaw
nosetests -sv || print_errors
