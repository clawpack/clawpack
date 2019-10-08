#!/usr/bin/env bash

set -e

cd $CLAW

cat /proc/cpuinfo
gfortran -v
echo "DEBUGGING --------------------"
echo "CLAW="$CLAW
echo "PYTHONPATH="$PYTHONPATH
yolk -l
echo "------------------------------"

echo
echo "Beginning tests"
echo "------------------------------"

echo "testing pyclaw"
cd $CLAW/pyclaw/examples
nosetests || for failed_test_path in *_output ; do cat $failed_test_path/run_output.txt ; cat $failed_test_path/error_output.txt ; done

echo "testing classic"
cd $CLAW/classic
nosetests || for failed_test_path in *_output ; do cat $failed_test_path/run_output.txt ; cat $failed_test_path/error_output.txt ; done


echo "testing amrclaw"
cd $CLAW/amrclaw
nosetests || for failed_test_path in *_output ; do cat $failed_test_path/run_output.txt ; cat $failed_test_path/error_output.txt ; done


echo "testing geoclaw"
cd $CLAW/geoclaw
nosetests || for failed_test_path in *_output ; do cat $failed_test_path/run_output.txt ; cat $failed_test_path/error_output.txt ; done
