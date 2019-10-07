#!/usr/bin/env bash

set -e

cd /clawpack

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
nosetests

echo "testing classic"
cd $CLAW/classic
nosetests

echo "testing amrclaw"
cd $CLAW/amrclaw
nosetests

echo "testing geoclaw"
cd $CLAW/geoclaw
nosetests

cd $CLAW
