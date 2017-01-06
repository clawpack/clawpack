#!/bin/sh

# Execute this script from the top level clawpack directory, 
# assumed to be a git clone of git://github.com/clawpack/clawpack.

# Before running this you might want to do:
#    python $CLAW/clawutil/src/python/clawutil/claw_git_status.py
# and then check the files
#    claw_git_status.txt  and  claw_git_diffs.txt
# to make sure you don't have uncommitted changes in these repositories.

# The commands below check out the master branch of each repository
# and pull any recent changes.

git checkout master
git pull

cd pyclaw
git checkout master
git pull

cd ../classic
git checkout master
git pull

cd ../amrclaw
git checkout master
git pull

cd ../geoclaw
git checkout master
git pull

cd ../clawutil
git checkout master
git pull

cd ../visclaw
git checkout master
git pull

cd ../riemann
git checkout master
git pull

cd ..
