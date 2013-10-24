"""
Execute nosetests in all subdirectories, to run a series of quick
regression tests.

Sends output and result/errors to separate files to simplify checking
results and looking for errors.
"""

import os,sys
import subprocess

# Determine directory:
try:
    CLAW = os.environ['CLAW']
except:
    raise Exception("Need to set CLAW environment variable")

# Remove library files for a clean build:

lib = CLAW + '/amrclaw/src/2d'
cmd = "rm -f %s/*.o %s/*.mod" % (lib,lib)
os.system(cmd)

lib = CLAW + '/geoclaw/src/2d/shallow'
cmd = "rm -f %s/*.o %s/*.mod" % (lib,lib)
os.system(cmd)


def list_tests(tests_dir='.'):
    """
    Searches all subdirectories of tests_dir for tests and prints out a list.
    """
    import os

    current_dir = os.getcwd()
    os.chdir(tests_dir)
    
    dirlist = []
    applist = []

    # Traverse directories depth-first (topdown=False) to insure e.g. that code in
    #    amrclaw/tests/acoustics_2d_radial/1drad 
    # is run before code in
    #    amrclaw/tests/acoustics_2d_radial

    for (dirpath, subdirs, files) in os.walk('.',topdown=False):

        # By convention we assume that a setrun.py file indicates this is an
        # example directory.
        files = os.listdir(os.path.abspath(dirpath))
        if 'setrun.py' in files:
            dirlist.append(os.path.abspath(dirpath))

    os.chdir(current_dir)

    return dirlist
        

def run_tests(tests_dir = '.'):
    import os,sys

    current_dir = os.getcwd()

    dir_list = list_tests(tests_dir)
    # dir_list = ['advection_2d_square']  ## just one for debugging
    print "Found the following test subdirectories:"
    for d in dir_list:
        print "    ", d
 
    print "Will run code and run regression_tests in above subdirectories\n"
    
    fname_output = 'run_tests_output.txt'
    fout = open(fname_output, 'w')
    fout.write("ALL OUTPUT FROM NOSETESTS\n\n")

    fname_errors = 'run_tests_results.txt'
    ferr = open(fname_errors, 'w')
    ferr.write("ALL RESULTS/ERRORS FROM NOSETESTS\n\n")

    os.chdir(tests_dir)

    goodlist_run = []
    badlist_run = []
    
    for directory in dir_list:

        fout.write("\n"+60*"="+"\n")
        fout.write(directory)
        fout.write("\n"+60*"="+"\n")
        ferr.write("\n"+60*"="+"\n")
        ferr.write(directory)
        ferr.write("\n"+60*"="+"\n")

        os.chdir(directory)
        print "Running nosetests in ", directory

        # flush I/O buffers:
        fout.flush()
        ferr.flush()
    
        job = subprocess.Popen(['nosetests'], stdout=fout, stderr=ferr)
        return_code = job.wait()

        if return_code == 0:
            print "No errors encountered\n"
            goodlist_run.append(directory)
        else:
            print "*** nosetest errors encountered: see %s\n" % fname_errors
            badlist_run.append(directory)


    print '------------------------------------------------------------- '
    print ' '
    print 'Successfully ran nosetests in directories:'
    if len(goodlist_run) == 0:
        print '    none'
    for d in goodlist_run:
        print '   ',d
    print ' '
    
    print 'Errors encountered in the following directories:'
    if len(badlist_run) == 0:
        print '    none'
    for d in badlist_run:
        print '   ',d
    print ' '
    
    fout.close()
    ferr.close()
    print 'For all output see ', fname_output
    print 'For results/errors see ', fname_errors

    os.chdir(current_dir)

if __name__=='__main__':
    run_tests()
