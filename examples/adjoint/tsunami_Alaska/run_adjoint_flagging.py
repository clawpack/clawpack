"""
    Running example using adjoint method
    """

import os

currentdir = os.getcwd()
adjointdir = currentdir + '/adjoint'

#-------------------------------------------
# Compute solution for adjoint problem
#-------------------------------------------

os.chdir(adjointdir)
os.system('python maketopo.py')
os.system('make new')
os.system('make .plots')

#--------------------------------------------------------------
# Compute solution for forward problem using adjoint-flagging
#--------------------------------------------------------------

os.chdir(currentdir)
os.system('make new')
os.system('python maketopo.py')
os.system('make .plots')
