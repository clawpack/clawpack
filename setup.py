"""Clawpack: Python-based Clawpack installer

This installer provides:

* automated clawpack subpackage developer environment setup 
* installation of Python modules from clawpack subpackages 
  in the clawpack.package namespace

Please send an email to claw-dev@googlegroups.com for any general questions
or raise installation-related issues and pull requests to our GitHub repository:

http://github.com/clawpack/clawpack
"""

# some of the functionality of this file is reused from the SciPy setup.py script.

from __future__ import absolute_import
from __future__ import print_function
DOCLINES = __doc__.split("\n")

import os
import sys
import warnings
import subprocess
import shutil
import re
import time
from contextlib import contextmanager

join = os.path.join

# Specify top-level subpackages in SUBPACKAGES.  SUBPACKAGES is a
# dictionary specifying which packages you would like installed.  By
# default, the installer will download all of these packages for you.
# You can disable packages by deleting from the dictionary, but you
# will probably need at least pyclaw, visclaw, clawutil, and riemann.

# ADVICE TO DEVELOPERS:
# The 'python_src_dir' dictionary specifies the symlink directory
# structure provided to enable the clawpack.xxx namespace.
# For example, the pyclaw Python package lives in pyclaw/src/pyclaw.
# It is made available in clawpack/pyclaw by creating a symlink
# specified by 'src' (the package name is implicit at the beginning
# and end).  Getting the pyclaw/examples directory requires a small feat
# of gymnastics, see the code in dev_setup() if you would like to
# refactor this.

SUBPACKAGES = {
    'amrclaw': {
        'python_src_dir': [('amrclaw', join('src', 'python'))]
    },
    'clawutil': {
        'python_src_dir': [('clawutil', join('src', 'python'))]
    },                
    'geoclaw': {
        'python_src_dir': [('geoclaw', join('src', 'python'))]
    },                
    'classic': {
        'python_src_dir': [('classic', join('src', 'python'))]
    },
    'pyclaw': {
        'python_src_dir': [('pyclaw', 'src'),
                           ('petclaw', 'src'),
                           ('forestclaw', 'src'),
                           (join('pyclaw','examples'), '..')]
    },
    'riemann': {
        'python_src_dir': [(None,'src')]
    },
    'visclaw': {
        'python_src_dir': [('visclaw', join('src', 'python'))]
    },            
}

#########################
### BEGIN BOILERPLATE ### 
#########################
    
CLASSIFIERS = """\
Development Status :: 4 - Beta
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved
Programming Language :: C
Programming Language :: Python
Topic :: Software Development
Topic :: Scientific/Engineering
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS

"""  

MAJOR               = 5
MINOR               = 5
MICRO               = 0
TYPE                = ''
VERSION             = '%d.%d.%d%s' % (MAJOR, MINOR, MICRO, TYPE)

package_path       = os.path.join(os.path.dirname(__file__),'clawpack')

version_file_path  = os.path.join(package_path,'version.py')

# Return the git revision as a string
def git_version():
    def _minimal_ext_cmd(cmd):
        # construct minimal environment
        env = {}
        for k in ['SYSTEMROOT', 'PATH']:
            v = os.environ.get(k)
            if v is not None:
                env[k] = v
        # LANGUAGE is used on win32
        env['LANGUAGE'] = 'C'
        env['LANG'] = 'C'
        env['LC_ALL'] = 'C'
        out = subprocess.Popen(cmd, stdout = subprocess.PIPE, env=env).communicate()[0]
        return out

    try:
        out = _minimal_ext_cmd(['git', 'rev-parse', 'HEAD'])
        GIT_REVISION = out.strip().decode('ascii')
    except OSError:
        GIT_REVISION = "Unknown"

    return GIT_REVISION

def write_version_py(filename=version_file_path):

    old_path = os.getcwd()
    local_path = os.path.dirname(os.path.abspath(sys.argv[0]))
    src_path = local_path

    os.chdir(local_path)
    sys.path.insert(0, local_path)
    sys.path.insert(0, os.path.join(local_path, 'clawpack'))  # to retrieve version

    old_path = os.getcwd()
    os.chdir(src_path)
    sys.path.insert(0, src_path)

    cnt = """
# THIS FILE IS GENERATED FROM CLAWPACK SETUP.PY
short_version = '%(version)s'
version = '%(version)s'
full_version = '%(full_version)s'
git_revision = '%(git_revision)s'
"""
    # Adding the git rev number needs to be done inside
    # write_version_py(), otherwise the import of clawpack.version messes
    # up the build under Python 3.

    FULLVERSION = VERSION
    if os.path.exists('.git'):
        GIT_REVISION = git_version()
    elif os.path.exists(version_file_path):
        # must be a source distribution, use existing version file
        from clawpack.version import git_revision as GIT_REVISION
    else:
        GIT_REVISION = "Unknown"

    with open(filename, 'w') as a:
        a.write(cnt % {'version': VERSION,
                       'full_version' : FULLVERSION,
                       'git_revision' : GIT_REVISION})

    del sys.path[0]
    os.chdir(old_path)

#########################
###  END BOILERPLATE  ### 
#########################
    
def symlink(src, target):
    """ symlinks src to target if target does not already exist

    Both paths may be relative (they are parsed through os.path.abspath)
    """
    src = os.path.abspath(src)
    target = os.path.abspath(target)

    if not os.path.exists(src):
        raise IOError("trying to symlink %s: which does not exist" % (src))
    
    if not os.path.exists(target):
        os.symlink(os.path.abspath(src), target)

def unsymlink(target):
    """ unsymlinks target if it exists
    """

    if os.path.exists(target):
        os.unlink(target)
    
def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration(None, parent_package, top_path)

    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True,
                       quiet=True)

    config.add_subpackage('clawpack')
    config.get_version(os.path.join('clawpack','version.py'))
    return config


def dev_setup(subpackages):
    """clawpack developer environment setup 
    
user has a .git subdirectory, assume they want us to set up submodules for them.

For each package in subpackages:
if the package directory does not exist or is empty, calls: 
    git submodule init <package> 
    git submodule update <package> 

with timeouts for update, which may be over a fickle remote connection

After each package is checked out, build symbolic links to ./clawpack/package
which allows for a consistent clawpack.package namespace.
"""
    if not os.path.exists('.git'):
        raise Exception("Developer setup requested but top-level clawpack" + \
                        " is not a git repository")

    for package, package_dict in subpackages.items():
        if not os.path.exists(package) or not (os.listdir(package)):
            subprocess.check_call(['git', 'submodule', 'init', package])

            fails = 0
            while fails < 20 and subprocess.call(['git', 'submodule', 'update', 
                                                  package]):
                fails = fails+1
                print("having difficulties updating submodules," + \
                  "waiting 5s and trying again [fail %d/20]" % fails)
                time.sleep(5)

        print("Git development environment initialized for:", package)


def make_symlinks(subpackages):
    for package, package_dict in subpackages.items():
        for subpackage, src_dir in package_dict['python_src_dir']:
            if subpackage:
                symlink(os.path.join(package, src_dir, subpackage),
                        os.path.join('clawpack', subpackage))
            else:
                symlink(os.path.join(package, src_dir),
                        os.path.join('clawpack', package))
                

@contextmanager
def stdout_redirected(new_stdout='install.log'):
    """This redirects stdout for this processes and those forked from it.
    Avoids many pages of warnings generated by f2py being printed to the screen.
    """

    old = os.dup(1)
    os.close(1)
    os.open(new_stdout, os.O_WRONLY | os.O_CREAT)

    try:
        yield None
    finally:
        os.close(1)
        os.dup2(old,1)
        os.close(old)


def setup_package(setup_dict, subpackages, symlink_only=False):
    from numpy.distutils.core import setup

    # Rewrite the version file every time we install
    write_version_py()

    with stdout_redirected(): # Don't print f2py warnings to screen

        # we may end up mucking with symbolic path links for the install 
        # to support a consistent clawpack.package namespace
        # the finally clause here undoes a potentially dangerous 
        # recursive symbolic link that is needed for the numpy.distutils
        # machinery to properly understand some Fortran source paths
        if os.path.exists('.git'):
            dev_setup(subpackages)
        make_symlinks(subpackages)
        if not symlink_only:
            setup(configuration=configuration,
                  **setup_dict)


if __name__ == '__main__':
    setup_dict = dict(
        name = 'clawpack',
        maintainer = "Clawpack Developers",
        maintainer_email = "claw-dev@googlegroups.com",
        description = DOCLINES[0],
        long_description = "\n".join(DOCLINES[2:]),
        url = "http://www.clawpack.org",
        download_url = "git+git://github.com/clawpack/clawpack.git#egg=clawpack-dev", 
        license = 'BSD',
        classifiers=[_f for _f in CLASSIFIERS.split('\n') if _f],
        platforms = ["Linux", "Solaris", "Mac OS-X", "Unix"],
        )

    # python setup.py git-dev sets up subpackages
    if 'git-dev' in sys.argv:
        # not a real install
        dev_setup(SUBPACKAGES)
        make_symlinks(SUBPACKAGES)
    # egg_info requests only provide install requirements
    # this is how "pip install clawpack" installs numpy correctly.
    elif 'egg_info' in sys.argv:
        # not a real install
        from setuptools import setup
        setuptools_dict = dict(
            install_requires = ['numpy >= 1.6',
                                'matplotlib >= 1.0.1',
                                ],                            
            extras_require = {'petclaw': ['petsc4py >= 1.2'],
                              'euler'  : ['scipy >= 0.10.0']},
            )
        setup_dict.update(setuptools_dict)
        setup(**setup_dict)
    else:
        # okay, real install
        if 'symlink-only' in sys.argv:
            setup_package(setup_dict, SUBPACKAGES, symlink_only=True) 
        else:
            setup_package(setup_dict, SUBPACKAGES) 
