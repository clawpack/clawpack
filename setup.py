"""Clawpack: Python-based Clawpack installer
"""

# much of the functionality of this file was taken from the scipy setup.py script.

DOCLINES = __doc__.split("\n")

import os
import sys
import warnings
import subprocess
import shutil
import re

if sys.version_info[0] < 3:
    import __builtin__ as builtins
else:
    import builtins

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

MAJOR               = 0
MINOR               = 1
MICRO               = 0
ISRELEASED          = False
VERSION             = '%d.%d.%d' % (MAJOR, MINOR, MICRO)

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


# BEFORE importing distutils, remove MANIFEST. distutils doesn't properly
# update it when the contents of directories change.
if os.path.exists('MANIFEST'):
    os.remove('MANIFEST')

def write_version_py(filename=version_file_path):
    cnt = """
# THIS FILE IS GENERATED FROM CLAWPACK SETUP.PY
short_version = '%(version)s'
version = '%(version)s'
full_version = '%(full_version)s'
git_revision = '%(git_revision)s'
release = %(isrelease)s

if not release:
    version = full_version
"""
    # Adding the git rev number needs to be done inside
    # write_version_py(), otherwise the import of scipy.version messes
    # up the build under Python 3.
    FULLVERSION = VERSION
    if os.path.exists('.git'):
        GIT_REVISION = git_version()
    elif os.path.exists(version_file_path):
        # must be a source distribution, use existing version file
        from clawpack.version import git_revision as GIT_REVISION
    else:
        GIT_REVISION = "Unknown"

    if not ISRELEASED:
        FULLVERSION += '.dev-' + GIT_REVISION[:7]

    a = open(filename, 'w')
    try:
        a.write(cnt % {'version': VERSION,
                       'full_version' : FULLVERSION,
                       'git_revision' : GIT_REVISION,
                       'isrelease': str(ISRELEASED)})
    finally:
        a.close()    

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration(None, parent_package, top_path)

    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True)

    config.add_subpackage('clawpack')
    config.get_version(os.path.join('clawpack','version.py'))
    return config


def setup_package():
    import sys
    # Rewrite the version file everytime

    old_path = os.getcwd()
    local_path = os.path.dirname(os.path.abspath(sys.argv[0]))
    src_path = local_path

    os.chdir(local_path)
    sys.path.insert(0, local_path)
    sys.path.insert(0, os.path.join(local_path, 'clawpack'))  # to retrieve version

    old_path = os.getcwd()
    os.chdir(src_path)
    sys.path.insert(0, src_path)

    write_version_py()

    setup_dict = dict(
        name = 'clawpack',
        maintainer = "Clawpack Developers",
        maintainer_email = "claw-dev@googlegroups.com",
        description = DOCLINES[0],
        long_description = "\n".join(DOCLINES[2:]),
        url = "http://www.clawpack.org",
        download_url = "http://github.com/clawpack/clawpack/tarball/master#egg=clawpack-dev", 
        license = 'BSD',
        classifiers=[_f for _f in CLASSIFIERS.split('\n') if _f],
        platforms = ["Linux", "Solaris", "Mac OS-X", "Unix"],
        )

    try:
        if 'egg_info' in sys.argv:
            # only egg information for downloading requirements
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
            return

        if os.path.exists('.git'):
            from numpy.distutils.exec_command import exec_command
            exec_command(['git', 'submodule', 'init'])
            fails = 0
            while fails < 20 and exec_command(['git', 'submodule', 'update'])[1]:
                fails = fails+1
                import time
                print "having difficulties updating submodules, waiting 5s and trying again [fail %d/20]" % fails
                time.sleep(5)
            # *always* need these
            # now build symbolic links to repositories
            if not os.path.exists('clawpack/clawutil'):
                os.symlink(os.path.abspath('clawutil/src/python/clawutil'),
                           'clawpack/clawutil')
            if not os.path.exists('clawpack/riemann'):
                os.symlink(os.path.abspath('riemann/src/python/riemann'),
                           'clawpack/riemann')
                # need this one to build Fortran sources naturally
            if not os.path.exists('clawpack/riemann/src'):
                os.symlink(os.path.abspath('riemann/src'),
                           'clawpack/riemann/src')
            if not os.path.exists('clawpack/visclaw'):
                os.symlink(os.path.abspath('visclaw/src/python/visclaw'),
                           'clawpack/visclaw')
            if not os.path.exists('clawpack/pyclaw'):
                os.symlink(os.path.abspath('pyclaw/src/pyclaw'),
                           'clawpack/pyclaw')
            if not os.path.exists('clawpack/petclaw'):
                os.symlink(os.path.abspath('pyclaw/src/petclaw'),
                           'clawpack/petclaw')

            from numpy.distutils.core import setup
            setup(configuration=configuration,
                  **setup_dict)

    except Exception as err:
        print err
        raise err
    finally:
        del sys.path[0]
        os.chdir(old_path)
    return

if __name__ == '__main__':
    setup_package() 
