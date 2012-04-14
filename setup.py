#!/usr/bin/env python

try:
    from numpy.distutils.core import setup
except Exception, err:
    sys.stderr.write('Unable to import numpy, please install numpy first before installing clawpack!')

from numpy.distutils.command.build import build
from numpy.distutils.exec_command import exec_command

import os.path

class build_with_submodules(build):
    def run(self):
        if os.path.exists('.git'):
            exec_command(['git', 'submodule', 'init'])
            exec_command(['git', 'submodule', 'update'])
        build.run(self)

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration(None, parent_package, top_path)
    config.add_subpackage('clawutil')
    config.add_subpackage('riemann')
    config.add_subpackage('visclaw')
    config.add_subpackage('pyclaw')
    return config

if __name__ == '__main__':
    setup(cmdclass={"build": build_with_submodules},
          **configuration(top_path='').todict())
