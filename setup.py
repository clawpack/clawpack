#!/usr/bin/env python

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration(None, parent_package, top_path)
    config.add_subpackage('clawutil')
    config.add_subpackage('riemann')
    config.add_subpackage('visclaw')
    config.add_subpackage('pyclaw')
    return config

if __name__ == '__main__':
    try:
        from numpy.distutils.core import setup
    except Exception, err:
        sys.stderr.write('Unable to import numpy, please install numpy first before installing clawpack!')
    setup(**configuration(top_path='').todict())
