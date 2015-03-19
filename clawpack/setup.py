
def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('clawpack',parent_package,top_path)
    config.add_subpackage('clawutil')
    config.add_subpackage('riemann')
    config.add_subpackage('visclaw')
    config.add_subpackage('pyclaw')
    config.add_subpackage('petclaw')
    config.add_subpackage('classic')
    config.add_subpackage('amrclaw')
    config.add_subpackage('geoclaw')
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
