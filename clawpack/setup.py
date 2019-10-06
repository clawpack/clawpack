
def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('clawpack',parent_package,top_path)
    config.add_subpackage('clawutil',  subpackage_path='clawutil/src/python/clawutil')
    config.add_subpackage('riemann', subpackage_path='riemann/riemann')
    config.add_subpackage('visclaw',  subpackage_path='visclaw/src/python/visclaw')
    config.add_subpackage('pyclaw', subpackage_path='pyclaw/src/pyclaw')
    config.add_subpackage('petclaw', subpackage_path='pyclaw/src/petclaw')
    config.add_subpackage('forestclaw', subpackage_path='pyclaw/src/forestclaw')
    config.add_subpackage('classic',  subpackage_path='classic/src/python/classic')
    config.add_subpackage('amrclaw',  subpackage_path='amrclaw/src/python/amrclaw')
    config.add_subpackage('geoclaw',  subpackage_path='geoclaw/src/python/geoclaw')
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
