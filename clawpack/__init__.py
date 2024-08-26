import os as _os

_subpackages = {
    'clawutil':   'clawutil/src/python',
    'riemann':    'riemann',
    'visclaw':    'visclaw/src/python',
    'pyclaw':     'pyclaw/src',
    'petclaw':    'pyclaw/src',
    'forestclaw': 'pyclaw/src',
    'classic':    'classic/src/python',
    'amrclaw':    'amrclaw/src/python',
    'geoclaw':    'geoclaw/src/python',
    'dclaw':      'dclaw/src/python',
}

__all__ = list(_subpackages.keys())

_init = _os.path.abspath(__file__)
_root = _os.path.dirname(_os.path.dirname(_init))
_path = (
    _os.path.join(_root, *_sdir.split('/'))
    for _sdir in set(_subpackages.values())
)
__path__.extend(filter(_os.path.isdir, _path))
del _os, _init, _root, _path

__version__ = '5.11.0'   # must also be changed in setup.py
