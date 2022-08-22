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
}

__all__ = list(_subpackages.keys())

_root = _os.path.dirname(__path__[0])
_path = (_os.path.join(_root, *_sdir.split('/'))
         for _sdir in set(_subpackages.values()))
__path__.extend(_sdir for _sdir in _path
                if _os.path.isdir(_sdir))
del _root, _path

__version__ = '5.9.0'   # must also be changed in setup.py
