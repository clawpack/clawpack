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

CLAW_dir = _os.environ.get('CLAW')
if CLAW_dir == None:
    raise Exception('You must set the CLAW environment to use an editable install.')

_root = _os.path.dirname(_os.path.join(CLAW_dir, 'clawpack'))
_path = (_os.path.join(_root, *_sdir.split('/'))
         for _sdir in set(_subpackages.values()))
__path__.extend(_sdir for _sdir in _path
                if _os.path.isdir(_sdir))
del _root, _path

__version__ = '5.9.1'   # must also be changed in setup.py
