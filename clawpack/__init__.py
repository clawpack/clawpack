
try:
    import clawutil
    import pyclaw
    import petclaw
    import riemann
    import visclaw

    __all__ = ['clawutil','pyclaw','petclaw','riemann','visclaw']
except ImportError:
    pass

# need top-level try/except protection because this package may be accessed mid-install
