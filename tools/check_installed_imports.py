#!/usr/bin/env python
# encoding: utf-8
r"""
Assert that every module a Clawpack *install* ships can actually be imported.

WHY THIS EXISTS
---------------
``clawpack/__init__.py`` extends ``__path__`` to the sibling source trees, so
with ``$CLAW`` on ``PYTHONPATH`` every file on disk is importable and a
packaging mistake is invisible.  An install is different: meson-python ships
exactly the files listed in each subpackage's ``meson.build``, and nothing
else exists.  A module that is tracked, imported and *not* listed therefore
works for every developer and fails for every user.

That is not hypothetical.  riemann ``0724f8c`` added
``euler_mapgrid_3D_constants.py`` plus an unconditional import of it in
``riemann/__init__.py``, but not the ``meson.build`` entry, so
``import clawpack.riemann`` raised ``ImportError`` on master for months for
anyone who pip-installed.  Separately, ``riemann_tools.py`` kept importing
``clawpack.visclaw.JSAnimation`` after visclaw deleted it.  Neither was caught:
riemann's workflows only trigger on one topic branch, the doc build imports
through the source shim and never installs, and this repository had no CI at
all.

WHAT IT CHECKS
--------------
One invariant, and only one: *every library module in the installed package
imports, using nothing but the install.*  Third-party dependencies that are
simply absent are reported as skips -- they are an environment property, not a
packaging bug.  A missing ``clawpack.*`` module is a failure, because it can
only mean the install is incomplete.  So is a pre-namespace import such as
``import pyclaw.limiters...``, which resolves in a checkout with ``$CLAW`` on
``PYTHONPATH`` and nowhere else.

Example modules are excluded by default (``--include-examples`` opts in).
They are scripts, not library code: importing
``clawpack.pyclaw.examples.euler_2d.euler_2d`` runs a simulation, and the
example *test* modules load an ``expected_sols.npy`` fixture that no wheel
ships.  Both are worth fixing in pyclaw, neither is the invariant here, and
sweeping them would mean a CI job that solves PDEs to answer a packaging
question.

Run it from a directory that is not a Clawpack checkout, against a real
(non-editable) install::

    pip install .
    cd "$(mktemp -d)"
    python /path/to/clawpack/tools/check_installed_imports.py

``--drift`` is a separate, advisory mode: it lists files that are tracked in
the source tree but shipped by no install.  The import sweep cannot see those
-- a file nobody imports fails nothing -- and that blind spot is how seven
library modules across amrclaw, clawutil and geoclaw stayed unshipped.  It
never fails the build, because some files (``setup.py`` shims, ``conversion/``
templates, test drivers) are deliberately not shipped.
"""

from __future__ import annotations

import argparse
import ast
import contextlib
import importlib
import io
import json
import os
import re
import subprocess
import sys
import tempfile
import traceback


# Must precede any matplotlib import: several visclaw modules call
# matplotlib.use()/rc() at import time and would otherwise want a display.
os.environ.setdefault('MPLBACKEND', 'Agg')

TOOLS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(TOOLS_DIR)
SCRIPT_NAME = os.path.basename(os.path.abspath(__file__))

# clawpack.* modules that an install is not expected to provide.  Empty by
# design: everything in the meson project ships, and the one genuinely
# external subpackage (clawpack.dclaw, its own repository) is imported inside a
# function, so a module-level sweep never reaches it.  An entry here is a
# standing admission that an installed module cannot be imported, so it needs
# the same justification as a new line in the doc warning baseline.
OPTIONAL_CLAWPACK: dict[str, str] = {}

# Modules that are shipped and cannot be imported, with the reason.  Same
# bargain as tools/doc_warnings_baseline.txt in the doc repository: existing
# breakage is listed explicitly so that *new* breakage still fails the build.
# A module that starts importing must be removed from here -- a baseline
# nobody prunes stops describing anything.
KNOWN_BROKEN = {
    'clawpack.pyclaw.limiters.reconstruct':
        'its one import is pre-namespace: `import pyclaw.limiters.weno.'
        'reconstruct`, which resolves nowhere since the clawpack namespace.  '
        'The C extension it wants IS built and installed, so the import line '
        'is the whole defect -- but sharpclaw/solver.py falls back to '
        'limiters.recon on ImportError, and only that fallback defines '
        'weno5_wave(), so simply prefixing `clawpack.` would break '
        'char_decomp=1.  A pyclaw decision, and a numerical one.',
}

# Example modules are scripts: importing one can run a simulation.  Excluded
# unless --include-examples; see WHAT IT CHECKS above.
_EXAMPLE_MARKER = '.examples.'


def _installed_clawpack():
    """Import clawpack, or explain why this check cannot be trusted.

    Returns the module, or ``None`` after printing why we refuse to run.  The
    refusal matters more than the sweep: if the source tree is importable the
    sweep passes unconditionally and guards nothing, which is precisely how
    the original bug survived in riemann's own test workflows.
    """
    try:
        import clawpack
    except ImportError as exc:
        print(f'FAIL: cannot import clawpack at all: {exc}\n\n'
              'Install it first (`pip install .`); this check inspects an '
              'install, not a source tree.', file=sys.stderr)
        return None

    path = os.path.abspath(clawpack.__file__)
    if os.path.commonpath([path, REPO_ROOT]) == REPO_ROOT:
        print(f'FAIL: `import clawpack` resolved to\n  {path}\n'
              f'which is inside the source tree at\n  {REPO_ROOT}\n\n'
              'This check is meaningless against the source tree: the shim in '
              'clawpack/__init__.py\nextends __path__ to the sibling '
              'checkouts, so every file on disk imports and no\nmissing '
              'meson.build entry can ever be detected.  Install the project '
              'and re-run\nfrom a directory outside the checkout:\n\n'
              '    pip install .\n'
              '    cd "$(mktemp -d)"\n'
              f'    python {os.path.join(TOOLS_DIR, SCRIPT_NAME)}',
              file=sys.stderr)
        return None
    return clawpack


def installed_modules(clawpack, include_examples: bool = False) -> list[str]:
    """Every importable dotted name the install actually contains.

    Walks the installed package directory rather than reading the meson
    install plan: it is the user's-eye view, and it works for a plain
    ``pip install .`` with no build directory left behind.
    """
    root = os.path.dirname(os.path.abspath(clawpack.__file__))
    names = ['clawpack']
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if d != '__pycache__']
        rel = os.path.relpath(dirpath, root)
        parts = [] if rel == '.' else rel.split(os.sep)
        if not all(p.isidentifier() for p in parts):
            dirnames[:] = []          # not reachable by import; skip subtree
            continue
        for name in sorted(filenames):
            if not name.endswith('.py') or name == '__init__.py':
                continue
            stem = name[:-3]
            if not stem.isidentifier():
                continue
            names.append('.'.join(['clawpack'] + parts + [stem]))
        if parts and os.path.exists(os.path.join(dirpath, '__init__.py')):
            names.append('.'.join(['clawpack'] + parts))
    if not include_examples:
        names = [n for n in names if _EXAMPLE_MARKER not in n]
    return sorted(set(names))


def missing_name(exc: BaseException) -> str | None:
    """The dotted name Python could not find, or ``None``.

    Four shapes, in order of reliability:
    ``ModuleNotFoundError.name``; ``ImportError.name`` + ``.name_from`` for a
    failed ``from . import X``; the ``__cause__``/``__context__`` chain, for a
    hand-raised ImportError that re-reports a missing optional dependency
    (geoclaw's xarray_backends does this, which is why it chains the original);
    and finally the message text.
    """
    name = getattr(exc, 'name', None)
    if (isinstance(exc, ImportError)
            and not isinstance(exc, ModuleNotFoundError)):
        name_from = getattr(exc, 'name_from', None)
        if name and name_from:
            return f'{name}.{name_from}'
    if name:
        return name

    for chained in (exc.__cause__, exc.__context__):
        if isinstance(chained, ImportError):
            found = missing_name(chained)
            if found:
                return found

    text = str(exc)
    match = re.search(r"cannot import name '([^']+)' from '?([\w.]+)'?", text)
    if match:
        return f'{match.group(2)}.{match.group(1)}'
    match = re.search(r"No module named '([\w.]+)'", text)
    if match:
        return match.group(1)
    return None


def _subpackage_dirs() -> dict[str, str]:
    """``{'riemann': 'riemann', 'geoclaw': 'geoclaw/src/python', ...}``.

    Read from the shim in the checkout by parsing, not importing -- the whole
    point of this script is to keep the source tree out of ``sys.modules``.
    """
    shim = os.path.join(REPO_ROOT, 'clawpack', '__init__.py')
    try:
        with open(shim, encoding='utf-8') as fh:
            tree = ast.parse(fh.read())
    except (OSError, SyntaxError):
        return {}
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(
                getattr(t, 'id', None) == '_subpackages'
                for t in node.targets):
            try:
                return ast.literal_eval(node.value)
            except ValueError:
                return {}
    return {}


def _meson_hint(dotted: str) -> str | None:
    """The ``meson.build`` most likely needing a ``python_sources`` line."""
    parts = dotted.split('.')
    if len(parts) < 2 or parts[0] != 'clawpack':
        return None
    subdir = _subpackage_dirs().get(parts[1])
    if subdir is None:
        return None
    return '/'.join([subdir, parts[1], 'meson.build'])


def classify(exc: BaseException) -> tuple[str, str | None]:
    """``(kind, missing name)`` for an import failure.

    ``missing-clawpack``  a module the install should contain but does not.
    ``legacy-namespace``  an unqualified ``pyclaw.``/``geoclaw.``/... import.
    ``third-party``       an absent optional dependency; not our problem.
    ``unknown``           nothing identifiable; reported with a traceback.
    """
    found = missing_name(exc)
    if found is None:
        return 'unknown', None
    head = found.split('.')[0]
    if head == 'clawpack':
        if any(found == prefix or found.startswith(prefix + '.')
               for prefix in OPTIONAL_CLAWPACK):
            return 'third-party', found
        return 'missing-clawpack', found
    if head in _subpackage_dirs():
        return 'legacy-namespace', found
    return 'third-party', found


def sweep(names: list[str], verbose: bool = False) -> int:
    """Import everything; report; return an exit code."""
    failures: list[tuple[str, BaseException, str, str, str | None]] = []
    skips: dict[str, list[str]] = {}
    known: list[str] = []
    stale: list[str] = []
    imported = 0

    # A module that reads stdin (clawutil.run_examples used to, at import)
    # must fail fast rather than hang a CI job forever.
    devnull = open(os.devnull)
    stdin, sys.stdin = sys.stdin, devnull

    try:
        for name in names:
            noise = io.StringIO()
            try:
                with contextlib.redirect_stdout(noise), \
                        contextlib.redirect_stderr(noise):
                    importlib.import_module(name)
            except (Exception, SystemExit) as exc:
                kind, found = classify(exc)
                if name in KNOWN_BROKEN:
                    known.append(name)
                elif kind == 'third-party':
                    skips.setdefault(found.split('.')[0], []).append(name)
                else:
                    failures.append((name, exc, noise.getvalue(), kind, found))
            else:
                imported += 1
                if name in KNOWN_BROKEN:
                    stale.append(name)
                if verbose:
                    print(f'  ok    {name}')
    finally:
        sys.stdin = stdin
        devnull.close()

    print(f'\nok   {imported} module(s) imported from the install')
    if skips:
        print(f'skip {sum(len(v) for v in skips.values())} module(s) needing '
              'absent third-party packages:')
        for dep in sorted(skips):
            mods = skips[dep]
            shown = ', '.join(mods[:3]) + (f', +{len(mods) - 3} more'
                                           if len(mods) > 3 else '')
            print(f'       {dep:<12} {shown}')
    for name in known:
        print(f'known-broken  {name}\n              {KNOWN_BROKEN[name]}')

    if stale:
        print('\nFAIL: KNOWN_BROKEN lists module(s) that now import fine.\n'
              'Delete them from the list -- a baseline nobody prunes stops '
              'describing anything.\n', file=sys.stderr)
        for name in stale:
            print(f'  {name}', file=sys.stderr)

    if not failures:
        return 1 if stale else 0

    print(f'\nFAIL: {len(failures)} installed module(s) cannot be imported '
          'from the install.\n'
          'Every one of these works from a source checkout, so only users see '
          'them.\n', file=sys.stderr)
    for name, exc, noise, kind, found in failures:
        print(f'  {name}', file=sys.stderr)
        print('      ' + ''.join(traceback.format_exception_only(
            type(exc), exc)).strip(), file=sys.stderr)
        if kind == 'missing-clawpack':
            hint = _meson_hint(found)
            print(f'      missing module: {found}', file=sys.stderr)
            if hint:
                print(f'      -> add it to python_sources in {hint}, or drop '
                      'the import if the\n         module is gone',
                      file=sys.stderr)
        elif kind == 'legacy-namespace':
            print(f'      pre-namespace import of {found!r}\n'
                  f'      -> it resolves only with $CLAW on PYTHONPATH; write '
                  f'clawpack.{found}', file=sys.stderr)
        else:
            print('      could not identify a missing module; full traceback:',
                  file=sys.stderr)
            for line in traceback.format_exception(
                    type(exc), exc, exc.__traceback__)[-4:]:
                print('      ' + line.rstrip(), file=sys.stderr)
        if noise.strip():
            print('      module output: ' + noise.strip()[:200],
                  file=sys.stderr)
    return 1


def _install_plan(build_dir: str | None) -> set[str]:
    """Absolute source paths of every ``.py`` meson would install."""
    if build_dir is None:
        build_dir = tempfile.mkdtemp(prefix='clawpack-drift-')
        subprocess.run(['meson', 'setup', build_dir], cwd=REPO_ROOT,
                       check=True, stdout=subprocess.DEVNULL)
    plan = json.loads(subprocess.run(
        ['meson', 'introspect', build_dir, '--install-plan'],
        cwd=REPO_ROOT, check=True, capture_output=True, text=True).stdout)
    return {os.path.abspath(os.path.join(REPO_ROOT, src))
            for entries in plan.values() for src, info in entries.items()
            if info.get('destination', '').endswith('.py')}


def drift(build_dir: str | None) -> int:
    """Report tracked-but-unshipped modules.  Advisory: always returns 0."""
    declared = _install_plan(build_dir)
    print('Tracked .py files that no install ships:\n')
    total = 0
    for pkg, subdir in sorted(_subpackage_dirs().items()):
        pkg_dir = os.path.join(REPO_ROOT, *subdir.split('/'), pkg)
        if not os.path.isdir(pkg_dir):
            continue                    # optional repo, e.g. dclaw
        repo = subdir.split('/')[0]
        listed = subprocess.run(['git', 'ls-files', '--', '*.py'],
                                cwd=pkg_dir, capture_output=True, text=True)
        undeclared = sorted(
            rel for rel in listed.stdout.split()
            if os.path.abspath(os.path.join(pkg_dir, rel)) not in declared)
        if undeclared:
            total += len(undeclared)
            print(f'  {repo}/  ({pkg})')
            for rel in undeclared:
                print(f'      {rel}')
    print(f'\n{total} file(s).  Advisory only -- setup.py shims, conversion '
          'templates and test\ndrivers are deliberately unshipped.  A library '
          'module here is a bug: it is\nimportable for every developer and '
          'missing for every user.')
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--verbose', action='store_true',
                        help='print every module as it is imported')
    parser.add_argument('--include-examples', action='store_true',
                        help='also import example modules (they run '
                             'simulations; see the module docstring)')
    parser.add_argument('--drift', action='store_true',
                        help='instead, list tracked files no install ships')
    parser.add_argument('--build-dir',
                        help='existing meson build directory for --drift')
    args = parser.parse_args(argv)

    if args.drift:
        return drift(args.build_dir)

    clawpack = _installed_clawpack()
    if clawpack is None:
        return 2

    print(f'clawpack {getattr(clawpack, "__version__", "?")} from '
          f'{os.path.dirname(os.path.abspath(clawpack.__file__))}')
    names = installed_modules(clawpack, include_examples=args.include_examples)
    print(f'checking {len(names)} installed module(s)'
          + ('' if args.include_examples else ', examples excluded'))
    return sweep(names, verbose=args.verbose)


if __name__ == '__main__':
    sys.exit(main())
