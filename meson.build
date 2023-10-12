project(
  'clawpack',
  'c', 'fortran',
  version: run_command(
    find_program('python3'),
    files('metadata.py'), 'version',
    check: true,
  ).stdout().strip(),
  license: 'BSD-3-Clause',
  license_files: ['LICENSE'],
  meson_version: '>=1.1.0',
)

fs = import('fs')
py = import('python').find_installation('python3', pure: false)
py_dep = py.dependency()

incdir_numpy = run_command(py,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

incdir_f2py = run_command(py,
    ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()

subpackages = {
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

subdir('clawpack')
foreach subpkg, subdir: subpackages
  subdir(subdir / subpkg)
endforeach
