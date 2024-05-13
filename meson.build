project(
  'clawpack',
  'c', 'fortran',
  version: run_command([find_program('metadata.py'), 'version'], check: true).stdout().strip(),
  license: 'BSD-3-Clause',
  license_files: ['LICENSE'],
  meson_version: '>=1.1.0',
)

fs = import('fs')
py = import('python').find_installation(pure: false)
py_dep = py.dependency()

incdir_numpy = run_command(py,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

incdir_f2py = run_command(py,
    ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()

inc_np = include_directories(incdir_numpy, incdir_f2py)

f2py_version = run_command([py, '-m', 'numpy.f2py', '-v'], check: true).stdout().strip()
message(f'f2py version: @f2py_version@')

f2py = [
  py, '-m', 'numpy.f2py',
  '--lower', '--quiet',
  '--build-dir', '@OUTDIR@',
  '@INPUT@',
  '-m', # <extension-name>
]

fc = meson.get_compiler('fortran')
cc = meson.get_compiler('c')

add_project_arguments(
  fc.get_supported_arguments(
    '-Wno-unused-dummy-argument',
    '-Wno-unused-variable',
    '-Wno-unused-label',
  ),
  language : 'fortran',
)
add_project_arguments(
#  '-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION',
  cc.get_supported_arguments(
    '-Wno-misleading-indentation',
    '-Wno-unused-but-set-variable',
  ),
  language : 'c',
)

add_project_link_arguments(
  cc.get_supported_link_arguments('-Wl,-no_compact_unwind'),
  language: 'c',
)

subpackages = {
  'clawutil':   'clawutil/src/python',
  'pyclaw':     'pyclaw/src',
  'riemann':    'riemann',
  'visclaw':    'visclaw/src/python',
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
