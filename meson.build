project('pyraysum', 'c', 'fortran',
  version : '1.1.0',
  meson_version: '>= 1.1.0',
  default_options : [
    'warning_level=2',
    'buildtype=debugoptimized',
    'fortran_std=legacy',
  ])

py = import('python').find_installation(pure: false)
cc = meson.get_compiler('c')

# Compiler flags
add_project_arguments('-Ofast', '-fallow-argument-mismatch', language : 'fortran')

# Check for numpy
incdir_numpy = run_command(py,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

incdir_f2py = run_command(py,
  ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
  check : true
).stdout().strip()

# Sources
fraysum_srcs = [
  'src/pyraysum/buildmodel.f',
  'src/pyraysum/eigenvec.f',
  'src/pyraysum/eispack-cg.f',
  'src/pyraysum/matrixops.f',
  'src/pyraysum/phaselist.f',
  'src/pyraysum/raysum.f',
  'src/pyraysum/readwrite.f',
  'src/pyraysum/seis-spread.f',
  'src/pyraysum/trace.f'
]

# f2py generation
# We use a custom target to generate the module source
# Note: We assume f2py generates both .c and -f2pywrappers.f
# If it doesn't, we might need to adjust outputs.
f2py_target = custom_target('fraysummodule.c',
  input : fraysum_srcs,
  output : ['fraysummodule.c', 'fraysum-f2pywrappers.f'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', 'fraysum', '--lower'],
)

# Extension
py.extension_module('fraysum',
  [fraysum_srcs, f2py_target, incdir_f2py + '/fortranobject.c'],
  include_directories : include_directories(incdir_numpy, incdir_f2py),
  install : true
)

# Install python sources
install_subdir('pyraysum', install_dir : py.get_install_dir())
