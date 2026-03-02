cmake_minimum_required(VERSION 3.17.2...3.29)
project(${SKBUILD_PROJECT_NAME} LANGUAGES C Fortran)

find_package(
  Python
  COMPONENTS Interpreter Development.Module NumPy
  REQUIRED)

# F2PY headers
execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" -c
          "import numpy.f2py; print(numpy.f2py.get_include())"
  OUTPUT_VARIABLE F2PY_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE)

add_library(fortranobject OBJECT "${F2PY_INCLUDE_DIR}/fortranobject.c")

target_link_libraries(fortranobject PUBLIC  Python::NumPy)
target_include_directories(fortranobject PUBLIC "${F2PY_INCLUDE_DIR}")
set_property(TARGET fortranobject PROPERTY POSITION_INDEPENDENT_CODE ON)

add_custom_command(
  OUTPUT fraysummodule.c fraysum-f2pywrappers.f
  DEPENDS src/pyraysum/raysum.f
          src/pyraysum/buildmodel.f
          src/pyraysum/eigenvec.f
          src/pyraysum/eispack-cg.f
          src/pyraysum/matrixops.f
          src/pyraysum/readwrite.f
          src/pyraysum/seis-spread.f
          src/pyraysum/trace.f
          src/pyraysum/params.h
  VERBATIM
  COMMAND "${Python_EXECUTABLE}" -m numpy.f2py
          "${CMAKE_CURRENT_SOURCE_DIR}/src/pyraysum/seis-spread.f" -m fraysum --lower)

python_add_library(
  fraysum MODULE 
  "${CMAKE_CURRENT_BINARY_DIR}/fraysummodule.c"
  "${CMAKE_CURRENT_BINARY_DIR}/fraysum-f2pywrappers.f"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/pyraysum/raysum.f" WITH_SOABI
  "${CMAKE_CURRENT_SOURCE_DIR}/src/pyraysum/buildmodel.f"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/pyraysum/eigenvec.f"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/pyraysum/eispack-cg.f"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/pyraysum/matrixops.f"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/pyraysum/readwrite.f"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/pyraysum/seis-spread.f"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/pyraysum/trace.f"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/pyraysum/params.h")

target_link_libraries(fraysum PRIVATE fortranobject)

install(TARGETS fraysum DESTINATION .)
