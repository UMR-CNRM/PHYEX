#-----------------------------------------------------------------------
# pyfortool, used at configure or build time by the source code
# transformation (phyex_transformation.cmake), the python binding
# generation and the coding norms check (top-level cmake directory)
#
# Usage:
#   phyex_find_pyfortool()
# pyfortool is searched in PATH. If it is not found, or if its version is
# older than the minimal version of PHYEX_PYFORTOOL_SPEC (pyfortool>=x.y.z),
# PHYEX_PYFORTOOL_SPEC is installed in a python virtual environment inside
# the build directory (${CMAKE_CURRENT_BINARY_DIR}/pyfortool). The result is
# kept in the PHYEX_PYFORTOOL cache variable (its version is checked again at
# each configuration). The function defines, in the caller scope:
#   PHYEX_PYFORTOOL        : pyfortool executable (cache variable)
#   PHYEX_PYFORTOOL_VERSION: its version
#   PHYEX_PYFORTOOL_BINDIR : directory of the pyfortool executable (where the
#                            other pyfortool scripts, e.g. pyfortool_parallel,
#                            are found)
#   PHYEX_PYFORTOOL_PYTHON : python interpreter able to import pyfortool
#   PHYEX_PYFORTOOL_IN_BUILD : TRUE if pyfortool is installed in the build
#                              directory (PYTHONPATH must then be unset when
#                              using it; otherwise PYTHONPATH must be kept,
#                              it can be needed to import pyfortool, e.g. in
#                              the isolated environment of a pip build)
#   PHYEX_PYFORTOOL_ENV    : command prefix to run PHYEX_PYFORTOOL_PYTHON or
#                            pyfortool scripts (PATH starting with
#                            PHYEX_PYFORTOOL_BINDIR, PYTHONPATH unset if
#                            needed)
#-----------------------------------------------------------------------

# Default value can be given by an environment variable
# The pyfortool minimal version must be kept in sync with tools/pyproject.toml
# and with the build requirements of the top-level pyproject.toml
set(_default "pyfortool>=0.2.25")
if(DEFINED ENV{PHYEX_PYFORTOOL_SPEC})
  set(_default "$ENV{PHYEX_PYFORTOOL_SPEC}")
endif()
set(PHYEX_PYFORTOOL_SPEC "${_default}" CACHE STRING
    "pip requirement used to install pyfortool if it is not found (or too old)")

#-----------------------------------------------------------------------
# Internal: python interpreter, environment and version of the pyfortool
# executable _exe
#-----------------------------------------------------------------------

function(_phyex_pyfortool_info _exe)
  get_filename_component(_bindir ${_exe} DIRECTORY)
  string(FIND "${_exe}" "${CMAKE_BINARY_DIR}/" _pos)
  if(_pos EQUAL 0)
    set(_in_build TRUE)
    set(_env ${CMAKE_COMMAND} -E env --unset=PYTHONPATH "PATH=${_bindir}:$ENV{PATH}")
  else()
    set(_in_build FALSE)
    set(_env ${CMAKE_COMMAND} -E env "PATH=${_bindir}:$ENV{PATH}")
  endif()

  # Python interpreter able to import pyfortool: the one of the pyfortool
  # script (shebang), the one of its installation directory (virtual
  # environment), or the one found by ecbuild
  set(_candidates)
  file(STRINGS ${_exe} _shebang LIMIT_COUNT 1)
  if(_shebang MATCHES "^#!([^ \t]+)$")
    list(APPEND _candidates ${CMAKE_MATCH_1})
  endif()
  list(APPEND _candidates ${_bindir}/python ${_bindir}/python3 ${Python_EXECUTABLE})
  set(_python)
  set(_version)
  foreach(_candidate ${_candidates})
    if(EXISTS ${_candidate})
      execute_process(COMMAND ${_env} ${_candidate} -c
                              "import pyfortool, importlib.metadata as m; print(m.version('pyfortool'))"
                      RESULT_VARIABLE _result OUTPUT_VARIABLE _version
                      OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
      if(_result EQUAL 0)
        set(_python ${_candidate})
        break()
      endif()
    endif()
  endforeach()
  if(NOT _python)
    ecbuild_critical("No python interpreter able to import pyfortool (${_exe}) found")
  endif()

  set(PHYEX_PYFORTOOL_VERSION ${_version} PARENT_SCOPE)
  set(PHYEX_PYFORTOOL_BINDIR ${_bindir} PARENT_SCOPE)
  set(PHYEX_PYFORTOOL_PYTHON ${_python} PARENT_SCOPE)
  set(PHYEX_PYFORTOOL_IN_BUILD ${_in_build} PARENT_SCOPE)
  set(PHYEX_PYFORTOOL_ENV ${_env} PARENT_SCOPE)
endfunction()

#-----------------------------------------------------------------------
# phyex_find_pyfortool()
#-----------------------------------------------------------------------

function(phyex_find_pyfortool)
  ecbuild_find_python(VERSION 3.8 NO_LIBS REQUIRED)

  # Minimal version required (none if the spec does not give one)
  set(_min_version)
  if(PHYEX_PYFORTOOL_SPEC MATCHES ">=[ \t]*([0-9][0-9.]*)")
    set(_min_version ${CMAKE_MATCH_1})
  endif()

  # pyfortool: found in PATH (or already in cache)
  set(_pyft_venv ${CMAKE_CURRENT_BINARY_DIR}/pyfortool)
  if(PHYEX_PYFORTOOL AND NOT EXISTS ${PHYEX_PYFORTOOL})
    unset(PHYEX_PYFORTOOL CACHE)
  endif()
  find_program(PHYEX_PYFORTOOL pyfortool DOC "pyfortool executable")
  set(_install FALSE)
  if(NOT PHYEX_PYFORTOOL)
    ecbuild_info("pyfortool not found, installing ${PHYEX_PYFORTOOL_SPEC} in ${_pyft_venv}")
    set(_install TRUE)
  else()
    _phyex_pyfortool_info(${PHYEX_PYFORTOOL})
    if(_min_version AND PHYEX_PYFORTOOL_VERSION VERSION_LESS _min_version)
      ecbuild_warn("pyfortool ${PHYEX_PYFORTOOL_VERSION} (${PHYEX_PYFORTOOL}) is too old "
                   "(${PHYEX_PYFORTOOL_SPEC} needed), installing ${PHYEX_PYFORTOOL_SPEC} in ${_pyft_venv}")
      set(_install TRUE)
    endif()
  endif()

  # Otherwise, installation in a python virtual environment inside the
  # build directory
  if(_install)
    file(REMOVE_RECURSE ${_pyft_venv})
    execute_process(COMMAND ${Python_EXECUTABLE} -m venv ${_pyft_venv}
                    RESULT_VARIABLE _result)
    if(_result EQUAL 0)
      execute_process(COMMAND ${CMAKE_COMMAND} -E env --unset=PYTHONPATH
                              ${_pyft_venv}/bin/python -m pip install ${PHYEX_PYFORTOOL_SPEC}
                      RESULT_VARIABLE _result)
    endif()
    if(NOT _result EQUAL 0 OR NOT EXISTS ${_pyft_venv}/bin/pyfortool)
      file(REMOVE_RECURSE ${_pyft_venv})
      ecbuild_critical("Installation of pyfortool failed")
    endif()
    set(PHYEX_PYFORTOOL ${_pyft_venv}/bin/pyfortool CACHE FILEPATH "pyfortool executable" FORCE)
    _phyex_pyfortool_info(${PHYEX_PYFORTOOL})
    if(_min_version AND PHYEX_PYFORTOOL_VERSION VERSION_LESS _min_version)
      ecbuild_critical("pyfortool ${PHYEX_PYFORTOOL_VERSION} installed in ${_pyft_venv} "
                       "does not satisfy ${PHYEX_PYFORTOOL_SPEC}")
    endif()
  endif()

  set(PHYEX_PYFORTOOL_VERSION ${PHYEX_PYFORTOOL_VERSION} PARENT_SCOPE)
  set(PHYEX_PYFORTOOL_BINDIR ${PHYEX_PYFORTOOL_BINDIR} PARENT_SCOPE)
  set(PHYEX_PYFORTOOL_PYTHON ${PHYEX_PYFORTOOL_PYTHON} PARENT_SCOPE)
  set(PHYEX_PYFORTOOL_IN_BUILD ${PHYEX_PYFORTOOL_IN_BUILD} PARENT_SCOPE)
  set(PHYEX_PYFORTOOL_ENV ${PHYEX_PYFORTOOL_ENV} PARENT_SCOPE)
endfunction()
