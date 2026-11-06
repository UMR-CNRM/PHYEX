#-----------------------------------------------------------------------
# PHYEX top-level CMakeLists.txt
#
# Builds the PHYEX library (src/phyex) together with its dependencies:
#   - ecbuild : used if found, otherwise downloaded at configure time
#   - fiat    : used if found, otherwise downloaded and built as a subproject
#
# Dependency versions and repositories can be overridden with environment
# variables (read at first configure) or with -D<VAR>=... on the command line:
#   PHYEX_ECBUILD_VERSION, PHYEX_ECBUILD_REPO, PHYEX_FIAT_VERSION, PHYEX_FIAT_REPO
#
# Usage example:
#   cmake -S . -B build/cmake -DCMAKE_TOOLCHAIN_FILE=<arch>.ecbuild
#   cmake --build build/cmake -j
#-----------------------------------------------------------------------

cmake_minimum_required(VERSION 3.18)

#-----------------------------------------------------------------------
# Dependency versions and source repositories
#-----------------------------------------------------------------------

# Defines a cache variable whose default value comes from the environment
# variable of the same name when it exists
macro(phyex_dependency_setting VAR DEFAULT DOC)
  if(DEFINED ENV{${VAR}})
    set(${VAR} "$ENV{${VAR}}" CACHE STRING "${DOC}")
  else()
    set(${VAR} "${DEFAULT}" CACHE STRING "${DOC}")
  endif()
endmacro()

phyex_dependency_setting(PHYEX_ECBUILD_VERSION "3.12.0"
                         "ecbuild version (git tag) to download if ecbuild is not found")
phyex_dependency_setting(PHYEX_ECBUILD_REPO "https://github.com/ecmwf/ecbuild"
                         "ecbuild git repository")
phyex_dependency_setting(PHYEX_FIAT_VERSION "1.6.2"
                         "fiat version (git tag) to download if fiat is not found")
phyex_dependency_setting(PHYEX_FIAT_REPO "https://github.com/ecmwf-ifs/fiat.git"
                         "fiat git repository")

#-----------------------------------------------------------------------
# Step 1: ecbuild
#-----------------------------------------------------------------------

# Look for an existing installation (ecbuild_DIR, CMAKE_PREFIX_PATH, PATH...)
find_package(ecbuild ${PHYEX_ECBUILD_VERSION} QUIET)

if(NOT ecbuild_FOUND)
  # ecbuild is not available yet: plain CMake must be used to fetch it.
  # Only the sources are needed (SOURCE_SUBDIR points to a non-existing
  # directory so that FetchContent does not call add_subdirectory)
  message(STATUS "ecbuild ${PHYEX_ECBUILD_VERSION} not found, fetching it from ${PHYEX_ECBUILD_REPO}")
  include(FetchContent)
  FetchContent_Declare(ecbuild
    GIT_REPOSITORY ${PHYEX_ECBUILD_REPO}
    GIT_TAG        ${PHYEX_ECBUILD_VERSION}
    GIT_SHALLOW    TRUE
    SOURCE_SUBDIR  _no_add_subdirectory_
  )
  FetchContent_MakeAvailable(ecbuild)
  find_package(ecbuild ${PHYEX_ECBUILD_VERSION} REQUIRED
               HINTS ${ecbuild_SOURCE_DIR} NO_DEFAULT_PATH)
endif()
message(STATUS "Using ecbuild ${ecbuild_VERSION} from ${ecbuild_DIR}")

#-----------------------------------------------------------------------
# Bundle project (ecbuild is now available)
#-----------------------------------------------------------------------

project(phyex_bundle VERSION 1.0 LANGUAGES Fortran)

include(ecbuild_bundle)
ecbuild_bundle_initialize()

#-----------------------------------------------------------------------
# Step 2: fiat
#-----------------------------------------------------------------------

# Once fiat has been fetched, it must stay a subproject on subsequent
# configurations (otherwise its build tree would be found as an installation)
if(NOT PHYEX_FIAT_FROM_SOURCE)
  ecbuild_find_package(NAME fiat VERSION ${PHYEX_FIAT_VERSION} QUIET)
endif()

if(fiat_FOUND)
  ecbuild_info("Using installed fiat ${fiat_VERSION} from ${fiat_DIR}")
else()
  set(PHYEX_FIAT_FROM_SOURCE ON CACHE INTERNAL "fiat is built as a subproject")
  ecbuild_info("fiat ${PHYEX_FIAT_VERSION} not found or already fetched, using sources from ${PHYEX_FIAT_REPO}")
  # Cloned in the build directory to keep the source tree clean
  set(PHYEX_FIAT_SOURCE_DIR ${CMAKE_BINARY_DIR}/_deps/fiat-src)
  ecbuild_git(PROJECT fiat
              DIR ${PHYEX_FIAT_SOURCE_DIR}
              URL ${PHYEX_FIAT_REPO}
              TAG ${PHYEX_FIAT_VERSION})
  add_subdirectory(${PHYEX_FIAT_SOURCE_DIR} ${CMAKE_BINARY_DIR}/fiat)
endif()

#-----------------------------------------------------------------------
# Step 3: python tools (pyphyextools), needed to build the python wheel
#-----------------------------------------------------------------------

ecbuild_add_option(FEATURE PYTHON_WHEEL
                   DESCRIPTION "Build the PHYEX python wheel (installs the python tools in the build directory)"
                   DEFAULT OFF)

if(HAVE_PYTHON_WHEEL)
  ecbuild_find_python(VERSION 3.8 NO_LIBS REQUIRED)

  # The tools are installed, with their dependencies, in a python virtual
  # environment inside the build directory so that each build can use its own
  # versions. They are installed at configure time because they are needed to
  # generate source code.
  # The pyphyextools package itself must be installed in editable mode (it
  # uses files from the repository): its source code is used directly from the
  # tools directory, a modification of the tools does not need a
  # re-installation.
  set(PHYEX_TOOLS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/tools)
  set(PHYEX_TOOLS_VENV ${CMAKE_BINARY_DIR}/pyphyextools)
  set(PHYEX_TOOLS_BINDIR ${PHYEX_TOOLS_VENV}/bin)
  set(PHYEX_TOOLS_PYTHON ${PHYEX_TOOLS_BINDIR}/python)

  # The virtual environment must not be polluted by a PYTHONPATH (set, for
  # example, by pip when the build is driven by pip in an isolated environment)
  unset(ENV{PYTHONPATH})

  # A modification of the dependencies triggers a new configuration and a
  # re-installation
  set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
               ${PHYEX_TOOLS_SOURCE_DIR}/pyproject.toml)
  file(SHA256 ${PHYEX_TOOLS_SOURCE_DIR}/pyproject.toml _deps_hash)

  if(NOT EXISTS ${PHYEX_TOOLS_BINDIR}/phyex-prep_code OR
     NOT _deps_hash STREQUAL PHYEX_TOOLS_DEPS_HASH)
    ecbuild_info("Installing python tools and their dependencies in ${PHYEX_TOOLS_VENV}")
    file(REMOVE_RECURSE ${PHYEX_TOOLS_VENV})
    execute_process(COMMAND ${Python_EXECUTABLE} -m venv ${PHYEX_TOOLS_VENV}
                    RESULT_VARIABLE _result)
    if(NOT _result EQUAL 0)
      ecbuild_critical("Creation of the python virtual environment ${PHYEX_TOOLS_VENV} failed")
    endif()
    execute_process(COMMAND ${PHYEX_TOOLS_PYTHON} -m pip install
                            --editable ${PHYEX_TOOLS_SOURCE_DIR}
                    RESULT_VARIABLE _result)
    if(NOT _result EQUAL 0 OR NOT EXISTS ${PHYEX_TOOLS_BINDIR}/phyex-prep_code)
      file(REMOVE_RECURSE ${PHYEX_TOOLS_VENV})
      ecbuild_critical("Installation of the python tools failed")
    endif()
    set(PHYEX_TOOLS_DEPS_HASH ${_deps_hash} CACHE INTERNAL "Hash of tools/pyproject.toml")
  else()
    ecbuild_info("Python tools are up to date in ${PHYEX_TOOLS_VENV}")
  endif()

  # Scripts installed in the virtual environment can be called directly
  # (e.g. ${PHYEX_TOOLS_BINDIR}/phyex-prep_code), but they also call other
  # tools by name (e.g. pyfortool): the virtual environment must come first
  # in PATH.
  #  - at configure time, PATH is set in the cmake environment
  #  - at build time, prefix the commands with ${PHYEX_TOOLS_ENV}, e.g.
  #    add_custom_command(... COMMAND ${PHYEX_TOOLS_ENV} phyex-prep_code ...)
  set(_path "${PHYEX_TOOLS_BINDIR}:$ENV{PATH}")
  set(PHYEX_TOOLS_ENV ${CMAKE_COMMAND} -E env --unset=PYTHONPATH "PATH=${_path}")
  set(ENV{PATH} "${_path}")
  ecbuild_info("Python tools available in ${PHYEX_TOOLS_BINDIR}")

  #---------------------------------------------------------------------
  # Python binding generation
  #---------------------------------------------------------------------

  # Fortran routines to bind: pairs of (source file relative to src/phyex, scope)
  set(PYPHYEX_BINDINGS
      micro/ice_adjust.F90          sub:ICE_ADJUST
      micro/rain_ice.F90            sub:RAIN_ICE
      micro/mode_ice4_sedimentation.F90
                                    module:MODE_ICE4_SEDIMENTATION/sub:ICE4_SEDIMENTATION
      micro/rain_ice_old.F90        sub:RAIN_ICE_OLD
      turb/shallow_mf.F90           sub:SHALLOW_MF
      turb/turb.F90                 sub:TURB
      aux/ini_phyex.F90             sub:INI_PHYEX
      micro/lima_adjust_split.F90   sub:LIMA_ADJUST_SPLIT
      micro/lima.F90                sub:LIMA
  )

  # Generated files:
  #  - FORTRAN interface, compiled in the PHYEX libraries
  #  - python module, next to the executables, loading the shared library
  #    (paths given relative to the python module, first found is used):
  #      - in the wheel, the module is pyphyex/__init__.py and the libraries
  #        are in pyphyex/lib (see install rules below)
  #      - in the build tree, the module is bin/pyphyex.py and the libraries
  #        are in lib (ecbuild puts executables in ${CMAKE_BINARY_DIR}/bin and
  #        libraries in ${CMAKE_BINARY_DIR}/lib)
  set(PYPHYEX_DIR ${CMAKE_BINARY_DIR}/pyphyex)
  set(PYPHYEX_F90 ${PYPHYEX_DIR}/pyphyex.F90)
  set(PYPHYEX_PY ${CMAKE_BINARY_DIR}/bin/pyphyex.py)
  set(PYPHYEX_WHEEL_LIBDIR lib)
  file(RELATIVE_PATH _libdir ${CMAKE_BINARY_DIR}/bin ${CMAKE_BINARY_DIR}/lib)
  set(_libso)
  foreach(_dir ${PYPHYEX_WHEEL_LIBDIR} ${_libdir})
    foreach(_precision dp sp)
      list(APPEND _libso
           ${_dir}/${CMAKE_SHARED_LIBRARY_PREFIX}phyex_${_precision}${CMAKE_SHARED_LIBRARY_SUFFIX})
    endforeach()
  endforeach()

  # phyex-pybinding appends to its output files: they are rebuilt from scratch
  set(_raw ${PYPHYEX_DIR}/pyphyex_raw.F90)
  file(MAKE_DIRECTORY ${PYPHYEX_DIR} ${CMAKE_BINARY_DIR}/bin)
  file(REMOVE ${_raw} ${PYPHYEX_PY})
  set(_inputs)
  list(LENGTH PYPHYEX_BINDINGS _n)
  math(EXPR _last "${_n} - 1")
  foreach(_i RANGE 0 ${_last} 2)
    math(EXPR _j "${_i} + 1")
    list(GET PYPHYEX_BINDINGS ${_i} _file)
    list(GET PYPHYEX_BINDINGS ${_j} _scope)
    set(_file ${CMAKE_CURRENT_SOURCE_DIR}/src/phyex/${_file})
    list(APPEND _inputs ${_file})
    execute_process(
      COMMAND ${PHYEX_TOOLS_BINDIR}/phyex-pybinding ${_file} ${_scope} ${_raw} ${PYPHYEX_PY} ${_libso}
      RESULT_VARIABLE _result)
    if(NOT _result EQUAL 0)
      ecbuild_critical("Python binding generation failed for ${_scope} (${_file})")
    endif()
  endforeach()

  # The generated routines are put in a module. The final file is only
  # written if its content changed, to avoid useless recompilations
  file(READ ${_raw} _content)
  file(WRITE ${_raw} "MODULE PYPHYEX\nCONTAINS\n${_content}\nEND MODULE PYPHYEX\n")
  configure_file(${_raw} ${PYPHYEX_F90} COPYONLY)

  # Any modification of the bound routines or of the generator triggers a
  # new configuration (and a new generation)
  set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
               ${_inputs} ${PHYEX_TOOLS_SOURCE_DIR}/src/pyphyextools/pybinding.py)
  ecbuild_info("Python binding generated in ${PYPHYEX_F90} and ${PYPHYEX_PY}")
endif()

#-----------------------------------------------------------------------
# Step 4: PHYEX
#-----------------------------------------------------------------------

add_subdirectory(src/phyex ${CMAKE_BINARY_DIR}/phyex)

if(HAVE_PYTHON_WHEEL)
  # The python binding is compiled in the PHYEX libraries
  set(_libs)
  foreach(_precision dp sp)
    if(TARGET phyex_${_precision})
      target_sources(phyex_${_precision} PRIVATE ${PYPHYEX_F90})
      list(APPEND _libs phyex_${_precision} parkind_${_precision})
    endif()
  endforeach()

  # Content of the python wheel, installed in the 'pyphyex' component only
  # (the wheel build installs only this component, see pyproject.toml):
  #   pyphyex/__init__.py : the generated python module
  #   pyphyex/lib/        : the shared libraries (their RPATH, $ORIGIN/../lib,
  #                         allows them to find each other)
  install(FILES ${PYPHYEX_PY} DESTINATION pyphyex RENAME __init__.py
          COMPONENT pyphyex)
  install(TARGETS ${_libs}
          LIBRARY DESTINATION pyphyex/${PYPHYEX_WHEEL_LIBDIR}
          COMPONENT pyphyex)
  get_target_property(_fiat_imported fiat IMPORTED)
  if(_fiat_imported)
    # fiat comes from an existing installation
    if(CMAKE_VERSION VERSION_LESS 3.21)
      ecbuild_critical("CMake >= 3.21 is needed to put an installed fiat in the python wheel")
    endif()
    install(IMPORTED_RUNTIME_ARTIFACTS fiat
            LIBRARY DESTINATION pyphyex/${PYPHYEX_WHEEL_LIBDIR}
            COMPONENT pyphyex)
  else()
    install(TARGETS fiat
            LIBRARY DESTINATION pyphyex/${PYPHYEX_WHEEL_LIBDIR}
            COMPONENT pyphyex)
  endif()
endif()

ecbuild_bundle_finalize()
