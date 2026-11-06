#-----------------------------------------------------------------------
# PHYEX top-level CMakeLists.txt
#
# Builds the PHYEX library (src/phyex) together with its dependencies:
#   - ecbuild : used if found, otherwise downloaded at configure time
#   - fiat    : used if found, otherwise downloaded and built as a subproject
#
# Dependency versions and repositories can be overridden with environment
# variables (read at first configure) or with -D<VAR>=... on the command line:
#   PHYEX_ECBUILD_VERSION, PHYEX_ECBUILD_REPO, PHYEX_FIAT_VERSION, PHYEX_FIAT_REPO,
#   PHYEX_PYFORTOOL_SPEC
#
# The PHYEX source code can be transformed at build time by pyfortool with
# -DPHYEX_TRANSFORMATION=<name> (or the environment variable of the same
# name), see step 4. Example: -DPHYEX_TRANSFORMATION=IAL
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
phyex_dependency_setting(PHYEX_PYFORTOOL_SPEC "pyfortool>=0.2.24"
                         "pip requirement used to install pyfortool if it is not found")

# Transformation of the PHYEX source code at build time (see step 4)
phyex_dependency_setting(PHYEX_TRANSFORMATION ""
                         "Transformation applied to the PHYEX source code at build time (empty or NONE: no transformation)")

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
# Step 4: source code transformation at build time (optional)
#-----------------------------------------------------------------------

# When PHYEX_TRANSFORMATION is set (e.g. -DPHYEX_TRANSFORMATION=IAL), each
# Fortran source file of the PHYEX targets, and each included file (.h) of
# src/phyex, is transformed by pyfortool at build
# time (script cmake/phyex_transform.py, which defines the pyfortool arguments
# for each transformation), and the PHYEX targets compile the transformed
# files, kept in ${CMAKE_BINARY_DIR}/phyex_transformed.
# Additional pyfortool options can be given, at build time, with the PYFT_OPTS
# environment variable (a modification of PYFT_OPTS does not trigger a new
# transformation: remove the phyex_transformed directory).
# The transformed files are compiled (instead of using a compiler wrapper) so
# that the build system computes the compilation order from the transformed
# code (a transformation can add USE statements).
set(PHYEX_TRANSFORM_SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/cmake/phyex_transform.py)
set(PHYEX_TRANSFORM_COMMAND)
if(PHYEX_TRANSFORMATION AND NOT PHYEX_TRANSFORMATION STREQUAL "NONE")
  ecbuild_find_python(VERSION 3.8 NO_LIBS REQUIRED)
  # -I: the script only needs the standard library, and must not be
  # influenced by PYTHONPATH
  set(_script ${Python_EXECUTABLE} -I ${PHYEX_TRANSFORM_SCRIPT})

  # Available transformations are defined in the script
  execute_process(COMMAND ${_script} --list
                  OUTPUT_VARIABLE _transformations OUTPUT_STRIP_TRAILING_WHITESPACE
                  RESULT_VARIABLE _result)
  if(NOT _result EQUAL 0)
    ecbuild_critical("Unable to execute ${PHYEX_TRANSFORM_SCRIPT}")
  endif()
  set_property(CACHE PHYEX_TRANSFORMATION PROPERTY STRINGS "NONE;${_transformations}")
  if(NOT PHYEX_TRANSFORMATION IN_LIST _transformations)
    string(REPLACE ";" ", " _values "NONE;${_transformations}")
    ecbuild_critical("Unknown transformation PHYEX_TRANSFORMATION=${PHYEX_TRANSFORMATION}, available values are: ${_values}")
  endif()

  # pyfortool: found in PATH (the python tools virtual environment comes first
  # when ENABLE_PYTHON_WHEEL is ON), otherwise installed in a python virtual
  # environment inside the build directory
  set(_pyft_venv ${CMAKE_BINARY_DIR}/pyfortool)
  if(PHYEX_PYFORTOOL AND NOT EXISTS ${PHYEX_PYFORTOOL})
    unset(PHYEX_PYFORTOOL CACHE)
  endif()
  find_program(PHYEX_PYFORTOOL pyfortool DOC "pyfortool executable")
  if(NOT PHYEX_PYFORTOOL)
    ecbuild_info("pyfortool not found, installing ${PHYEX_PYFORTOOL_SPEC} in ${_pyft_venv}")
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
  endif()
  ecbuild_info("Source code transformation ${PHYEX_TRANSFORMATION} applied at build time, using ${PHYEX_PYFORTOOL}")

  set(PHYEX_TRANSFORM_COMMAND ${_script}
      --transformation ${PHYEX_TRANSFORMATION}
      --pyfortool ${PHYEX_PYFORTOOL})
  # pyfortool installed in the build directory (by this file) must not be
  # polluted by a PYTHONPATH
  string(FIND "${PHYEX_PYFORTOOL}" "${CMAKE_BINARY_DIR}/" _pos)
  if(_pos EQUAL 0)
    list(APPEND PHYEX_TRANSFORM_COMMAND --unset-pythonpath)
  endif()
endif()

# Coding norms of the transformed source code (as done by the former build
# scripts through phyex-prep_code). The check is executed at build time, once
# all the files are transformed, in the phyex_check_norms target (part of the
# default build target; the compilation does not wait for it).
# It uses the coding_norms module of the python tools, directly from the tools
# directory: only pyfortool is needed (not all the python tools dependencies)
ecbuild_add_option(FEATURE CHECK_NORMS
                   DESCRIPTION "Check the coding norms of the transformed source code (PHYEX_TRANSFORMATION)"
                   DEFAULT ON)
set(PHYEX_CHECK_NORMS_COMMAND)
if(HAVE_CHECK_NORMS AND PHYEX_TRANSFORM_COMMAND)
  # A python interpreter able to import pyfortool: the one of the pyfortool
  # installation (virtual environment) or the one found by ecbuild
  get_filename_component(_pyft_bindir ${PHYEX_PYFORTOOL} DIRECTORY)
  set(_norms_python)
  foreach(_python ${_pyft_bindir}/python ${_pyft_bindir}/python3 ${Python_EXECUTABLE})
    if(EXISTS ${_python})
      execute_process(COMMAND ${CMAKE_COMMAND} -E env --unset=PYTHONPATH
                              ${_python} -c "import pyfortool"
                      RESULT_VARIABLE _result OUTPUT_QUIET ERROR_QUIET)
      if(_result EQUAL 0)
        set(_norms_python ${_python})
        break()
      endif()
    endif()
  endforeach()
  if(NOT _norms_python)
    ecbuild_critical("No python interpreter able to import pyfortool found for the coding norms check (use -DENABLE_CHECK_NORMS=OFF to disable the check)")
  endif()
  if(NOT EXISTS ${_pyft_bindir}/pyfortool_parallel)
    ecbuild_critical("pyfortool_parallel not found in ${_pyft_bindir}, needed for the coding norms check (use -DENABLE_CHECK_NORMS=OFF to disable the check)")
  endif()
  # pyfortool_parallel is called by name; no bytecode written in the tools
  # directory (source tree)
  set(PHYEX_CHECK_NORMS_COMMAND ${CMAKE_COMMAND} -E env
      "PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/tools/src"
      "PATH=${_pyft_bindir}:$ENV{PATH}"
      PYTHONDONTWRITEBYTECODE=1
      ${_norms_python} -m pyphyextools.testing.coding_norms -v)
  ecbuild_info("Coding norms of the transformed source code checked at build time, using ${_norms_python}")
elseif(HAVE_CHECK_NORMS)
  ecbuild_debug("No source code transformation: CHECK_NORMS has no effect")
endif()

#-----------------------------------------------------------------------
# Step 5: PHYEX
#-----------------------------------------------------------------------

add_subdirectory(src/phyex ${CMAKE_BINARY_DIR}/phyex)

if(PHYEX_TRANSFORM_COMMAND)
  # Source transformation (step 4): the Fortran source files of the PHYEX
  # targets are replaced by their transformed version
  set(_srcdir ${CMAKE_CURRENT_SOURCE_DIR}/src/phyex)
  set(_outdir ${CMAKE_BINARY_DIR}/phyex_transformed)
  # Changing the transformation triggers a new transformation
  set(_stamp ${CMAKE_BINARY_DIR}/phyex_transformation.txt)
  file(WRITE ${_stamp}.tmp "${PHYEX_TRANSFORMATION}\n")
  configure_file(${_stamp}.tmp ${_stamp} COPYONLY)

  set(_outputs)
  set(_transformed_targets)

  # Included files (.h) are transformed too. Their transformed version is
  # found first: it is in the same directory as the transformed source file
  # including it, or in a directory put first in the include path
  file(GLOB_RECURSE _headers ${_srcdir}/*.h)
  set(_header_dirs)
  foreach(_src ${_headers})
    file(RELATIVE_PATH _rel ${_srcdir} ${_src})
    set(_out ${_outdir}/${_rel})
    add_custom_command(OUTPUT ${_out}
                       COMMAND ${PHYEX_TRANSFORM_COMMAND} ${_src} ${_out}
                       DEPENDS ${_src} ${PHYEX_TRANSFORM_SCRIPT} ${_stamp}
                       COMMENT "Transforming (${PHYEX_TRANSFORMATION}) ${_rel}"
                       VERBATIM)
    list(APPEND _outputs ${_out})
    get_filename_component(_dir ${_out} DIRECTORY)
    list(APPEND _header_dirs ${_dir})
  endforeach()
  list(REMOVE_DUPLICATES _header_dirs)

  get_property(_targets DIRECTORY ${_srcdir} PROPERTY BUILDSYSTEM_TARGETS)
  foreach(_target ${_targets})
    get_target_property(_type ${_target} TYPE)
    if(NOT _type MATCHES "^(STATIC_LIBRARY|SHARED_LIBRARY|MODULE_LIBRARY|OBJECT_LIBRARY|EXECUTABLE)$")
      continue()
    endif()
    get_target_property(_sources ${_target} SOURCES)
    get_target_property(_target_srcdir ${_target} SOURCE_DIR)
    set(_new_sources)
    set(_incdirs)
    set(_changed FALSE)
    foreach(_src ${_sources})
      if(_src MATCHES "\\$<")
        # Generator expression (e.g. $<TARGET_OBJECTS:...>): kept unchanged
        list(APPEND _new_sources ${_src})
        continue()
      endif()
      if(NOT IS_ABSOLUTE ${_src})
        set(_src ${_target_srcdir}/${_src})
      endif()
      file(RELATIVE_PATH _rel ${_srcdir} ${_src})
      if(_src MATCHES "\\.[fF](90|95|03|08)?$" AND NOT _rel MATCHES "^\\.\\./")
        set(_out ${_outdir}/${_rel})
        if(NOT _out IN_LIST _outputs)
          # One command per file, shared by all the targets (dp, sp...)
          add_custom_command(OUTPUT ${_out}
                             COMMAND ${PHYEX_TRANSFORM_COMMAND} ${_src} ${_out}
                             DEPENDS ${_src} ${PHYEX_TRANSFORM_SCRIPT} ${_stamp}
                             COMMENT "Transforming (${PHYEX_TRANSFORMATION}) ${_rel}"
                             VERBATIM)
          list(APPEND _outputs ${_out})
        endif()
        list(APPEND _new_sources ${_out})
        # Files included by the original source code are also searched in its
        # directory (after the transformed headers)
        get_filename_component(_dir ${_src} DIRECTORY)
        list(APPEND _incdirs ${_dir})
        set(_changed TRUE)
      else()
        list(APPEND _new_sources ${_src})
      endif()
    endforeach()
    if(_changed)
      set_property(TARGET ${_target} PROPERTY SOURCES ${_new_sources})
      list(REMOVE_DUPLICATES _incdirs)
      target_include_directories(${_target} BEFORE PRIVATE ${_header_dirs} ${_incdirs})
      list(APPEND _transformed_targets ${_target})
    endif()
  endforeach()

  # The transformed files are generated by a dedicated target, built before
  # the PHYEX targets (which are defined in another directory)
  set_source_files_properties(${_outputs} DIRECTORY ${_srcdir} PROPERTIES GENERATED TRUE)
  add_custom_target(phyex_transformation DEPENDS ${_outputs})

  # Coding norms of the transformed source code (see step 4). The stamp file
  # is only created if the check succeeds
  if(PHYEX_CHECK_NORMS_COMMAND)
    set(_norms_stamp ${CMAKE_BINARY_DIR}/phyex_check_norms.stamp)
    add_custom_command(OUTPUT ${_norms_stamp}
                       COMMAND ${CMAKE_COMMAND} -E remove -f ${_norms_stamp}
                       COMMAND ${PHYEX_CHECK_NORMS_COMMAND} ${_outdir}
                       COMMAND ${CMAKE_COMMAND} -E touch ${_norms_stamp}
                       DEPENDS ${_outputs}
                               ${CMAKE_CURRENT_SOURCE_DIR}/tools/src/pyphyextools/testing/coding_norms.py
                       COMMENT "Checking coding norms of the transformed source code"
                       VERBATIM)
    add_custom_target(phyex_check_norms ALL DEPENDS ${_norms_stamp})
    add_dependencies(phyex_check_norms phyex_transformation)
  endif()
  foreach(_target ${_transformed_targets})
    add_dependencies(${_target} phyex_transformation)
  endforeach()
  ecbuild_info("${PHYEX_TRANSFORMATION} transformation of the PHYEX source code in ${_outdir}")
endif()

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
