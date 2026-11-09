#-----------------------------------------------------------------------
# Transformation of the PHYEX source code at build time
#
# When PHYEX_TRANSFORMATION is set (cache or environment variable, e.g.
# -DPHYEX_TRANSFORMATION=IAL), each Fortran source file given to
# phyex_transform_sources, and each included file (.h) of the PHYEX source
# directory, is transformed by pyfortool at build time (script
# phyex_transform.py, which defines the pyfortool arguments for each
# transformation). The targets compile the transformed files, kept in
# ${CMAKE_CURRENT_BINARY_DIR}/phyex_transformed.
# With an empty value (default) or NONE, nothing is changed.
#
# Additional pyfortool options can be given, at build time, with the PYFT_OPTS
# environment variable (a modification of PYFT_OPTS does not trigger a new
# transformation: remove the phyex_transformed directory).
#
# The transformed files are compiled (instead of using a compiler wrapper) so
# that the build system computes the compilation order from the transformed
# code (a transformation can add USE statements).
#
# Usage (in the PHYEX CMakeLists.txt):
#   phyex_transform_sources(<var> <sources>...) : <var> receives the list of
#                                                 sources to compile
#   phyex_transformation_target(<target>)        : for each target compiling
#                                                 transformed sources
#   phyex_transformation_finalize()              : once all the targets are
#                                                 defined
# When a transformation is active, phyex_transformation_finalize defines the
# phyex_transformation target (generation of all the transformed files) and
# the PHYEX_TRANSFORMED_DIR and PHYEX_TRANSFORMED_FILES directory properties
# (used, for example, to check the coding norms of the transformed code).
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
# Settings (default values can be given by environment variables)
#-----------------------------------------------------------------------

foreach(_setting
        "PHYEX_TRANSFORMATION;;Transformation applied to the PHYEX source code at build time (empty or NONE: no transformation)"
        "PHYEX_PYFORTOOL_SPEC;pyfortool>=0.2.24;pip requirement used to install pyfortool if it is not found")
  list(GET _setting 0 _var)
  list(GET _setting 1 _default)
  list(GET _setting 2 _doc)
  if(DEFINED ENV{${_var}})
    set(_default "$ENV{${_var}}")
  endif()
  set(${_var} "${_default}" CACHE STRING "${_doc}")
endforeach()

set(PHYEX_TRANSFORM_SCRIPT ${CMAKE_CURRENT_LIST_DIR}/phyex_transform.py)
set(PHYEX_TRANSFORM_COMMAND)
set(PHYEX_TRANSFORMED_SRCDIR ${CMAKE_CURRENT_SOURCE_DIR})
set(PHYEX_TRANSFORMED_DIR ${CMAKE_CURRENT_BINARY_DIR}/phyex_transformed)
set(PHYEX_TRANSFORMED_FILES)
set(PHYEX_TRANSFORMED_INCLUDEDIRS)

if(PHYEX_TRANSFORMATION AND NOT PHYEX_TRANSFORMATION STREQUAL "NONE")

  #---------------------------------------------------------------------
  # Transformation command (and pyfortool)
  #---------------------------------------------------------------------

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

  # pyfortool: found in PATH, otherwise installed in a python virtual
  # environment inside the build directory
  set(_pyft_venv ${CMAKE_CURRENT_BINARY_DIR}/pyfortool)
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
  # pyfortool installed in the build directory must not be polluted by a
  # PYTHONPATH
  string(FIND "${PHYEX_PYFORTOOL}" "${CMAKE_BINARY_DIR}/" _pos)
  if(_pos EQUAL 0)
    list(APPEND PHYEX_TRANSFORM_COMMAND --unset-pythonpath)
  endif()

  # Changing the transformation triggers a new transformation
  set(PHYEX_TRANSFORM_STAMP ${CMAKE_CURRENT_BINARY_DIR}/phyex_transformation.txt)
  file(WRITE ${PHYEX_TRANSFORM_STAMP}.tmp "${PHYEX_TRANSFORMATION}\n")
  configure_file(${PHYEX_TRANSFORM_STAMP}.tmp ${PHYEX_TRANSFORM_STAMP} COPYONLY)
endif()

#-----------------------------------------------------------------------
# Internal: transformation of one file (relative path _rel) into _out
#-----------------------------------------------------------------------

macro(_phyex_transform_file _src _rel _out)
  add_custom_command(OUTPUT ${_out}
                     COMMAND ${PHYEX_TRANSFORM_COMMAND} ${_src} ${_out}
                     DEPENDS ${_src} ${PHYEX_TRANSFORM_SCRIPT} ${PHYEX_TRANSFORM_STAMP}
                     COMMENT "Transforming (${PHYEX_TRANSFORMATION}) ${_rel}"
                     VERBATIM)
  list(APPEND PHYEX_TRANSFORMED_FILES ${_out})
endmacro()

#-----------------------------------------------------------------------
# Included files (.h) are transformed too. Their transformed version is found
# first: it is in the same directory as the transformed source file including
# it, or in a directory put first in the include path
#-----------------------------------------------------------------------

if(PHYEX_TRANSFORM_COMMAND)
  file(GLOB_RECURSE _headers ${PHYEX_TRANSFORMED_SRCDIR}/*.h)
  foreach(_src ${_headers})
    file(RELATIVE_PATH _rel ${PHYEX_TRANSFORMED_SRCDIR} ${_src})
    set(_out ${PHYEX_TRANSFORMED_DIR}/${_rel})
    _phyex_transform_file(${_src} ${_rel} ${_out})
    get_filename_component(_dir ${_out} DIRECTORY)
    list(APPEND PHYEX_TRANSFORMED_INCLUDEDIRS ${_dir})
  endforeach()
  list(REMOVE_DUPLICATES PHYEX_TRANSFORMED_INCLUDEDIRS)
endif()

#-----------------------------------------------------------------------
# phyex_transform_sources(<var> <sources>...)
# <var> receives the sources to compile: the Fortran files of the PHYEX
# source directory are replaced by their transformed version (one
# transformation per file, shared by all the targets: dp, sp...)
#-----------------------------------------------------------------------

function(phyex_transform_sources _var)
  set(_sources)
  foreach(_src ${ARGN})
    if(PHYEX_TRANSFORM_COMMAND AND _src MATCHES "\\.[fF](90|95|03|08)?$")
      get_filename_component(_src ${_src} ABSOLUTE)
      file(RELATIVE_PATH _rel ${PHYEX_TRANSFORMED_SRCDIR} ${_src})
      if(NOT _rel MATCHES "^\\.\\./")
        set(_out ${PHYEX_TRANSFORMED_DIR}/${_rel})
        if(NOT _out IN_LIST PHYEX_TRANSFORMED_FILES)
          _phyex_transform_file(${_src} ${_rel} ${_out})
        endif()
        set(_src ${_out})
      endif()
    endif()
    list(APPEND _sources ${_src})
  endforeach()
  set(${_var} ${_sources} PARENT_SCOPE)
  set(PHYEX_TRANSFORMED_FILES ${PHYEX_TRANSFORMED_FILES} PARENT_SCOPE)
endfunction()

#-----------------------------------------------------------------------
# phyex_transformation_target(<target>)
# The transformed headers are found first, and the transformed files are
# generated before the target is compiled
#-----------------------------------------------------------------------

function(phyex_transformation_target _target)
  if(PHYEX_TRANSFORM_COMMAND)
    target_include_directories(${_target} BEFORE PRIVATE ${PHYEX_TRANSFORMED_INCLUDEDIRS})
    add_dependencies(${_target} phyex_transformation)
  endif()
endfunction()

#-----------------------------------------------------------------------
# phyex_transformation_finalize()
# The transformed files are generated by a dedicated target: several targets
# (dp, sp...) use the same transformed files, they must not generate them
# concurrently
#-----------------------------------------------------------------------

function(phyex_transformation_finalize)
  if(PHYEX_TRANSFORM_COMMAND)
    add_custom_target(phyex_transformation DEPENDS ${PHYEX_TRANSFORMED_FILES})
    set_property(DIRECTORY PROPERTY PHYEX_TRANSFORMED_DIR ${PHYEX_TRANSFORMED_DIR})
    set_property(DIRECTORY PROPERTY PHYEX_TRANSFORMED_FILES ${PHYEX_TRANSFORMED_FILES})
    ecbuild_info("${PHYEX_TRANSFORMATION} transformation of the PHYEX source code in ${PHYEX_TRANSFORMED_DIR}")
  endif()
endfunction()
