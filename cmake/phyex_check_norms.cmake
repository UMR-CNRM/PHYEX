#-----------------------------------------------------------------------
# Coding norms of the transformed PHYEX source code (included by the
# top-level CMakeLists.txt when ENABLE_CHECK_NORMS is ON and a transformation
# is active, i.e. when src/phyex defines the phyex_transformation target)
#
# The check is done as by the former build scripts (through phyex-prep_code),
# at build time, once all the files are transformed, in the phyex_check_norms
# target (part of the default build target; the compilation does not wait for
# it). It uses the coding_norms module of the python tools, directly from the
# tools directory: only pyfortool is needed (not all the python tools
# dependencies). This is why it is done here and not in src/phyex (the tools
# directory is not part of the PHYEX source directory).
#
# Uses, from src/phyex: the phyex_transformation target, the
# PHYEX_TRANSFORMED_DIR and PHYEX_TRANSFORMED_FILES directory properties and
# the PHYEX_PYFORTOOL cache variable.
#-----------------------------------------------------------------------

set(_phyex_dir ${CMAKE_CURRENT_SOURCE_DIR}/src/phyex)
get_property(_outdir DIRECTORY ${_phyex_dir} PROPERTY PHYEX_TRANSFORMED_DIR)
get_property(_outputs DIRECTORY ${_phyex_dir} PROPERTY PHYEX_TRANSFORMED_FILES)

# A python interpreter able to import pyfortool: the one of the pyfortool
# installation (virtual environment) or the one found by ecbuild
ecbuild_find_python(VERSION 3.8 NO_LIBS REQUIRED)
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
# directory (source tree). The stamp file is only created if the check succeeds
set(_norms_stamp ${CMAKE_BINARY_DIR}/phyex_check_norms.stamp)
add_custom_command(OUTPUT ${_norms_stamp}
                   COMMAND ${CMAKE_COMMAND} -E remove -f ${_norms_stamp}
                   COMMAND ${CMAKE_COMMAND} -E env
                           "PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/tools/src"
                           "PATH=${_pyft_bindir}:$ENV{PATH}"
                           PYTHONDONTWRITEBYTECODE=1
                           ${_norms_python} -m pyphyextools.testing.coding_norms -v ${_outdir}
                   COMMAND ${CMAKE_COMMAND} -E touch ${_norms_stamp}
                   DEPENDS ${_outputs}
                           ${CMAKE_CURRENT_SOURCE_DIR}/tools/src/pyphyextools/testing/coding_norms.py
                   COMMENT "Checking coding norms of the transformed source code"
                   VERBATIM)
# The transformed files are generated (in src/phyex) before the check
add_custom_target(phyex_check_norms ALL DEPENDS ${_norms_stamp})
add_dependencies(phyex_check_norms phyex_transformation)
ecbuild_info("Coding norms of the transformed source code checked at build time, using ${_norms_python}")
