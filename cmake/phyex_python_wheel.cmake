#-----------------------------------------------------------------------
# PHYEX python wheel (included by the top-level CMakeLists.txt when
# ENABLE_PYTHON_WHEEL is ON, after the PHYEX targets are defined)
#
#  - generation of the python binding (pyphyex), compiled in the PHYEX
#    libraries, by src/phyex/cmake/pybinding.py (pyfortool is needed, see
#    src/phyex/cmake/phyex_pyfortool.cmake)
#  - install rules of the 'pyphyex' component, used to build the wheel
#
# Variables defined for the rest of the configuration:
#   PYPHYEX_F90, PYPHYEX_PY
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
# Python binding generator (src/phyex/cmake/pybinding.py), only needs
# pyfortool (found or installed as for the source code transformation)
#-----------------------------------------------------------------------

include(${CMAKE_CURRENT_SOURCE_DIR}/src/phyex/cmake/phyex_pyfortool.cmake)
phyex_find_pyfortool()
set(PYPHYEX_BINDING_SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/src/phyex/cmake/pybinding.py)
ecbuild_info("Python binding generated with ${PYPHYEX_BINDING_SCRIPT}, using ${PHYEX_PYFORTOOL_PYTHON}")

#-----------------------------------------------------------------------
# Python binding generation
#-----------------------------------------------------------------------

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

# pybinding.py appends to its output files: they are rebuilt from scratch
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
    COMMAND ${PHYEX_PYFORTOOL_ENV} PYTHONDONTWRITEBYTECODE=1
            ${PHYEX_PYFORTOOL_PYTHON} ${PYPHYEX_BINDING_SCRIPT}
            ${_file} ${_scope} ${_raw} ${PYPHYEX_PY} ${_libso}
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
             ${_inputs} ${PYPHYEX_BINDING_SCRIPT})
ecbuild_info("Python binding generated in ${PYPHYEX_F90} and ${PYPHYEX_PY}")

#-----------------------------------------------------------------------
# Compilation of the binding and content of the wheel
#-----------------------------------------------------------------------

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
