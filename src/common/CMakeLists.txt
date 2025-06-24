cmake_minimum_required(VERSION 3.10)
project(PHYEX)
enable_language(Fortran)
# set(CMAKE_VERBOSE_MAKEFILE on)

#################################################################################
# Progs compilation
set(PHYEX_BUILD_PROGS AUTO CACHE STRING "Build the program executables AUTO|ON|OFF")
set_property(CACHE PHYEX_BUILD_PROGS PROPERTY STRINGS AUTO ON OFF)
if(${PHYEX_BUILD_PROGS} STREQUAL "AUTO")
  if(EXISTS ${CMAKE_SOURCE_DIR}/progs)
    set(PHYEX_BUILD_PROGS ON)
  else()
    set(PHYEX_BUILD_PROGS OFF)
  endif()
endif()

#################################################################################
# Library building
set(ENABLE_SHARED_LIB BOTH CACHE STRING "Build shared lib ON|OFF|BOTH|NONE")
set_property(CACHE ENABLE_SHARED_LIB PROPERTY STRINGS ON OFF BOTH NONE)
if(${ENABLE_SHARED_LIB} STREQUAL "ON" OR ${ENABLE_SHARED_LIB} STREQUAL "BOTH")
  set(BUILD_SHARED ON)
else()
  set(BUILD_SHARED OFF)
endif()
if(${ENABLE_SHARED_LIB} STREQUAL "OFF" OR ${ENABLE_SHARED_LIB} STREQUAL "BOTH")
  set(BUILD_STATIC ON)
else()
  set(BUILD_STATIC OFF)
endif()

#################################################################################
# Fiat
set(fiat_ROOT "" CACHE FILEPATH "Path to fiat installation directory")
if(NOT EXISTS ${fiat_ROOT})
  message(FATAL_ERROR "fiat_ROOT is not set correctly!")
endif()
include_directories(${fiat_ROOT}/module/fiat)
include_directories(${fiat_ROOT}/include/fiat)

#################################################################################
# Default real type
set(REALDEFAULTSIZE 0)
foreach(ONESIZE 4 8)
  set(TEST_FILE ${CMAKE_SOURCE_DIR}/test_fortan.F90)
  file(WRITE ${TEST_FILE}
  "
  PROGRAM test_size
    REAL(KIND=${ONESIZE}), TARGET :: X
    REAL, POINTER :: PX
    PX => X
  END PROGRAM
  ")
  try_compile(RES ${CMAKE_CURRENT_BINARY_DIR} ${TEST_FILE})
  if(RES)
    message(STATUS "Default real size is ${ONESIZE}")
    set(REALDEFAULTSIZE ${ONESIZE})
  endif()
endforeach()
if(${REALDEFAULTSIZE} EQUAL 4)
  include_directories(${fiat_ROOT}/module/parkind_sp)
elseif(${REALDEFAULTSIZE} EQUAL 8)
  include_directories(${fiat_ROOT}/module/parkind_dp)
else()
  message(FATAL_ERROR "Unable to determine default real size")
endif()
file(REMOVE ${TEST_FILE})

#################################################################################
# Group module files to allow their installation
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Mod)

#################################################################################
# Build and install shared and static libraries (with module and include files)

# Source files: directories
set(LIB_DIRECTORY_LIST aux turb micro)

# Source files: entry points
set(LIB_ENTRYPOINTS ${CMAKE_SOURCE_DIR}/aux/ini_phyex.F90
                    ${CMAKE_SOURCE_DIR}/micro/rain_ice.F90
                    ${CMAKE_SOURCE_DIR}/turb/shallow_mf.F90
                    ${CMAKE_SOURCE_DIR}/turb/turb.F90
                    ${CMAKE_SOURCE_DIR}/micro/ice_adjust.F90
                    ${CMAKE_SOURCE_DIR}/micro/lima.F90
                    ${CMAKE_SOURCE_DIR}/micro/lima_adjust_split.F90)
if(EXISTS ${CMAKE_SOURCE_DIR}/pyphyex.F90)
   list(APPEND LIB_ENTRYPOINTS ${CMAKE_SOURCE_DIR}/pyphyex.F90)
endif()

# Source files: all other source files (and specific list for header files)
set(LIB_FILES "")
set(HEADER_FILES "")
foreach(D ${LIB_DIRECTORY_LIST})
  file(GLOB DFILES ${D}/*)
  list(APPEND LIB_FILES ${DFILES})
  file(GLOB DFILES ${D}/*.h)
  list(APPEND HEADER_FILES ${DFILES})
endforeach()
list(REMOVE_ITEM LIB_FILES ${LIB_ENTRYPOINTS})

# This is the objlib target; this intermediary library is used
# by both the static and the shared libs. This way, the sources are
# compiled only once.
# Because the same objects are used for static and shared, they must be
# compiled with the PIC option.
add_library(objlib OBJECT ${LIB_FILES})
target_include_directories(objlib PRIVATE ${LIB_DIRECTORY_LIST})
set_property(TARGET objlib PROPERTY POSITION_INDEPENDENT_CODE 1)

add_library(entrypointlib OBJECT ${LIB_ENTRYPOINTS})
add_dependencies(entrypointlib objlib)
target_include_directories(entrypointlib PRIVATE ${LIB_DIRECTORY_LIST})
set_property(TARGET entrypointlib PROPERTY POSITION_INDEPENDENT_CODE 1)

# Static library
if(${BUILD_STATIC} OR ${PHYEX_BUILD_PROGS})
  add_library(phyexStatic STATIC $<TARGET_OBJECTS:entrypointlib> $<TARGET_OBJECTS:objlib>)
  target_include_directories(phyexStatic PRIVATE ${LIB_DIRECTORY_LIST})
endif()
if(${BUILD_STATIC})
  set_target_properties(phyexStatic PROPERTIES OUTPUT_NAME phyex)
  install(TARGETS phyexStatic LIBRARY DESTINATION lib)
endif()

# Shared library
if(${BUILD_SHARED})
  add_library(phyexShared SHARED $<TARGET_OBJECTS:entrypointlib>)
  set_property(TARGET phyexShared PROPERTY POSITION_INDEPENDENT_CODE 1)
  target_link_directories(phyexShared PRIVATE ${fiat_ROOT}/lib)
  target_link_directories(phyexShared PRIVATE ${fiat_ROOT}/lib64)
  target_link_libraries(phyexShared PRIVATE fiat)
  target_link_libraries(phyexShared PRIVATE mpi_serial)
  target_link_libraries(phyexShared PRIVATE $<TARGET_OBJECTS:objlib>)
  target_include_directories(phyexShared PRIVATE ${LIB_DIRECTORY_LIST})
  set_target_properties(phyexShared PROPERTIES OUTPUT_NAME phyex)
  install(TARGETS phyexShared LIBRARY DESTINATION lib)
endif()

# Installation
if(${BUILD_SHARED} OR ${BUILD_STATIC})
  install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/
          DESTINATION module/phyex)
  install(FILES ${HEADER_FILES}
          DESTINATION include/phyex)
endif()

#################################################################################
# Build and install progs

if(${PHYEX_BUILD_PROGS})
  # Source file directories
  set(PROGS_DIRECTORY_LIST support progs)

  # All source files
  set(PROG_FILES "")
  foreach(D ${PROGS_DIRECTORY_LIST})
    file(GLOB DFILES ${D}/*)
    list(APPEND PROG_FILES ${DFILES})
  endforeach()

  # Entry points
  file(GLOB PROGS_ENTRYPOINTS progs/main_*.F90)

  # We must excludes from PROG_FILES objects used as entry points
  list(REMOVE_ITEM PROG_FILES ${PROGS_ENTRYPOINTS})

  # Buid a lib with common files to all executables
  add_library(progsobj OBJECT ${PROG_FILES})
  if(DEFINED FPP_FLAGS_TESTPROGS)
    target_compile_definitions(progsobj PUBLIC ${FPP_FLAGS_TESTPROGS})
  endif()
  
  foreach(EP ${PROGS_ENTRYPOINTS})
    get_filename_component(NAME ${EP} NAME_WLE)
    add_executable(${NAME})
    target_sources(${NAME} PRIVATE ${EP})
    target_link_directories(${NAME} PRIVATE ${fiat_ROOT}/lib)
    target_link_directories(${NAME} PRIVATE ${fiat_ROOT}/lib64)
    target_link_libraries(${NAME} PRIVATE phyexStatic progsobj fiat mpi_serial)
    if(DEFINED FPP_FLAGS_TESTPROGS)
      target_compile_definitions(${NAME} PUBLIC ${FPP_FLAGS_TESTPROGS})
    endif()
    set_target_properties(${NAME} PROPERTIES OUTPUT_NAME ${NAME}.exe)
    install(TARGETS ${NAME} RUNTIME DESTINATION bin)
  endforeach()
endif()
