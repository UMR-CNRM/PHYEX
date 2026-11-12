cmake_minimum_required(VERSION 3.10)
find_package(ecbuild 3.11 REQUIRED HINTS ${CMAKE_SOURCE_DIR})

# PHYEX version: the VERSION file is the only place where it is stored (also
# used by the top-level CMakeLists.txt, the python wheel and the python tools)
if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/VERSION)
  message(FATAL_ERROR "${CMAKE_CURRENT_SOURCE_DIR}/VERSION (PHYEX version) not found")
endif()
file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/VERSION PHYEX_VERSION
     LIMIT_COUNT 1 REGEX "^[0-9]+(\\.[0-9]+)*$")
if(NOT PHYEX_VERSION)
  message(FATAL_ERROR "No valid version (x.y.z) found in ${CMAKE_CURRENT_SOURCE_DIR}/VERSION")
endif()
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/VERSION)
project(phyex VERSION ${PHYEX_VERSION} LANGUAGES Fortran)

include(phyex_macros)
include(phyex_transformation)
include(phyex_modi)

#-----------------------------------------------------------------------
# Options
#-----------------------------------------------------------------------

ecbuild_add_option(FEATURE DOUBLE_PRECISION
                   DESCRIPTION "Compile on double precision"
                   DEFAULT ON)

ecbuild_add_option(FEATURE SINGLE_PRECISION
                   DESCRIPTION "Compile on single precision"
                   DEFAULT OFF)

ecbuild_add_option(FEATURE IFS_FLAGS
                   DESCRIPTION "Compile with IFS flags"
                   DEFAULT OFF)

ecbuild_add_option(FEATURE PHYEX_BUILD_PROGS
                   DESCRIPTION "Build PHYEX program executables"
                   DEFAULT OFF)

ecbuild_add_option(FEATURE OMP
                   DESCRIPTION "Support for OpenMP shared memory parallelism"
                   REQUIRED_PACKAGES "OpenMP COMPONENTS Fortran")

ecbuild_add_option(FEATURE ACC
                   DESCRIPTION "Support for OpenACC"
                   REQUIRED_PACKAGES "OpenACC COMPONENTS Fortran"
                   DEFAULT OFF)

if(HAVE_DOUBLE_PRECISION)
  list(APPEND precision_list dp) 
endif()

if(HAVE_SINGLE_PRECISION)
  list(APPEND precision_list sp)
endif()

#-----------------------------------------------------------------------
# Dependencies 
#-----------------------------------------------------------------------

ecbuild_find_package(NAME fiat REQUIRED)

#-----------------------------------------------------------------------
# Source structure
#-----------------------------------------------------------------------

# Common sources
file(GLOB PHYEX_SOURCES
     ${CMAKE_CURRENT_SOURCE_DIR}/aux/*F90
     ${CMAKE_CURRENT_SOURCE_DIR}/turb/*F90
     ${CMAKE_CURRENT_SOURCE_DIR}/micro/*F90
     ${CMAKE_CURRENT_SOURCE_DIR}/conv/*F90
)

# Source for the python binding
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/pyphyex.F90)
  list(APPEND PHYEX_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/pyphyex.F90)
endif()

# Sources for the programs
if(HAVE_PHYEX_BUILD_PROGS)
  # Entry points
  file(GLOB PHYEXPROGS_ENTRYPOINTS progs/main_*.F90)

  # Source files
  file(GLOB PHYEXPROGS_SOURCES
       ${CMAKE_CURRENT_SOURCE_DIR}/support/*90
       ${CMAKE_CURRENT_SOURCE_DIR}/progs/*90
  ) 
  list(REMOVE_ITEM PHYEXPROGS_SOURCES ${PHYEXPROGS_ENTRYPOINTS})
  list(APPEND PHYEX_SOURCES ${PHYEXPROGS_SOURCES})
endif()

# Optional transformation of the source code (PHYEX_TRANSFORMATION, see
# cmake/phyex_transformation.cmake); without transformation, the list is unchanged
phyex_transform_sources(PHYEX_SOURCES ${PHYEX_SOURCES})

# Interface modules (MODI_*) generated at build time, from the transformed
# source code if any (see cmake/phyex_modi.cmake and cmake/modi_sources.txt)
phyex_generate_modi(PHYEX_SOURCES)

# Headers
set(GLOB PHYEX_HEADERS
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/phyex/turb>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/phyex/micro>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/phyex/conv>
)

# Include directories
set(PHYEX_INCLUDEDIRS
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/aux>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/micro>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/turb>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/conv>
)

#---------------
#  Loop on precisions
#---------------

foreach(PRECISION ${precision_list})
  #---------------
  #  Determine suitable compilation option to use for unqualified REALs
  #---------------

  if(${PRECISION} STREQUAL "dp")
    set(PHYEX_COMPILE_OPTIONS_REAL ${ECBUILD_Fortran_COMPILE_OPTIONS_REAL8})
  else()
    set(PHYEX_COMPILE_OPTIONS_REAL ${ECBUILD_Fortran_COMPILE_OPTIONS_REAL4})
  endif()
  
  #---------------
  #  Public library
  #---------------

  ecbuild_add_library(
    TARGET phyex_${PRECISION}
    SOURCES ${PHYEX_SOURCES}
    PUBLIC_LIBS
      fiat
      parkind_${PRECISION}
    PUBLIC_INCLUDES ${PHYEX_HEADERS}
    INSTALL_HEADERS ALL
    HEADER_DESTINATION include/${PRECISION}
    PUBLIC_INCLUDES ${PHYEX_INCLUDEDIRS}
  ) 
  phyex_transformation_target(phyex_${PRECISION})
  phyex_modi_target(phyex_${PRECISION})

  if(HAVE_OMP)
    target_link_libraries(phyex_${PRECISION} PRIVATE OpenMP::OpenMP_Fortran)
  endif()

  if(HAVE_ACC)
    target_link_libraries(phyex_${PRECISION} PRIVATE OpenACC::OpenACC_Fortran)
  endif()

  # Compilation option to use for unqualified REALs
  target_compile_options(phyex_${PRECISION} PRIVATE
    $<$<COMPILE_LANGUAGE:Fortran>:${PHYEX_COMPILE_OPTIONS_REAL}>)

  if(HAVE_IFS_FLAGS)
    target_compile_options(phyex_${PRECISION} PRIVATE
    ${ECBUILD_IFS_Fortran_FLAGS})
  endif()

  # Install Library 
  install(
    TARGETS phyex_${PRECISION}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )
  
  # Install Modules 
  fortran_modules_install(
    TARGET phyex_${PRECISION}
    MODULE_DIR ${CMAKE_CURRENT_BINARY_DIR}/modules_${PRECISION}
  	INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/modules/${PRECISION}
  )
  
  #---------------
  #  Progs
  #---------------
  
  if(HAVE_PHYEX_BUILD_PROGS)
    # Build all executables
    foreach(EP ${PHYEXPROGS_ENTRYPOINTS})
      get_filename_component(NAME ${EP} NAME_WLE)
      phyex_transform_sources(EP_SOURCES ${EP})
      ecbuild_add_executable(
        TARGET ${NAME}_${PRECISION}
        SOURCES ${EP_SOURCES}
        OBJECTS phyex_${PRECISION}
        INCLUDES ${CMAKE_CURRENT_BINARY_DIR}/modules_${PRECISION}
        LIBS
          fiat
          mpi_serial
          parkind_${PRECISION}
          $<$<BOOL:${HAVE_OMP}>:OpenMP::OpenMP_Fortran>
          $<$<BOOL:${HAVE_ACC}>:OpenACC::OpenACC_Fortran>
        OUTPUT_NAME ${NAME}_${PRECISION}.exe
      )
      phyex_transformation_target(${NAME}_${PRECISION})

      if(DEFINED FPP_FLAGS_TESTPROGS)
        target_compile_definitions(${NAME}_${PRECISION} PUBLIC ${FPP_FLAGS_TESTPROGS})
      endif()

      # Compilation option to use for unqualified REALs
      target_compile_options(${NAME}_${PRECISION} PRIVATE
        $<$<COMPILE_LANGUAGE:Fortran>:${PHYEX_COMPILE_OPTIONS_REAL}>)
    endforeach()
  endif()

endforeach()

phyex_transformation_finalize()

#-----------------------------------------------------------------------
# Final packaging
#-----------------------------------------------------------------------

ecbuild_install_project(NAME PHYEX)
ecbuild_print_summary()
