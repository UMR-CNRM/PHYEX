cmake_minimum_required(VERSION 3.10)
find_package(ecbuild 3.11 REQUIRED HINTS ${CMAKE_SOURCE_DIR})
project(PHYEX VERSION 1.0 LANGUAGES Fortran)

include(phyex_macros)

#-----------------------------------------------------------------------
# Options
#-----------------------------------------------------------------------

ecbuild_add_option(FEATURE DOUBLE_PRECISION
                   DESCRIPTION "Compile on double precision"
                   DEFAULT ON)

ecbuild_add_option(FEATURE SINGLE_PRECISION
                   DESCRIPTION "Compile on single precision"
                   DEFAULT OFF)

ecbuild_add_option(FEATURE PHYEX_BUILD_PROGS
                   DESCRIPTION "Build PHYEX program executables"
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

# Fetch fiat from GitHub
include(FetchContent)
FetchContent_Declare(
  fiat
  GIT_REPOSITORY https://github.com/ecmwf-ifs/fiat
  GIT_TAG        main
)
FetchContent_MakeAvailable(fiat)

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
      mpi_serial
      parkind_${PRECISION}
    PUBLIC_INCLUDES ${PHYEX_HEADERS}
    INSTALL_HEADERS ALL
    HEADER_DESTINATION include/${PRECISION}
    PUBLIC_INCLUDES ${PHYEX_INCLUDEDIRS}
  ) 

  # Compilation option to use for unqualified REALs
  target_compile_options(phyex_${PRECISION} PUBLIC
    $<$<COMPILE_LANGUAGE:Fortran>:${PHYEX_COMPILE_OPTIONS_REAL}>)

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
      ecbuild_add_executable(
        TARGET ${NAME}_${PRECISION}
        SOURCES ${EP}
        OBJECTS phyex_${PRECISION}
        INCLUDES ${CMAKE_CURRENT_BINARY_DIR}/modules_${PRECISION}
        LIBS
          fiat
          mpi_serial
          parkind_${PRECISION}
        OUTPUT_NAME ${NAME}_${PRECISION}.exe
      )
      if(DEFINED FPP_FLAGS_TESTPROGS)
        target_compile_definitions(${NAME}_${PRECISION} PUBLIC ${FPP_FLAGS_TESTPROGS})
      endif()

      # Compilation option to use for unqualified REALs
      target_compile_options(${NAME}_${PRECISION} PUBLIC
        $<$<COMPILE_LANGUAGE:Fortran>:${PHYEX_COMPILE_OPTIONS_REAL}>)
    endforeach()
  endif()

endforeach()

#-----------------------------------------------------------------------
# Final packaging
#-----------------------------------------------------------------------

ecbuild_install_project(NAME PHYEX)
ecbuild_print_summary()
