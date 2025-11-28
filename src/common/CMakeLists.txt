cmake_minimum_required(VERSION 3.10)
find_package(ecbuild 3.7 REQUIRED)
project(PHYEX VERSION 1.0 LANGUAGES Fortran)

include(phyex_macros)

#-----------------------------------------------------------------------
# Options
#-----------------------------------------------------------------------

ecbuild_add_option( FEATURE DOUBLE_PRECISION
                    DESCRIPTION "Compile on double precision"
                    DEFAULT ON )

ecbuild_add_option( FEATURE SINGLE_PRECISION
                    DESCRIPTION "Compile on single precision"
                    DEFAULT OFF ) 

ecbuild_add_option( FEATURE PHYEX_BUILD_PROGS
                    DESCRIPTION "Build PHYEX program executables"
                    DEFAULT OFF )

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

file(GLOB PHYEX_SOURCES
    ${CMAKE_SOURCE_DIR}/aux/*90
    ${CMAKE_SOURCE_DIR}/turb/*90
    ${CMAKE_SOURCE_DIR}/micro/*90
) 
set(PHYEX_ENTRYPOINTS
    aux/ini_phyex.F90
    micro/rain_ice.F90
    turb/shallow_mf.F90
    turb/turb.F90
    micro/ice_adjust.F90
    micro/lima.F90
    micro/lima_adjust_split.F90
)
if(EXISTS ${CMAKE_SOURCE_DIR}/pyphyex.F90)
  list(APPEND PHYEX_ENTRYPOINTS pyphyex.F90)
endif()

# diry headers and include dirs : 
set(GLOB PHYEX_HEADERS
    ${CMAKE_SOURCE_DIR}/turb/*h
    ${CMAKE_SOURCE_DIR}/micro/*h
)
if(EXISTS ${CMAKE_SOURCE_DIR}/conv)
  file(GLOB CONV_HEADERS ${CMAKE_SOURCE_DIR}/conv/*h)
  list(APPEND PHYEX_HEADERS ${CONV_HEADERS})
endif()

set(PHYEX_INCLUDEDIRS
  ${CMAKE_SOURCE_DIR}/aux
  ${CMAKE_SOURCE_DIR}/micro
  ${CMAKE_SOURCE_DIR}/turb
)
if(EXISTS ${CMAKE_SOURCE_DIR}/conv)
  list(APPEND PHYEX_INCLUDEDIRS ${CMAKE_SOURCE_DIR}/conv)
endif()

if(HAVE_PHYEX_BUILD_PROGS)
  # Source file directories
  file(GLOB PHYEXPROGS_SOURCES
      ${CMAKE_SOURCE_DIR}/support/*90
      ${CMAKE_SOURCE_DIR}/progs/*90
  ) 

  # Entry points
  file(GLOB PHYEXPROGS_ENTRYPOINTS progs/main_*.F90)

  # We must excludes from PROG_FILES objects used as entry points
  list(REMOVE_ITEM PHYEXPROGS_SOURCES ${PHYEXPROGS_ENTRYPOINTS})
endif()

#---------------
#  Loop on precisions
#---------------

foreach(PRECISION ${precision_list})
  
  #---------------
  #  Internal library (to prevent multiple compilations)
  #---------------

  ecbuild_add_library(
    TARGET internal_phyex_${PRECISION}
    SOURCES ${PHYEX_SOURCES}
    PUBLIC_LIBS
      fiat
      mpi_serial
      parkind_${PRECISION}
    TYPE STATIC
    NOINSTALL
    PRIVATE_INCLUDES ${PHYEX_INCLUDEDIRS}
  )

  # For module files to be available
  fortran_modules_install(
      TARGET internal_phyex_${PRECISION}
      MODULE_DIR ${CMAKE_CURRENT_BINARY_DIR}/modules_${PRECISION}
  )

  #---------------
  #  Public library
  #---------------

  ecbuild_add_library(
    TARGET phyex_${PRECISION}
    SOURCES ${PHYEX_ENTRYPOINTS}
    PRIVATE_LIBS internal_phyex_${PRECISION}
    PUBLIC_LIBS
      fiat
      mpi_serial
      parkind_${PRECISION}
    PUBLIC_INCLUDES ${PHYEX_HEADERS}
    INSTALL_HEADERS ALL
    HEADER_DESTINATION include/${PRECISION}
    PRIVATE_INCLUDES ${PHYEX_INCLUDEDIRS}
  ) 
  
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
    # Buid a lib with files common to all the executables
    ecbuild_add_library(
      TARGET phyexprogs_${PRECISION}
      SOURCES ${PHYEXPROGS_SOURCES}
      PUBLIC_LIBS
        fiat
        mpi_serial
        parkind_${PRECISION}
      TYPE STATIC
      NOINSTALL
    ) 
    if(DEFINED FPP_FLAGS_TESTPROGS)
      target_compile_definitions(phyexprogs_${PRECISION} PUBLIC ${FPP_FLAGS_TESTPROGS})
    endif()
    
    # Build all executables
    foreach(EP ${PHYEXPROGS_ENTRYPOINTS})
      get_filename_component(NAME ${EP} NAME_WLE)
      ecbuild_add_executable(
        TARGET ${NAME}_${PRECISION}
        SOURCES ${EP}
        LIBS
          internal_phyex_${PRECISION}
          phyexprogs_${PRECISION}
          fiat
          mpi_serial
          parkind_${PRECISION}
        OUTPUT_NAME ${NAME}_${PRECISION}.exe
      )
      if(DEFINED FPP_FLAGS_TESTPROGS)
        target_compile_definitions(${NAME}_${PRECISION} PUBLIC ${FPP_FLAGS_TESTPROGS})
      endif()
    endforeach()
  endif()

endforeach()

#-----------------------------------------------------------------------
# Final packaging
#-----------------------------------------------------------------------

ecbuild_install_project(NAME PHYEX)
ecbuild_print_summary()
