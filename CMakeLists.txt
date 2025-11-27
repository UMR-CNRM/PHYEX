cmake_minimum_required(VERSION 3.10)
find_package(ecbuild 3.7 REQUIRED)
project(PHYEX VERSION 1.0 LANGUAGES Fortran)

include( phyex_macros )
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

set(GLOB PHYEX_HEADERS
    ${CMAKE_SOURCE_DIR}/turb/*h
    ${CMAKE_SOURCE_DIR}/micro/*h
)
if(EXISTS ${CMAKE_SOURCE_DIR}/conv)
  file(GLOB CONV_HEADERS ${CMAKE_SOURCE_DIR}/conv/*h)
  list(APPEND PHYEX_HEADERS ${CONV_HEADERS})
endif()
if(EXISTS ${CMAKE_SOURCE_DIR}/pyphyex.F90)
  list(APPEND PHYEX_ENTRYPOINTS pyphyex.F90)
endif()

set(PHYEX_INCLUDEDIRS
  ${CMAKE_SOURCE_DIR}/aux
  ${CMAKE_SOURCE_DIR}/micro
  ${CMAKE_SOURCE_DIR}/turb
)
if(EXISTS ${CMAKE_SOURCE_DIR}/conv)
  list(APPEND PHYEX_INCLUDEDIRS ${CMAKE_SOURCE_DIR}/conv)
endif()

#-----------------------------------------------------------------------
# Library
#-----------------------------------------------------------------------
foreach(PRECISION ${precision_list})
  
ecbuild_add_library(
  TARGET phyex_${PRECISION}
  SOURCES
    ${PHYEX_ENTRYPOINTS}
    ${PHYEX_SOURCES}
  PUBLIC_LIBS
    fiat
    mpi_serial
    parkind_${PRECISION}
  PUBLIC_INCLUDES     ${PHYEX_HEADERS}
  INSTALL_HEADERS     ALL
  HEADER_DESTINATION include/${PRECISION}

  PRIVATE_INCLUDES
    ${PHYEX_INCLUDEDIRS}

) 

install(TARGETS phyex_${PRECISION}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )

fortran_modules_install(
        TARGET phyex_${PRECISION}
        MODULE_DIR ${CMAKE_CURRENT_BINARY_DIR}/modules_${PRECISION}
	INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/modules/${PRECISION}
  )

endforeach()



###############################################################################
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

#-----------------------------------------------------------------------
# Final packaging
#-----------------------------------------------------------------------
ecbuild_install_project(NAME PHYEX)
ecbuild_print_summary()
