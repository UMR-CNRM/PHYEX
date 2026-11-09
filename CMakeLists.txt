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
# Optional parts:
#   - python wheel (-DENABLE_PYTHON_WHEEL=ON): cmake/phyex_python_wheel.cmake
#   - transformation of the PHYEX source code at build time by pyfortool,
#     with -DPHYEX_TRANSFORMATION=<name> (or the environment variable of the
#     same name), e.g. -DPHYEX_TRANSFORMATION=IAL: done by src/phyex (see
#     src/phyex/cmake/phyex_transformation.cmake)
#   - check of the coding norms of the transformed code
#     (-DENABLE_CHECK_NORMS=ON, default): cmake/phyex_check_norms.cmake
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
# Options
#-----------------------------------------------------------------------

ecbuild_add_option(FEATURE PYTHON_WHEEL
                   DESCRIPTION "Build the PHYEX python wheel (installs the python tools in the build directory)"
                   DEFAULT OFF)
ecbuild_add_option(FEATURE CHECK_NORMS
                   DESCRIPTION "Check the coding norms of the transformed source code (PHYEX_TRANSFORMATION)"
                   DEFAULT ON)

#-----------------------------------------------------------------------
# Step 3: PHYEX (including the optional transformation of the source code)
#-----------------------------------------------------------------------

add_subdirectory(src/phyex ${CMAKE_BINARY_DIR}/phyex)

#-----------------------------------------------------------------------
# Step 4: python wheel (optional)
#-----------------------------------------------------------------------

if(HAVE_PYTHON_WHEEL)
  include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/phyex_python_wheel.cmake)
endif()

#-----------------------------------------------------------------------
# Step 5: coding norms of the transformed source code (optional)
#-----------------------------------------------------------------------

# The phyex_transformation target exists only if a transformation is active
if(HAVE_CHECK_NORMS AND TARGET phyex_transformation)
  include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/phyex_check_norms.cmake)
elseif(HAVE_CHECK_NORMS)
  ecbuild_debug("No source code transformation: CHECK_NORMS has no effect")
endif()

ecbuild_bundle_finalize()
