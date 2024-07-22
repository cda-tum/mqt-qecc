# Code adapted from https://fossies.org/linux/llvm/cmake/modules/FindZ3.cmake
include(CheckCXXSourceRuns)
include(FindPackageHandleStandardArgs)

# Function to check Z3's version
function(check_z3_version z3_include z3_lib)
  try_run(
    Z3_RUN_RESULT Z3_COMPILE_RESULT ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/try_z3.cpp
    CMAKE_FLAGS -DINCLUDE_DIRECTORIES:STRING=${z3_include} LINK_LIBRARIES ${z3_lib}
    COMPILE_OUTPUT_VARIABLE COMPILE_OUTPUT RUN_OUTPUT_STDOUT_VARIABLE RUN_OUTPUT
                            RUN_OUTPUT_STDERR_VARIABLE RUN_OUTPUT_STDERR)

  if(NOT Z3_COMPILE_RESULT OR RUN_OUTPUT_STDERR)
    if(NOT Z3_COMPILE_RESULT)
      message(STATUS "Could not compile test program for Z3 version check. Compile output: "
                     ${COMPILE_OUTPUT})
    else()
      message(STATUS "Could not run test program for Z3 version check. Run output: " ${RUN_OUTPUT}
                     " RUN_OUTPUT_STDERR: " ${RUN_OUTPUT_STDERR})
    endif()
    if(z3_include AND EXISTS "${z3_include}/z3_version.h")
      file(STRINGS "${z3_include}/z3_version.h" z3_version_str
           REGEX "^#define[\t ]+Z3_MAJOR_VERSION[\t ]+.*")
      string(REGEX REPLACE "^.*Z3_MAJOR_VERSION[\t ]+([0-9]*).*$" "\\1" Z3_MAJOR
                           "${z3_version_str}")

      file(STRINGS "${z3_include}/z3_version.h" z3_version_str
           REGEX "^#define[\t ]+Z3_MINOR_VERSION[\t ]+.*")
      string(REGEX REPLACE "^.*Z3_MINOR_VERSION[\t ]+([0-9]*).*$" "\\1" Z3_MINOR
                           "${z3_version_str}")

      file(STRINGS "${z3_include}/z3_version.h" z3_version_str
           REGEX "^#define[\t ]+Z3_BUILD_NUMBER[\t ]+.*")
      string(REGEX REPLACE "^.*Z3_BUILD_NUMBER[\t ]+([0-9]*).*$" "\\1" Z3_BUILD "${z3_version_str}")

      set(z3_version_string ${Z3_MAJOR}.${Z3_MINOR}.${Z3_BUILD})
    endif()

    if(NOT z3_version_string)
      message(STATUS "Could not determine Z3 version from z3_version.h")
      return()
    endif()
  else()
    string(REGEX MATCH "(Z3 )?([0-9]+.[0-9]+.[0-9]+.[0-9]+)" Z3_VERSION_STRING ${RUN_OUTPUT})
    set(z3_version_string ${CMAKE_MATCH_2})
  endif()

  find_package_check_version(${z3_version_string} suitable_version RESULT_MESSAGE_VARIABLE reason)

  if(suitable_version)
    set(FOUND_SUITABLE_VERSION
        TRUE
        PARENT_SCOPE)
    set(Z3_VERSION_STRING
        ${z3_version_string}
        PARENT_SCOPE)
  else()
    message(STATUS "${reason}")
  endif()

endfunction(check_z3_version)

set(Z3_ROOT
    ""
    CACHE PATH "Root of Z3 distribution.")
if(DEFINED ENV{Z3_ROOT})
  set(Z3_ROOT $ENV{Z3_ROOT})
  message(STATUS "Z3_ROOT from environment: ${Z3_ROOT}")
endif()

# if Z3_ROOT is provided, check there first
if(NOT ${Z3_ROOT} STREQUAL "")
  find_path(
    Z3_CXX_INCLUDE_DIRS
    NAMES z3.h z3++.h
    NO_DEFAULT_PATH
    PATHS ${Z3_ROOT}/include
    PATH_SUFFIXES libz3 z3)

  find_library(
    Z3_LIBRARIES
    NAMES z3 libz3
    NO_DEFAULT_PATH
    PATHS ${Z3_ROOT}
    PATH_SUFFIXES lib bin)

  if(Z3_CXX_INCLUDE_DIRS AND Z3_LIBRARIES)
    message(STATUS "Z3_ROOT provided and includes and libraries found.")
    message(VERBOSE "Z3_CXX_INCLUDE_DIRS: ${Z3_CXX_INCLUDE_DIRS}")
    message(VERBOSE "Z3_LIBRARIES: ${Z3_LIBRARIES}")
    check_z3_version(${Z3_CXX_INCLUDE_DIRS} ${Z3_LIBRARIES})
  endif()
endif()

# see if a config file is available
if(NOT FOUND_SUITABLE_VERSION)
  unset(Z3_CXX_INCLUDE_DIRS CACHE)
  unset(Z3_LIBRARIES CACHE)

  find_package(Z3 CONFIG QUIET)
  if(Z3_FOUND)
    message(STATUS "Found Z3 includes and libraries from config file")
    message(VERBOSE "Z3_CXX_INCLUDE_DIRS: ${Z3_CXX_INCLUDE_DIRS}")
    message(VERBOSE "Z3_LIBRARIES: ${Z3_LIBRARIES}")
    set(FOUND_SUITABLE_VERSION TRUE)
    set(Z3_VERSION_STRING ${Z3_VERSION})
  endif()
endif()

# if Z3 has not been found yet, look in the system paths
if(NOT FOUND_SUITABLE_VERSION)
  unset(Z3_CXX_INCLUDE_DIRS CACHE)
  unset(Z3_LIBRARIES CACHE)

  find_path(
    Z3_CXX_INCLUDE_DIRS
    NAMES z3.h z3++.h
    PATH_SUFFIXES libz3 z3)
  find_library(
    Z3_LIBRARIES
    NAMES z3 libz3
    PATH_SUFFIXES lib bin)

  if(Z3_CXX_INCLUDE_DIRS AND Z3_LIBRARIES)
    message(STATUS "Found Z3 includes and libraries in system paths.")
    message(VERBOSE "Z3_CXX_INCLUDE_DIRS: ${Z3_CXX_INCLUDE_DIRS}")
    message(VERBOSE "Z3_LIBRARIES: ${Z3_LIBRARIES}")
    check_z3_version(${Z3_CXX_INCLUDE_DIRS} ${Z3_LIBRARIES})
  endif()
endif()

# if it is still not found, try to find it with Python as a last resort
if(NOT FOUND_SUITABLE_VERSION)
  unset(Z3_CXX_INCLUDE_DIRS CACHE)
  unset(Z3_LIBRARIES CACHE)

  set(PYTHON_FIND_VIRTUALENV FIRST)
  find_package(Python COMPONENTS Interpreter Development.Module)
  if(Python_FOUND)
    execute_process(
      COMMAND ${Python_EXECUTABLE} -c "import os, z3; print(os.path.dirname(z3.__file__))"
      OUTPUT_VARIABLE Z3_PYTHON_ROOT)
    string(STRIP "${Z3_PYTHON_ROOT}" Z3_PYTHON_ROOT)
    message(STATUS "Z3_PYTHON_ROOT: ${Z3_PYTHON_ROOT}")

    if(Z3_PYTHON_ROOT)
      find_path(
        Z3_CXX_INCLUDE_DIRS
        NAMES z3.h z3++.h
        NO_DEFAULT_PATH
        PATHS ${Z3_PYTHON_ROOT}
        PATH_SUFFIXES libz3 z3 include)

      find_library(
        Z3_LIBRARIES
        NAMES z3 libz3
        NO_DEFAULT_PATH
        PATHS ${Z3_PYTHON_ROOT}
        PATH_SUFFIXES lib bin)

      if(Z3_CXX_INCLUDE_DIRS AND Z3_LIBRARIES)
        message(STATUS "Found Z3 includes and libraries from Python installation.")
        message(VERBOSE "Z3_CXX_INCLUDE_DIRS: ${Z3_CXX_INCLUDE_DIRS}")
        message(VERBOSE "Z3_LIBRARIES: ${Z3_LIBRARIES}")
        check_z3_version(${Z3_CXX_INCLUDE_DIRS} ${Z3_LIBRARIES})
      endif()
    endif()
  endif()
endif()

if(NOT FOUND_SUITABLE_VERSION)
  if(Z3_CXX_INCLUDE_DIRS AND Z3_LIBRARIES)
    message(STATUS "Found include and library directories but could not find a suitable Z3 version")
  endif()
  set(Z3_VERSION_STRING "0.0.0")
endif()

find_package_handle_standard_args(
  Z3
  REQUIRED_VARS Z3_LIBRARIES Z3_CXX_INCLUDE_DIRS
  VERSION_VAR Z3_VERSION_STRING)

if(Z3_FOUND)
  if(NOT TARGET z3::z3lib)
    add_library(z3::z3lib INTERFACE IMPORTED GLOBAL)
    target_include_directories(z3::z3lib INTERFACE ${Z3_CXX_INCLUDE_DIRS})
    target_link_libraries(z3::z3lib INTERFACE ${Z3_LIBRARIES})
  endif()
  add_compile_definitions(Z3_FOUND)
endif()

mark_as_advanced(Z3_CXX_INCLUDE_DIRS Z3_LIBRARIES Z3_VERSION_STRING)
