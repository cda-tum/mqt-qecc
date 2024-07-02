set(MPFR_ROOT
    ""
    CACHE PATH "Root of mpfr distribution.")
if(DEFINED ENV{MPFR_ROOT})
  set(MPFR_ROOT $ENV{MPFR_ROOT})
  message(STATUS "MPFR_ROOT from environment: ${MPFR_ROOT}")
endif()

# if MPFR_ROOT is provided, check there first
if(NOT ${MPFR_ROOT} STREQUAL "")
  message(STATUS "MPFR_ROOT provided.")
  find_path(
    MPFR_INCLUDE_DIR
    NAMES mpfr.h
    NO_DEFAULT_PATH
    PATHS ${MPFR_ROOT}
    PATH_SUFFIXES libmpfr mpfr)

  if(MPFR_INCLUDE_DIR)
    message(VERBOSE "MPFR_INCLUDE_DIR: ${MPFR_INCLUDE_DIR}")
  endif()

  find_library(
    MPFR_LIBRARY
    NAMES mpfr libmpfr
    NO_DEFAULT_PATH
    PATHS ${MPFR_ROOT}
    PATH_SUFFIXES lib bin)

  if(MPFR_LIBRARY)
    message(VERBOSE "MPFR_LIBRARY: ${MPFR_LIBRARY}")
  endif()
endif()

# if MPFR_ROOT is not provided or does not contain mpfr, check system
if(NOT MPFR_INCLUDE_DIR OR NOT MPFR_LIBRARY)
  # gather include directories
  set(HEADER "mpfr.h")
  if(DEFINED ENV{CPATH})
    # Check for all paths in CPATH
    string(REPLACE ":" ";" TMP_INCLUDE_DIRS $ENV{CPATH})
    find_path(
      MPFR_INCLUDE_DIR
      NAMES ${HEADER}
      PATHS ${TMP_INCLUDE_DIRS}
      NO_DEFAULT_PATH)

    # If not found, append "/include" and try again
    if(MPFR_INCLUDE_DIR STREQUAL "MPFR_INCLUDE_DIR-NOTFOUND")
      string(REPLACE ":" "/include;" TMP_INCLUDE_DIRS $ENV{CPATH})
      find_path(
        MPFR_INCLUDE_DIR
        NAMES ${HEADER}
        PATHS ${TMP_INCLUDE_DIRS}
        NO_DEFAULT_PATH)
    endif()

    # Still not found, check default directories
    if(MPFR_INCLUDE_DIR STREQUAL "MPFR_INCLUDE_DIR-NOTFOUND")
      find_path(MPFR_INCLUDE_DIR NAMES ${HEADER})
    endif()

  else()
    find_path(MPFR_INCLUDE_DIR NAMES ${HEADER})
  endif()

  if(DEFINED ENV{LD_LIBRARY_PATH})
    string(REPLACE ":" ";" TMP_LIBRARY_DIRS $ENV{LD_LIBRARY_PATH})
    find_library(
      MPFR_LIBRARY
      NAMES mpfr
      PATHS ${TMP_LIBRARY_DIRS}
      NO_DEFAULT_PATH)
    if(MPFR_LIBRARY STREQUAL "MPFR_LIBRARY-NOTFOUND")
      find_library(MPFR_LIBRARY NAMES mpfr)
    endif()
  else()
    find_library(MPFR_LIBRARY NAMES mpfr)
  endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MPFR DEFAULT_MSG MPFR_LIBRARY MPFR_INCLUDE_DIR)

if(NOT TARGET MPFR::mpfr)
  add_library(MPFR::mpfr UNKNOWN IMPORTED)
  set_target_properties(MPFR::mpfr PROPERTIES IMPORTED_LOCATION ${MPFR_LIBRARY})
  target_include_directories(MPFR::mpfr INTERFACE ${MPFR_INCLUDE_DIR})
endif()
