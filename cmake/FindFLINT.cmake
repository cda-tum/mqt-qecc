set(FLINT_ROOT
    ""
    CACHE PATH "Root of flint distribution.")
if(DEFINED ENV{FLINT_ROOT})
  set(FLINT_ROOT $ENV{FLINT_ROOT})
  message(STATUS "FLINT_ROOT from environment: ${FLINT_ROOT}")
endif()

# if FLINT_ROOT is provided, check there first
if(NOT ${FLINT_ROOT} STREQUAL "")
  message(STATUS "FLINT_ROOT provided.")
  find_path(
    FLINT_INCLUDE_DIR
    NAMES flint.h flintxx.h
    NO_DEFAULT_PATH
    PATHS ${FLINT_ROOT}
    PATH_SUFFIXES libflint flint)

  if(FLINT_INCLUDE_DIR)
    message(VERBOSE "FLINT_INCLUDE_DIR: ${FLINT_INCLUDE_DIR}")
  endif()

  find_library(
    FLINT_LIBRARY
    NAMES flint libflint
    NO_DEFAULT_PATH
    PATHS ${FLINT_ROOT}
    PATH_SUFFIXES lib bin)

  if(FLINT_LIBRARY)
    message(VERBOSE "FLINT_LIBRARY: ${FLINT_LIBRARY}")
  endif()
endif()

# if FLINT_ROOT is not provided or does not contain flint, check system
if(NOT FLINT_INCLUDE_DIR OR NOT FLINT_LIBRARY)
  # gather include directories
  set(HEADER "flint/flint.h")
  if(DEFINED ENV{CPATH})
    # Check for all paths in CPATH
    string(REPLACE ":" ";" TMP_INCLUDE_DIRS $ENV{CPATH})
    find_path(
      FLINT_INCLUDE_DIR
      NAMES ${HEADER}
      PATHS ${TMP_INCLUDE_DIRS}
      NO_DEFAULT_PATH)

    # If not found, append "/include" and try again
    if(FLINT_INCLUDE_DIR STREQUAL "FLINT_INCLUDE_DIR-NOTFOUND")
      string(REPLACE ":" "/include;" TMP_INCLUDE_DIRS $ENV{CPATH})
      find_path(
        FLINT_INCLUDE_DIR
        NAMES ${HEADER}
        PATHS ${TMP_INCLUDE_DIRS}
        NO_DEFAULT_PATH)
    endif()

    # Still not found, check default directories
    if(FLINT_INCLUDE_DIR STREQUAL "FLINT_INCLUDE_DIR-NOTFOUND")
      find_path(FLINT_INCLUDE_DIR NAMES ${HEADER})
    endif()

  else()
    find_path(FLINT_INCLUDE_DIR NAMES ${HEADER})
  endif()

  if(DEFINED ENV{LD_LIBRARY_PATH})
    string(REPLACE ":" ";" TMP_LIBRARY_DIRS $ENV{LD_LIBRARY_PATH})
    find_library(
      FLINT_LIBRARY
      NAMES flint
      PATHS ${TMP_LIBRARY_DIRS}
      NO_DEFAULT_PATH)
    if(FLINT_LIBRARY STREQUAL "FLINT_LIBRARY-NOTFOUND")
      find_library(FLINT_LIBRARY NAMES flint)
    endif()
  else()
    find_library(FLINT_LIBRARY NAMES flint)
  endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FLINT DEFAULT_MSG FLINT_LIBRARY
                                  FLINT_INCLUDE_DIR)

if(NOT TARGET flint)
  add_library(flint UNKNOWN IMPORTED)
  set_target_properties(flint PROPERTIES IMPORTED_LOCATION ${FLINT_LIBRARY})
  target_include_directories(flint INTERFACE ${FLINT_INCLUDE_DIR})
endif()
