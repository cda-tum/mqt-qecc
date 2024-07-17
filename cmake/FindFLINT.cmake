find_package(GMP REQUIRED)
find_package(MPFR REQUIRED)

# Add a feature summary for this package
include(FeatureSummary)
set_package_properties(
  flint PROPERTIES
  DESCRIPTION "FLINT: Fast Library for Number Theory"
  URL "https://flintlib.org")

# Try finding the package with pkg-config
find_package(PkgConfig QUIET)
pkg_check_modules(PKG QUIET flint)

find_path(
  FLINT_INCLUDE_DIR
  NAMES flint.h flintxx.h
  PATH_SUFFIXES libflint flint
  HINTS ${PKG_flint_INCLUDEDIR})

find_library(
  FLINT_LIBRARY
  NAMES flint libflint
  PATH_SUFFIXES lib bin
  HINTS ${PKG_flint_LIBDIR})

# Remove these variables from cache inspector
mark_as_advanced(FLINT_INCLUDE_DIR FLINT_LIBRARY)

# Report if package was found
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FLINT DEFAULT_MSG FLINT_LIBRARY FLINT_INCLUDE_DIR)

# Set targets
if(FLINT_FOUND)
  add_library(flint UNKNOWN IMPORTED)
  set_target_properties(flint PROPERTIES IMPORTED_LOCATION ${FLINT_LIBRARY})

  # Strip the trailing "flint" from the include directory
  if("${FLINT_INCLUDE_DIR}" MATCHES "flint$")
    get_filename_component(FLINT_INCLUDE_DIR ${FLINT_INCLUDE_DIR} DIRECTORY)
  endif()

  target_include_directories(flint INTERFACE ${FLINT_INCLUDE_DIR})
  target_link_libraries(flint INTERFACE GMP::gmp GMP::gmpxx MPFR::mpfr)
endif()
