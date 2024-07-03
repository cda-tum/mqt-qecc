# Add a feature summary for this package
include(FeatureSummary)
set_package_properties(
  mpfr PROPERTIES
  DESCRIPTION "MPFR: Multiple Precision Floating-Point Reliable Library"
  URL "https://www.mpfr.org")

# Try finding the package with pkg-config
find_package(PkgConfig QUIET)
pkg_check_modules(PKG QUIET mpfr)

find_path(
  MPFR_INCLUDE_DIR
  NAMES mpfr.h
  PATH_SUFFIXES libmpfr mpfr
  HINTS ${PKG_mpfr_INCLUDEDIR})

find_library(
  MPFR_LIBRARY
  NAMES mpfr libmpfr
  PATH_SUFFIXES lib bin
  HINTS ${PKG_mpfr_LIBDIR})

# Remove these variables from cache inspector
mark_as_advanced(MPFR_INCLUDE_DIR MPFR_LIBRARY)

# Report if package was found
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MPFR DEFAULT_MSG MPFR_LIBRARY MPFR_INCLUDE_DIR)

# Set targets
if(MPFR_FOUND)
  add_library(MPFR::mpfr UNKNOWN IMPORTED)
  set_target_properties(MPFR::mpfr PROPERTIES IMPORTED_LOCATION ${MPFR_LIBRARY})
  target_include_directories(MPFR::mpfr INTERFACE ${MPFR_INCLUDE_DIR})
endif()
