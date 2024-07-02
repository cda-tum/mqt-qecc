# Add a feature summary for this package
include(FeatureSummary)
set_package_properties(
  MPFR PROPERTIES
  DESCRIPTION
    "MPFR is a C library for multiple-precision floating-point computations with correct rounding."
  URL "https://www.mpfr.org/")

# Try finding the package with pkg-config
find_package(PkgConfig QUIET)
pkg_check_modules(PKG QUIET mpfr)

# Try to locate the libraries and their headers, using pkg-config hints
find_path(MPFR_INCLUDE_DIR mpfr.h HINTS ${PKG_mpfr_INCLUDEDIR})
find_library(MPFR_LIB mpfr HINTS ${PKG_mpfr_LIBDIR})

# Remove these variables from cache inspector
mark_as_advanced(MPFR_INCLUDE_DIR MPFR_LIB)

# Report if package was found
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MPFR DEFAULT_MSG MPFR_INCLUDE_DIR MPFR_LIB)

# Set target
if(MPFR_FOUND)
  if(NOT TARGET MPFR::mpfr)
    add_library(MPFR::mpfr UNKNOWN IMPORTED)
    set_target_properties(MPFR::mpfr PROPERTIES IMPORTED_LOCATION ${MPFR_LIB}
                                                INTERFACE_INCLUDE_DIRECTORIES ${MPFR_INCLUDE_DIR})
  endif()
endif()
