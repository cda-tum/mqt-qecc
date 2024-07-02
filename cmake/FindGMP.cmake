# File adapted from
# https://github.com/dune-project/dune-common/blob/master/cmake/modules/FindGMP.cmake

#[=======================================================================[.rst:
FindGMP
-------

Find the GNU MULTI-Precision Bignum (GMP) library
and the corresponding C++ bindings GMPxx.

This module searches for both libraries and only considers the package
found if both can be located. It then defines separate targets for the C
and the C++ library.

Imported Targets
^^^^^^^^^^^^^^^^

This module provides the following imported targets, if found:

``GMP::gmp``
  Library target of the C library.
``GMP::gmpxx``
  Library target of the C++ library, which also links to the C library.

Result Variables
^^^^^^^^^^^^^^^^

This will define the following variables:

``GMP_FOUND``
  True if the GMP library, the GMPxx headers and
  the GMPxx library were found.

Cache Variables
^^^^^^^^^^^^^^^

You may set the following variables to modify the behaviour of
this module:

``GMP_INCLUDE_DIR``
  The directory containing ``gmp.h``.
``GMP_LIB``
  The path to the gmp library.
``GMPXX_INCLUDE_DIR``
  The directory containing ``gmpxx.h``.
``GMPXX_LIB``
  The path to the gmpxx library.

#]=======================================================================]

# Add a feature summary for this package
include(FeatureSummary)
set_package_properties(
  GMP PROPERTIES
  DESCRIPTION "GNU multi-precision library"
  URL "https://gmplib.org")

# Try finding the package with pkg-config
find_package(PkgConfig QUIET)
pkg_check_modules(PKG QUIET gmp gmpxx)

# Try to locate the libraries and their headers, using pkg-config hints
find_path(GMP_INCLUDE_DIR gmp.h HINTS ${PKG_gmp_INCLUDEDIR})
find_library(GMP_LIB gmp HINTS ${PKG_gmp_LIBDIR})

find_path(GMPXX_INCLUDE_DIR gmpxx.h HINTS ${PKG_gmpxx_INCLUDEDIR})
find_library(GMPXX_LIB gmpxx HINTS ${PKG_gmpxx_LIBDIR})

# Remove these variables from cache inspector
mark_as_advanced(GMP_INCLUDE_DIR GMP_LIB GMPXX_INCLUDE_DIR GMPXX_LIB)

# Report if package was found
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GMP DEFAULT_MSG GMPXX_LIB GMPXX_INCLUDE_DIR GMP_INCLUDE_DIR
                                  GMP_LIB)

# Set targets
if(GMP_FOUND)
  # C library
  if(NOT TARGET GMP::gmp)
    add_library(GMP::gmp UNKNOWN IMPORTED)
    set_target_properties(GMP::gmp PROPERTIES IMPORTED_LOCATION ${GMP_LIB}
                                              INTERFACE_INCLUDE_DIRECTORIES ${GMP_INCLUDE_DIR})
  endif()

  # C++ library, which requires a link to the C library
  if(NOT TARGET GMP::gmpxx)
    add_library(GMP::gmpxx UNKNOWN IMPORTED)
    set_target_properties(
      GMP::gmpxx
      PROPERTIES IMPORTED_LOCATION ${GMPXX_LIB}
                 INTERFACE_INCLUDE_DIRECTORIES ${GMPXX_INCLUDE_DIR}
                 INTERFACE_LINK_LIBRARIES GMP::gmp)
  endif()
endif()
