# gather include directories
set(HEADER "flint/flint.h")
if(DEFINED ENV{CPATH})
	#Check for all paths in CPATH
	string(REPLACE ":" ";" TMP_INCLUDE_DIRS $ENV{CPATH})
	find_path(FLINT_INCLUDE_DIR NAMES ${HEADER} PATHS ${TMP_INCLUDE_DIRS} NO_DEFAULT_PATH)

	#If not found, append "/include" and try again
	if(FLINT_INCLUDE_DIR STREQUAL "FLINT_INCLUDE_DIR-NOTFOUND")
		string(REPLACE ":" "/include;" TMP_INCLUDE_DIRS $ENV{CPATH})
		find_path(FLINT_INCLUDE_DIR NAMES ${HEADER} PATHS ${TMP_INCLUDE_DIRS} NO_DEFAULT_PATH)
	endif()

	#Still not found, check default directories
	if(FLINT_INCLUDE_DIR STREQUAL "FLINT_INCLUDE_DIR-NOTFOUND")
		find_path(FLINT_INCLUDE_DIR NAMES ${HEADER})
	endif()

else()
	find_path(FLINT_INCLUDE_DIR NAMES ${HEADER})
endif()

if(DEFINED ENV{LD_LIBRARY_PATH})
	string(REPLACE ":" ";" TMP_LIBRARY_DIRS $ENV{LD_LIBRARY_PATH})
	find_library(FLINT_LIBRARY NAMES flint PATHS ${TMP_LIBRARY_DIRS} NO_DEFAULT_PATH)
	if(FLINT_LIBRARY STREQUAL "FLINT_LIBRARY-NOTFOUND")
		find_library(FLINT_LIBRARY NAMES flint)
	endif()
else()
	find_library(FLINT_LIBRARY NAMES flint)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FLINT DEFAULT_MSG FLINT_LIBRARY FLINT_INCLUDE_DIR)

if (NOT TARGET flint)
	add_library(flint UNKNOWN IMPORTED)
	set_target_properties(flint PROPERTIES IMPORTED_LOCATION ${FLINT_LIBRARY})
	target_include_directories(flint INTERFACE ${FLINT_INCLUDE_DIR})
endif()
