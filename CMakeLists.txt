cmake_minimum_required(VERSION 3.14...3.22)

project(qecc
        LANGUAGES CXX
        DESCRIPTION "QECC - A MQT Library for Quantum Error Correcting Codes"
        )

set(CMAKE_CXX_STANDARD 20)

#add_compile_options(-fsanitize=address)
#add_link_options(-fsanitize=address)
#add_compile_options(-fsanitize=undefined)
#add_link_options(-fsanitize=undefined)
# enable organization of targets into folders

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# configuration options
option(DEPLOY "Configure for deployment")
option(BINDINGS "Configure for building Python bindings")
option(COVERAGE "Configure for coverage report generation")
option(GENERATE_POSITION_INDEPENDENT_CODE "Generate position independent code")
option(BUILD_QECC_TESTS "Also build tests for project")

if (DEFINED ENV{DEPLOY})
    set(DEPLOY $ENV{DEPLOY} CACHE BOOL "Use deployment configuration from environment" FORCE)
    message(STATUS "Setting deployment configuration to '${DEPLOY}' from environment")
endif ()

# set deployment specific options
if (DEPLOY)
    # set the macOS deployment target appropriately
    set(CMAKE_OSX_DEPLOYMENT_TARGET "10.15" CACHE STRING "" FORCE)
endif ()

# build type settings
set(default_build_type "Release")
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif ()

macro(check_submodule_present MODULENAME)
    if (NOT EXISTS "${PROJECT_SOURCE_DIR}/extern/${MODULENAME}/CMakeLists.txt")
        message(FATAL_ERROR "${MODULENAME} submodule not cloned properly. Please run `git submodule update --init --recursive` from the main project directory")
    endif ()
endmacro()

check_submodule_present(json)
check_submodule_present(pybind11)
check_submodule_present(pybind11_json)

# Add path for custom modules
list (APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# add main library code
add_subdirectory(src)

# add test code
if (BUILD_QECC_TESTS)
    check_submodule_present(googletest)

    enable_testing()
    include(GoogleTest)
    add_subdirectory(test)
endif ()

# add Python binding code
if(BINDINGS)
    add_subdirectory(mqt/qecc)
endif()

# files containing codes for tests
FILE(COPY examples/test/ DESTINATION "${CMAKE_BINARY_DIR}/test/resources/codes/" PATTERN "*")