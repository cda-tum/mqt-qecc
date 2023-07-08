cmake_minimum_required(VERSION 3.19)

# This avoids googletest complaining that this (IPO) policy is not set
cmake_policy(SET CMP0069 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0069 NEW)

project(
  qecc
  LANGUAGES CXX
  DESCRIPTION "QECC - An MQT Library for Quantum Error Correcting Codes")

# check whether `modulename` is correctly cloned in the `/extern` directory.
macro(CHECK_SUBMODULE_PRESENT modulename)
  if(NOT EXISTS "${PROJECT_SOURCE_DIR}/extern/${modulename}/CMakeLists.txt")
    message(
      FATAL_ERROR
        "${modulename} submodule not cloned properly. \
        Please run `git submodule update --init --recursive` \
        from the main project directory")
  endif()
endmacro()

check_submodule_present(mqt-core)

# Add path for custom modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# add main library code
add_subdirectory(src)

# add test code
option(BUILD_QECC_TESTS "Also build tests for project")
if(BUILD_QECC_TESTS)
  enable_testing()
  include(GoogleTest)
  add_subdirectory(test)
endif()

# add Python binding code
option(BINDINGS "Configure for building Python bindings")
if(BINDINGS)
  add_subdirectory(mqt/qecc)
endif()
