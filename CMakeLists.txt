cmake_minimum_required(VERSION 3.19...3.27)

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
