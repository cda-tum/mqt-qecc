if(NOT SKBUILD)
  message(
    NOTICE
    "\
  This CMake file is meant to be executed using 'scikit-build'. Running
  it directly will almost certainly not produce the desired result. If
  you are a user trying to install this package, please use the command
  below, which will install all necessary build dependencies, compile
  the package in an isolated environment, and then install it.
  =====================================================================
   $ pip install .
  =====================================================================
  If you are a software developer, and this is your own package, then
  it is usually much more efficient to install the build dependencies
  in your environment once and use the following command that avoids
  a costly creation of a new virtual environment at every compilation:
  =====================================================================
   $ pip install 'scikit-build-core[pyproject]' setuptools_scm pybind11
   $ pip install --no-build-isolation -ve .
  =====================================================================
  You may optionally add -Ceditable.rebuild=true to auto-rebuild when
  the package is imported. Otherwise, you need to re-run the above
  after editing C++ files.")
endif()

if(NOT SKBUILD)
  # Manually detect the installed pybind11 package and import it into CMake.
  execute_process(
    COMMAND "${Python_EXECUTABLE}" -m pybind11 --cmakedir
    OUTPUT_STRIP_TRAILING_WHITESPACE
    OUTPUT_VARIABLE pybind11_DIR)
  list(APPEND CMAKE_PREFIX_PATH "${pybind11_DIR}")
endif()

# Import pybind11 through CMake's find_package mechanism
find_package(pybind11 CONFIG REQUIRED)

# We are now ready to compile the actual extension module
pybind11_add_module(py${PROJECT_NAME} bindings.cpp)
target_compile_definitions(py${PROJECT_NAME} PRIVATE Z3_FOUND)
target_link_libraries(py${PROJECT_NAME} PUBLIC ${PROJECT_NAME}_lib MQT::CoreECC MQT::CorePython)

# Install directive for scikit-build-core
install(TARGETS py${PROJECT_NAME} LIBRARY DESTINATION mqt/qecc)
