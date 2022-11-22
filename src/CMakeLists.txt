# main project library
add_library(
  ${PROJECT_NAME}_lib
  ${CMAKE_CURRENT_SOURCE_DIR}/UFHeuristic.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/UFDecoder.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/DecodingSimulator.cpp
  ${${PROJECT_NAME}_SOURCE_DIR}/include/Decoder.hpp
  ${${PROJECT_NAME}_SOURCE_DIR}/include/DecodingSimulator.hpp
  ${${PROJECT_NAME}_SOURCE_DIR}/include/Code.hpp
  ${${PROJECT_NAME}_SOURCE_DIR}/include/TreeNode.hpp
  ${${PROJECT_NAME}_SOURCE_DIR}/include/Codes.hpp
  ${${PROJECT_NAME}_SOURCE_DIR}/include/Utils.hpp
  ${${PROJECT_NAME}_SOURCE_DIR}/include/UFHeuristic.hpp
  ${${PROJECT_NAME}_SOURCE_DIR}/include/UFDecoder.hpp
  ${${PROJECT_NAME}_SOURCE_DIR}/include/DecodingRunInformation.hpp
  ${${PROJECT_NAME}_SOURCE_DIR}/include/QeccException.hpp)

# set include directories
target_include_directories(${PROJECT_NAME}_lib
                           PUBLIC $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/include>)

# set required C++ standard and disable compiler specific extensions
target_compile_features(${PROJECT_NAME}_lib PUBLIC cxx_std_20)
set_target_properties(${PROJECT_NAME}_lib PROPERTIES CMAKE_CXX_STANDARD_REQUIRED ON)
set_target_properties(${PROJECT_NAME}_lib PROPERTIES CXX_EXTENSIONS OFF)

add_subdirectory("${PROJECT_SOURCE_DIR}/extern/qfr" "extern/qfr")
target_link_libraries(${PROJECT_NAME}_lib PUBLIC MQT::qfr)

# enable interprocedural optimization if it is supported
include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported)
if(ipo_supported)
  set_target_properties(${PROJECT_NAME}_lib PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

# set compiler flags depending on compiler target_compile_options(${PROJECT_NAME}_lib PUBLIC
# -fsanitize=address -fsanitize=undefined)
if(MSVC)
  target_compile_options(${PROJECT_NAME}_lib PUBLIC /utf-8)
else()
  target_compile_options(${PROJECT_NAME}_lib PUBLIC -Wall -Wextra $<$<CONFIG:DEBUG>:-Og>)
  if(BINDINGS)
    # adjust visibility settings for building Python bindings
    target_compile_options(${PROJECT_NAME}_lib PUBLIC -fvisibility=hidden)
  endif()
  if(NOT DEPLOY)
    # only include machine-specific optimizations when building for the host machine
    target_compile_options(${PROJECT_NAME}_lib PUBLIC -mtune=native)
    include(CheckCXXCompilerFlag)
    check_cxx_compiler_flag(-march=native HAS_MARCH_NATIVE)
    if(HAS_MARCH_NATIVE)
      target_compile_options(${PROJECT_NAME}_lib PUBLIC -march=native)
    endif()
  endif()
endif()

if(GENERATE_POSITION_INDEPENDENT_CODE OR BINDINGS)
  include(CheckPIESupported)
  check_pie_supported()
  set_target_properties(${PROJECT_NAME}_lib PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
endif()

# add nlohmann::json library
set(JSON_BuildTests # cmake-lint: disable=C0103
    OFF
    CACHE INTERNAL "")
set(JSON_MultipleHeaders # cmake-lint: disable=C0103
    OFF
    CACHE INTERNAL "")
#add_subdirectory("${PROJECT_SOURCE_DIR}/extern/qfr/extern/json" "extern/qfr/extern/json" EXCLUDE_FROM_ALL)
target_link_libraries(${PROJECT_NAME}_lib PUBLIC nlohmann_json)

# add coverage compiler and linker flag to the library and all targets that link against it, if
# COVERAGE is set
if(COVERAGE)
  target_compile_options(${PROJECT_NAME}_lib PUBLIC --coverage)
  target_link_libraries(${PROJECT_NAME}_lib PUBLIC --coverage)
endif()

find_package(FLINT REQUIRED)
target_link_libraries(${PROJECT_NAME}_lib PUBLIC flint)

# add MQT alias
add_library(MQT::${PROJECT_NAME}_lib ALIAS ${PROJECT_NAME}_lib)
