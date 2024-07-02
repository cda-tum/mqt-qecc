# main project library
add_library(
  ${PROJECT_NAME}
  ${PROJECT_SOURCE_DIR}/include/Code.hpp
  ${PROJECT_SOURCE_DIR}/include/Codes.hpp
  ${PROJECT_SOURCE_DIR}/include/Decoder.hpp
  ${PROJECT_SOURCE_DIR}/include/DecodingRunInformation.hpp
  ${PROJECT_SOURCE_DIR}/include/DecodingSimulator.hpp
  ${PROJECT_SOURCE_DIR}/include/QeccException.hpp
  ${PROJECT_SOURCE_DIR}/include/TreeNode.hpp
  ${PROJECT_SOURCE_DIR}/include/UFDecoder.hpp
  ${PROJECT_SOURCE_DIR}/include/UFHeuristic.hpp
  ${PROJECT_SOURCE_DIR}/include/Utils.hpp
  DecodingSimulator.cpp
  UFDecoder.cpp
  UFHeuristic.cpp)

# set include directories
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/include
                                                  ${PROJECT_BINARY_DIR}/include)

# add options and warnings to the library
target_link_libraries(
  ${PROJECT_NAME}
  PRIVATE MQT::ProjectWarnings MQT::ProjectOptions
  PUBLIC nlohmann_json flint)

# add MQT alias
add_library(MQT::QECC ALIAS ${PROJECT_NAME})

if(BUILD_MQT_QECC_BINDINGS)
  add_subdirectory(python)
endif()
