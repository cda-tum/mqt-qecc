# main project library
add_library(
  ${PROJECT_NAME}
  ${MQT_QECC_INCLUDE_BUILD_DIR}/Code.hpp
  ${MQT_QECC_INCLUDE_BUILD_DIR}/Codes.hpp
  ${MQT_QECC_INCLUDE_BUILD_DIR}/Decoder.hpp
  ${MQT_QECC_INCLUDE_BUILD_DIR}/DecodingRunInformation.hpp
  ${MQT_QECC_INCLUDE_BUILD_DIR}/DecodingSimulator.hpp
  ${MQT_QECC_INCLUDE_BUILD_DIR}/QeccException.hpp
  ${MQT_QECC_INCLUDE_BUILD_DIR}/TreeNode.hpp
  ${MQT_QECC_INCLUDE_BUILD_DIR}/UFDecoder.hpp
  ${MQT_QECC_INCLUDE_BUILD_DIR}/UFHeuristic.hpp
  ${MQT_QECC_INCLUDE_BUILD_DIR}/Utils.hpp
  DecodingSimulator.cpp
  UFDecoder.cpp
  UFHeuristic.cpp
  Utils.cpp)

# set include directories
target_include_directories(${PROJECT_NAME} PUBLIC ${MQT_QECC_INCLUDE_BUILD_DIR})

# add options and warnings to the library
target_link_libraries(
  ${PROJECT_NAME}
  PRIVATE MQT::ProjectWarnings MQT::ProjectOptions
  PUBLIC nlohmann_json::nlohmann_json ldpc)

# add MQT alias
add_library(MQT::QECC ALIAS ${PROJECT_NAME})

# add the ECC framework code
add_subdirectory(ecc)

if(BUILD_MQT_QECC_BINDINGS)
  add_subdirectory(python)
endif()
