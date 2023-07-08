# add MQT::QFR library
add_subdirectory("${PROJECT_SOURCE_DIR}/extern/mqt-core" "extern/mqt-core"
                 EXCLUDE_FROM_ALL)
# main project library
add_library(
  ${PROJECT_NAME}_lib
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
target_include_directories(
  ${PROJECT_NAME}_lib PUBLIC ${PROJECT_SOURCE_DIR}/include
                             ${PROJECT_BINARY_DIR}/include)

# add options and warnings to the library
target_link_libraries(${PROJECT_NAME}_lib PUBLIC project_options
                                                 project_warnings)

target_link_libraries(${PROJECT_NAME}_lib PUBLIC nlohmann_json)

find_package(FLINT REQUIRED)
target_link_libraries(${PROJECT_NAME}_lib PUBLIC flint)

# add MQT alias
add_library(MQT::${PROJECT_NAME}_lib ALIAS ${PROJECT_NAME}_lib)
