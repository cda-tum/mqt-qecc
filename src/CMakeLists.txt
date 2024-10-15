# collect project headers
file(GLOB MQT_QECC_HEADERS CONFIGURE_DEPENDS ${MQT_QECC_INCLUDE_BUILD_DIR}/*.hpp)

# main project library
add_library(${PROJECT_NAME} ${MQT_QECC_HEADERS} DecodingSimulator.cpp UFDecoder.cpp UFHeuristic.cpp
                            Utils.cpp GF2.cpp)

# set include directories
target_include_directories(${PROJECT_NAME} PUBLIC ${MQT_QECC_INCLUDE_BUILD_DIR})

# add options and warnings to the library
target_link_libraries(
  ${PROJECT_NAME}
  PRIVATE MQT::ProjectWarnings MQT::ProjectOptions
  PUBLIC nlohmann_json::nlohmann_json)

# add MQT alias
add_library(MQT::QECC ALIAS ${PROJECT_NAME})

# add the ECC framework code
add_subdirectory(ecc)

if(BUILD_MQT_QECC_BINDINGS)
  add_subdirectory(python)
endif()
