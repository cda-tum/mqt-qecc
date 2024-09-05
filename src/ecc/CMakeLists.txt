if(NOT TARGET ${PROJECT_NAME}-ecc)
  file(GLOB_RECURSE ECC_HEADERS ${MQT_QECC_INCLUDE_BUILD_DIR}/ecc/*.hpp)

  # add ECC package library
  add_library(
    ${PROJECT_NAME}-ecc
    ${ECC_HEADERS}
    Ecc.cpp
    Q3Shor.cpp
    Q5Laflamme.cpp
    Q7Steane.cpp
    Q9Shor.cpp
    Q9Surface.cpp
    Q18Surface.cpp)

  target_link_libraries(
    ${PROJECT_NAME}-ecc
    PUBLIC MQT::CoreIR
    PRIVATE MQT::ProjectOptions MQT::ProjectWarnings)

  target_include_directories(${PROJECT_NAME}-ecc PUBLIC ${MQT_QECC_INCLUDE_BUILD_DIR})

  # add MQT alias
  add_library(MQT::QECCFramework ALIAS ${PROJECT_NAME}-ecc)
endif()
