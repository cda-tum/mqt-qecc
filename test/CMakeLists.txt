file(COPY ${PROJECT_SOURCE_DIR}/test/resources DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

package_add_test_with_working_dir(
  mqt-qecc-test
  MQT::QECC
  ${CMAKE_CURRENT_BINARY_DIR}
  test_ufdheuristic.cpp
  test_utils.cpp
  test_originalUfd.cpp
  test_simulator.cpp
  test_treenode.cpp)
