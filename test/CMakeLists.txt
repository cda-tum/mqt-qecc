if (NOT TARGET gtest OR NOT TARGET gmock)
    # Prevent overriding the parent project's compiler/linker settings on Windows
    set(gtest_force_shared_crt # cmake-lint: disable=C0103
            ON
            CACHE BOOL "" FORCE)
    add_subdirectory("${PROJECT_SOURCE_DIR}/extern/mqt-core/extern/googletest"
            "extern/mqt-core/extern/googletest" EXCLUDE_FROM_ALL)
    set_target_properties(gtest gtest_main gmock gmock_main PROPERTIES FOLDER extern)
endif ()

# create executable for /test/app.cpp
add_executable(${PROJECT_NAME}_app ${CMAKE_CURRENT_SOURCE_DIR}/app.cpp)
target_link_libraries(${PROJECT_NAME}_app PRIVATE ${PROJECT_NAME}_lib)
set_target_properties(
        ${PROJECT_NAME}_app
        PROPERTIES FOLDER tests
        CMAKE_CXX_STANDARD_REQUIRED ON
        CXX_EXTENSIONS OFF)

package_add_test(
        mqt-qecc-test
        MQT::QECC
        test_ufdheuristic.cpp
        test_utils.cpp
        test_originalUfd.cpp
        test_simulator.cpp
        test_treenode.cpp)

# files containing codes for tests
file(
        COPY "${PROJECT_SOURCE_DIR}/examples/test/"
        DESTINATION "${CMAKE_CURRENT_SOURCE_DIR}/resources/codes/"
        PATTERN "*")
