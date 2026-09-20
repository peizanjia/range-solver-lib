# Unit tests for range-solver-lib (QTest based).
#
# Enabled by configuring the project with -DRANGE_BUILD_TESTS=ON.
# Run with `ctest` from the build directory.

find_package(Qt6 REQUIRED COMPONENTS Test)

set(CMAKE_AUTOMOC ON)

set(range_solver_lib_tests
    tst_scales
    tst_local_rotation
    tst_iteration_info_value
    tst_matrix_solver
    tst_solver_acoustic
    tst_solver_stress
    tst_structural_solver_gate
    tst_eigen_value_solver
)

foreach(test_name IN LISTS range_solver_lib_tests)
    qt_add_executable(${test_name} ${test_name}.cpp)
    target_link_libraries(${test_name}
        PRIVATE
            range-solver-lib
            Qt6::Test
    )
    add_test(NAME ${test_name} COMMAND ${test_name})
endforeach()
