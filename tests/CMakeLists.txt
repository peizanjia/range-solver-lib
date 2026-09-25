# Unit tests for range-solver-lib (QTest based).
#
# Enabled by configuring the project with -DRANGE_BUILD_TESTS=ON.
# Run with `ctest` from the build directory.

find_package(Qt6 REQUIRED COMPONENTS Test)

if(WIN32 AND CMAKE_VERSION VERSION_LESS 3.22)
    message(FATAL_ERROR "Windows QTest execution requires CMake 3.22 or newer for ENVIRONMENT_MODIFICATION")
endif()

set(CMAKE_AUTOMOC ON)

set(range_solver_lib_tests
    tst_scales
    tst_local_rotation
    tst_iteration_info_value
    tst_matrix_solver
    tst_equation_book
    tst_solver_acoustic
    tst_solver_stress
    tst_structural_solver_gate
    tst_eigen_value_solver
    tst_modal_damping
    tst_parallel_execution_context
    tst_structural_dynamics_integrator
    tst_structural_cable_kernel
    tst_structural_cable_equilibrium
    tst_structural_cable_fe
    tst_structural_piezo_kernel
    tst_structural_piezo_fe
    tst_structural_control_force_map
    tst_structural_control_conditioner
    tst_structural_dynamics_fe
    tst_structural_dynamics_beam
    tst_structural_dynamics_loads
    tst_structural_dynamics_support
    tst_structural_dynamics_restart_io
    tst_structural_dynamics_failure_atomicity
    tst_structural_dynamics_rotational_initial
    tst_structural_dynamics_parallel_repro
    tst_structural_temperature_field
    tst_structural_heat_certification
    tst_structural_reference_thermal_mass
    tst_structural_cable_dynamics
    tst_structural_cable_restart
    tst_structural_piezo_coupled_1d
    tst_structural_heat_stress_coupling
    tst_structural_monitor_output_concurrency
    tst_structural_thermal_bending
    tst_structural_piezo_network_fe
    tst_structural_thermal_control
    tst_structural_thermo_coupled_1d
    tst_structural_thermo_coupled_fe_adapter
    tst_structural_thermo_monolithic_1d
)

foreach(test_name IN LISTS range_solver_lib_tests)
    qt_add_executable(${test_name} ${test_name}.cpp)
    target_link_libraries(${test_name}
        PRIVATE
            range-solver-lib
            Qt6::Test
    )
    add_test(NAME ${test_name} COMMAND ${test_name})
    if(WIN32)
        set_tests_properties(${test_name} PROPERTIES
            ENVIRONMENT_MODIFICATION "PATH=path_list_prepend:$<TARGET_FILE_DIR:Qt6::Core>")
    endif()
endforeach()
