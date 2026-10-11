# Create a static library
qt_add_library(range-solver-lib
    STATIC
        src/rconvection.cpp
        src/reigenvaluesolver.cpp
        src/rhemicube.cpp
        src/rhemicubepixel.cpp
        src/rhemicubesector.cpp
        src/riterationinfo.cpp
        src/riterationinfovalue.cpp
        src/rlocalrotation.cpp
        src/rmodal_damping.cpp
        src/rstructural_linear_model.cpp
        src/rparallel_execution_context.cpp
        src/rmatrixmanager.cpp
        src/rmatrixpreconditioner.cpp
        src/rmatrixsolver.cpp
        src/requationbook.cpp
        src/rscales.cpp
        src/rsolver.cpp
        src/rsolveracoustic.cpp
        src/rsolverelectrostatics.cpp
        src/rsolverfluid.cpp
        src/rsolverfluidheat.cpp
        src/rsolverfluidparticle.cpp
        src/rsolvergeneric.cpp
        src/rsolverheat.cpp
        src/rsolvermagnetostatics.cpp
        src/rsolvermesh.cpp
        src/rsolverradiativeheat.cpp
        src/rsolvershareddata.cpp
        src/rsolverstress.cpp
        src/rstructural_dynamics_integrator.cpp
        src/rstructural_cable_kernel.cpp
        src/rstructural_cable_control_input.cpp
        src/rstructural_cable_equilibrium.cpp
        src/rstructural_cable_dynamics_restart_adapter.cpp
        src/rstructural_cable_model_adapter.cpp
        src/rstructural_static_truss_model.cpp
        src/rstructural_truss_assembly.cpp
        src/rstructural_truss_dynamics.cpp
        src/rstructural_truss_dynamics_model_adapter.cpp
        src/rstructural_truss_dynamics_restart_adapter.cpp
        src/rstructural_truss_native_run.cpp
        src/rstructural_piezo_model_adapter.cpp
        src/rstructural_piezo_kernel.cpp
        src/rstructural_control_force_map.cpp
        src/rstructural_control_conditioner.cpp
        src/rstructural_control_runtime_parameters.cpp
        src/rstructural_typed_control.cpp
        src/rstructural_native_control.cpp
        src/rstructural_typed_fe_control.cpp
        src/rstructural_dynamics_support.cpp
        src/rstructural_temperature_field.cpp
        src/rstructural_coupled_field_snapshot.cpp
        src/rstructural_cable_dynamics.cpp
        src/rstructural_piezo_coupled_1d.cpp
        src/rstructural_piezo_beam_fe_adapter.cpp
        src/rstructural_thermal_control.cpp
        src/rstructural_thermal_fe_adapter.cpp
        src/rstructural_thermo_coupled_1d.cpp
        src/rstructural_thermo_coupled_fe_adapter.cpp
        src/rstructural_thermo_monolithic_1d.cpp
        src/rstructural_thermo_monolithic_1d_fe_adapter.cpp
        src/rsolverwave.cpp

        include/rconvection.h
        include/reigenvaluesolver.h
        include/rhemicube.h
        include/rhemicubepixel.h
        include/rhemicubesector.h
        include/riterationinfo.h
        include/riterationinfovalue.h
        include/rlocalrotation.h
        include/rmodal_damping.h
        include/rstructural_linear_model.h
        include/rparallel_execution_context.h
        include/rmatrixmanager.h
        include/rmatrixpreconditioner.h
        include/rmatrixsolver.h
        include/requationbook.h
        include/rscales.h
        include/rsolver.h
        include/rsolveracoustic.h
        include/rsolverelectrostatics.h
        include/rsolverfluid.h
        include/rsolverfluidheat.h
        include/rsolverfluidparticle.h
        include/rsolvergeneric.h
        include/raccepted_step_publication.h
        include/rsolverheat.h
        include/rsolvermagnetostatics.h
        include/rsolvermesh.h
        include/rsolverradiativeheat.h
        include/rsolvershareddata.h
        include/rsolverstress.h
        include/rstructural_dynamics_integrator.h
        include/rstructural_cable_kernel.h
        include/rstructural_cable_control_input.h
        include/rstructural_cable_equilibrium.h
        include/rstructural_cable_dynamics_restart_adapter.h
        include/rstructural_cable_model_adapter.h
        include/rstructural_static_truss_model.h
        include/rstructural_truss_assembly.h
        include/rstructural_truss_dynamics.h
        include/rstructural_truss_dynamics_model_adapter.h
        include/rstructural_truss_dynamics_restart_adapter.h
        include/rstructural_truss_native_run.h
        include/rstructural_piezo_model_adapter.h
        include/rstructural_piezo_kernel.h
        include/rstructural_control_force_map.h
        include/rstructural_control_conditioner.h
        include/rstructural_control_runtime_parameters.h
        include/rstructural_typed_control.h
        include/rstructural_native_control.h
        include/rstructural_typed_fe_control.h
        include/rstructural_control_protocol.h
        include/rstructural_dynamics_support.h
        include/rstructural_temperature_field.h
        include/rstructural_coupled_field_snapshot.h
        include/rstructural_cable_dynamics.h
        include/rstructural_piezo_coupled_1d.h
        include/rstructural_piezo_beam_fe_adapter.h
        include/rstructural_thermal_control.h
        include/rstructural_thermal_fe_adapter.h
        include/rstructural_thermo_coupled_1d.h
        include/rstructural_thermo_coupled_fe_adapter.h
        include/rstructural_thermo_monolithic_1d.h
        include/rstructural_thermo_monolithic_1d_fe_adapter.h
        include/rsolverwave.h
)

# Include directories
target_include_directories(range-solver-lib
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(range-solver-lib
    PUBLIC
        range-model-lib
        common_defines
)

# The persisted built-in provider identity describes the actual compiled source.
# Reconfigure when an input changes so an old identity cannot survive a rebuild.
set(range_native_control_hash_inputs
    "${CMAKE_CURRENT_SOURCE_DIR}/src/rstructural_native_control.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/rstructural_native_control.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../range-model-lib/include/rml_structural_native_control.h")
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${range_native_control_hash_inputs})
set(range_native_control_identity_material "range.native.pi.implementation/v1\n")
foreach(range_native_control_hash_input IN LISTS range_native_control_hash_inputs)
    file(SHA256 "${range_native_control_hash_input}" range_native_control_input_sha256)
    get_filename_component(range_native_control_input_name "${range_native_control_hash_input}" NAME)
    string(APPEND range_native_control_identity_material
        "${range_native_control_input_name}\n${range_native_control_input_sha256}\n")
endforeach()
string(SHA256 range_native_control_implementation_sha256 "${range_native_control_identity_material}")
target_compile_definitions(range-solver-lib PRIVATE
    RANGE_NATIVE_CONTROL_IMPLEMENTATION_SHA256="${range_native_control_implementation_sha256}")

# Unit tests (opt-in via -DRANGE_BUILD_TESTS=ON)
if(RANGE_BUILD_TESTS)
    add_subdirectory(tests)
endif()
