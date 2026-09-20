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
        src/rsolverwave.cpp

        include/rconvection.h
        include/reigenvaluesolver.h
        include/rhemicube.h
        include/rhemicubepixel.h
        include/rhemicubesector.h
        include/riterationinfo.h
        include/riterationinfovalue.h
        include/rlocalrotation.h
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
        include/rsolverheat.h
        include/rsolvermagnetostatics.h
        include/rsolvermesh.h
        include/rsolverradiativeheat.h
        include/rsolvershareddata.h
        include/rsolverstress.h
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

# Unit tests (opt-in via -DRANGE_BUILD_TESTS=ON)
if(RANGE_BUILD_TESTS)
    add_subdirectory(tests)
endif()
