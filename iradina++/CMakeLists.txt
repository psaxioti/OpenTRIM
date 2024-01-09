
find_package(Eigen3 3.4 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(HDF5 REQUIRED COMPONENTS C CXX)
include_directories(${HDF5_INCLUDE_DIRS})

# add_compile_definitions(CORTEO_6_BIT)

include_directories( ../extern/inifile-cpp/include )

add_executable(iradina++  
    iradina.cpp 
    corteo.h
    xs.h xs.cpp
    elements.h elements.cpp
    ion.h ion.cpp
    random_vars.h random_vars.cpp
    target.h target.cpp
    simulation.h simulation.cpp
    urbg.h
    geometry.h
    ion_beam.h ion_beam.cpp
    straggling.cpp
    arrays.h
    event_recorder.h event_recorder.cpp
    out_file.h out_file.cpp
    gencorteo.cpp)



add_executable(gendedx
    gendedx.cpp elements.cpp)

file(GLOB gendedx_files CONFIGURE_DEPENDS
    ${PROJECT_SOURCE_DIR}/dedx/*.cpp)

add_library(iondedx SHARED
    dedx.h dedx.cpp ${gendedx_files})

target_link_libraries(iradina++ "-lm"
    Eigen3::Eigen iondedx ${HDF5_LIBRARIES})

