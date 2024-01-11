
find_package(Eigen3 3.4 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(HDF5 REQUIRED COMPONENTS C CXX)
include_directories(${HDF5_INCLUDE_DIRS})

# add_compile_definitions(CORTEO_6_BIT)

include_directories( ../extern/inifile-cpp/include )

add_executable(gendedx
    gendedx.cpp
    elements.h elements.cpp)

file(GLOB gendedx_files CONFIGURE_DEPENDS
    ${PROJECT_SOURCE_DIR}/dedx/*.cpp)

add_library(iondedx SHARED
    dedx.h dedx.cpp ${gendedx_files})

add_executable(gencorteo
    gencorteo.cpp
    corteo.h
    xs.h xs.cpp)

add_library(corteo4bit SHARED
    corteo4bit.h corteo4bitdata.cpp)

add_library(corteo6bit SHARED
    corteo6bit.h corteo6bitdata.cpp)

add_executable(iradina++  
    iradina.cpp 
    corteo.h corteo4bit.h corteo6bit.h
    xs.h xs.cpp
    elements.h elements.cpp
    ion.h ion.cpp
    random_vars.h random_vars.cpp
    target.h target.cpp
    simulation.h simulation.cpp

    geometry.h
    ion_beam.h ion_beam.cpp
    straggling.cpp
    arrays.h
    event_recorder.h event_recorder.cpp
    out_file.h out_file.cpp
    simulation_base.cpp)

target_link_libraries(iradina++ "-lm"
    Eigen3::Eigen
    iondedx
    corteo4bit
    corteo6bit
    ${HDF5_LIBRARIES})

