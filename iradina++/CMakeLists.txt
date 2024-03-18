
find_package(Eigen3 3.4 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(HDF5 REQUIRED COMPONENTS C CXX)
include_directories(${HDF5_INCLUDE_DIRS})

# add_compile_definitions(CORTEO_6_BIT)

include_directories( ../extern/inifile-cpp/include )

include_directories( ../extern/json-3.11.3/single_include/)

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
    xs.h)

add_custom_command(
    OUTPUT  corteo4bitdata.cpp # corteo6bitdata.cpp
    COMMAND gencorteo
    DEPENDS gencorteo
)

add_library(corteo4bit SHARED
    ${CMAKE_CURRENT_BINARY_DIR}/corteo4bitdata.cpp)

#add_library(corteo6bit SHARED
#    corteo6bit.h
#    ${CMAKE_CURRENT_BINARY_DIR}/corteo6bitdata.cpp)

add_executable(iradina++  
    iradina.cpp 
    corteo.h
    xs.h
    elements.h elements.cpp
    ion.h ion.cpp
    random_vars.h
    target.h target.cpp
    simulation.h simulation.cpp
    geometry.h
    ion_beam.h ion_beam.cpp
    straggling.cpp
    arrays.h
    event_stream.h event_stream.cpp
    out_file.h out_file.cpp
    simulation_base.cpp
    tally.h
    options.h options.cpp
    parse_json.cpp
)

target_link_libraries(iradina++ "-lm"
    Eigen3::Eigen
    iondedx
    corteo4bit  # corteo6bit
    ${HDF5_LIBRARIES})



