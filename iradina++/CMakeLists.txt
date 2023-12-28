
find_package(Eigen3 3.4 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# add_compile_definitions(CORTEO_6_BIT)

add_executable(iradina++  
    iradina.cpp 
    corteo.h corteo.cpp
    xs.h xs.cpp
    elements.h elements.cpp
    ion.h ion.cpp
    random_vars.h random_vars.cpp
    dedx.h dedx.cpp)

add_executable(gendedx
    gendedx.cpp elements.cpp
    target.h target.cpp)

file(GLOB gendedx_files CONFIGURE_DEPENDS
    ${PROJECT_SOURCE_DIR}/dedx/*.cpp)

add_library(iondedx SHARED
    dedx.h dedx.cpp ${gendedx_files})

target_link_libraries(iradina++ "-lm" Eigen3::Eigen iondedx)

