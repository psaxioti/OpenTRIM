cmake_minimum_required(VERSION 3.5.0)
    
# here we specify that the project is C++ language only, so the default
# C++ compiler on the system will be used
project(iradina++ CXX)

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

target_link_libraries(iradina++ "-lm" Eigen3::Eigen)

