cmake_minimum_required(VERSION 3.5.0)
    
# here we specify that the project is C language only, so the default
# C compiler on the system will be used
project(iradina C)

set(CMAKE_C_FLAGS "-O1 -Wall")

add_executable(iradina 
    fileio.c fileio.h 
    fromcorteo.c fromcorteo.h 
    geometry.c geometry.h 
    iradina.c iradina.h 
    target.c target.h 
    transport.c transport.h 
    utils.c utils.h)

target_link_libraries(iradina "-lm")