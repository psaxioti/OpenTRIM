cmake_minimum_required(VERSION 3.5.0)
    
project(iradina++ C CXX)

set (IRADINAPP_MAJOR_VERSION "0")
set (IRADINAPP_MINOR_VERSION "1")
set (IRADINAPP_BUGFIX_VERSION "0")
set (IRADINAPP_VERSION_STRING ${IRADINAPP_MAJOR_VERSION}.${IRADINAPP_MINOR_VERSION}.${IRADINAPP_BUGFIX_VERSION})
message (STATUS "iradina++ version ${IRADINAPP_VERSION_STRING}")

include (GNUInstallDirs)

# by default, install into $HOME/.local (not /usr/local), so that no root access (and sudo!!) is needed
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  if((CMAKE_SYSTEM_NAME STREQUAL "Windows") AND (NOT CMAKE_CROSSCOMPILING))
    set(CMAKE_INSTALL_PREFIX "$ENV{USERPROFILE}/iradina++" CACHE PATH "Default install path" FORCE)
  else()
    set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/.local" CACHE PATH "Default install path" FORCE)
  endif()
endif()

message (STATUS "Install prefix =  ${CMAKE_INSTALL_PREFIX}") 

# use install rpath, no need to use LD_LIBRARY_PATH 
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH OFF)

find_package(Eigen3 3.4 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})
message(STATUS "Eigen3 include ${EIGEN3_INCLUDE_DIR}")

find_package(HDF5 REQUIRED COMPONENTS C CXX)
include_directories(${HDF5_INCLUDE_DIRS})
message(STATUS "HDF5 include: ${HDF5_INCLUDE_DIRS}")
message(STATUS "HDF5 libs: ${HDF5_LIBRARIES}")

add_executable(gendedx
    src/gendedx.cpp
    src/elements.h src/elements.cpp)

file(GLOB gendedx_files CONFIGURE_DEPENDS
    ${PROJECT_SOURCE_DIR}/dedx/*.cpp)

add_library(iondedx SHARED
    src/dedx.h src/dedx.cpp 
    ${gendedx_files})

add_executable(gencorteo
    src/gencorteo.cpp
    src/corteo.h
    src/xs.h)

add_custom_command(
    OUTPUT  corteo4bitdata.cpp corteo6bitdata.cpp
    COMMAND gencorteo
    DEPENDS gencorteo
)

add_library(corteo4bit SHARED
    ${CMAKE_CURRENT_BINARY_DIR}/corteo4bitdata.cpp)

add_library(corteo6bit SHARED
    ${CMAKE_CURRENT_BINARY_DIR}/corteo6bitdata.cpp)

set (IRADINAPP_HEADERS
    src/corteo.h
    src/xs.h
    src/elements.h
    src/ion.h 
    src/random_vars.h
    src/target.h 
    src/mccore.h 
    src/geometry.h
    src/ion_beam.h 
    src/arrays.h
    src/event_stream.h 
    src/out_file.h
    src/tally.h 
    src/mcdriver.h 
  )

add_library(iradinapp SHARED
    ${IRADINAPP_HEADERS}
    src/elements.cpp
    src/ion.cpp   
    src/target.cpp
    src/mccore.cpp   
    src/ion_beam.cpp
    src/straggling.cpp    
    src/event_stream.cpp
    src/out_file.cpp
    src/tally.cpp
    src/parse_json.cpp    
    src/mcdriver.cpp)

add_executable(iradina++  
    src/iradina.cpp 
)

target_link_libraries(iradinapp "-lm"
    Eigen3::Eigen
    iondedx
    corteo4bit  corteo6bit
    ${HDF5_LIBRARIES})

set(IRADINAPP_LIBS
    "-lm"
    Eigen3::Eigen
    iondedx
    corteo4bit  corteo6bit
    iradinapp
    ${HDF5_LIBRARIES})

target_link_libraries(iradina++ 
    ${IRADINAPP_LIBS})

target_include_directories(iradinapp PRIVATE
    extern/json-3.11.3/single_include/ )

# add_subdirectory(doc)

set(installable_libs
    iondedx
    corteo4bit
    corteo6bit
    iradinapp )
install(TARGETS ${installable_libs} DESTINATION lib)
install(FILES ${IRADINAPP_HEADERS} DESTINATION include)
install(TARGETS iradina++ DESTINATION bin)


