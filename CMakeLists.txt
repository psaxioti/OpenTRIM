cmake_minimum_required(VERSION 3.5.0)
    
project(iradina++ C CXX)

set (IRADINAPP_MAJOR_VERSION "0")
set (IRADINAPP_MINOR_VERSION "1")
set (IRADINAPP_BUGFIX_VERSION "4")
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

target_include_directories(gencorteo PRIVATE
    extern/cxxopts-3.1.1/include)

add_custom_command(OUTPUT xs_zbl_data.cpp COMMAND gencorteo -s zbl DEPENDS gencorteo)
add_library(xs_zbl SHARED ${CMAKE_CURRENT_BINARY_DIR}/xs_zbl_data.cpp)
add_custom_command(OUTPUT xs_lj_data.cpp COMMAND gencorteo -s lj DEPENDS gencorteo)
add_library(xs_lj SHARED ${CMAKE_CURRENT_BINARY_DIR}/xs_lj_data.cpp)
add_custom_command(OUTPUT xs_krc_data.cpp COMMAND gencorteo -s krc DEPENDS gencorteo)
add_library(xs_krc SHARED ${CMAKE_CURRENT_BINARY_DIR}/xs_krc_data.cpp)
add_custom_command(OUTPUT xs_moliere_data.cpp COMMAND gencorteo -s moliere DEPENDS gencorteo)
add_library(xs_moliere SHARED ${CMAKE_CURRENT_BINARY_DIR}/xs_moliere_data.cpp)

set (IRADINAPP_HEADERS
    src/dedx.h
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
    src/h5serialize.cpp
    src/tally.cpp
    src/parse_json.cpp    
    src/mcdriver.cpp)

set(IRADINAPP_LIBS
    "-lm"
    Eigen3::Eigen
    iondedx
    xs_zbl xs_lj xs_krc xs_moliere
    ${HDF5_LIBRARIES})

target_link_libraries(iradinapp ${IRADINAPP_LIBS})

target_include_directories(iradinapp PRIVATE
    extern/json-3.11.3/single_include/ )

add_executable(iradina++
    src/iradina.cpp
)

target_include_directories(iradina++ PRIVATE
    extern/cxxopts-3.1.1/include)

target_link_libraries(iradina++
    ${IRADINAPP_LIBS} iradinapp
    )

# add_subdirectory(doc)

set(installable_libs
    iondedx
    xs_zbl xs_lj xs_krc xs_moliere
    iradinapp )
install(TARGETS ${installable_libs} DESTINATION lib)
install(FILES ${IRADINAPP_HEADERS} DESTINATION include)
install(TARGETS iradina++ DESTINATION bin)


