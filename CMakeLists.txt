cmake_minimum_required(VERSION 3.12.0)
    
project(OpenTRIM 
        VERSION 0.4.0
        DESCRIPTION "Monte-Carlo simulation of ion transport in materials"
        HOMEPAGE_URL https://gitlab.com/ir2-lab/
        LANGUAGES C CXX)

set(PROJECT_NAME_LOWERCASE opentrim)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

option (BUILD_GUI "build gui" ON)

message (STATUS "${PROJECT_NAME} version ${PROJECT_VERSION}")

string(TIMESTAMP CMAKE_BUILD_TIME UTC)

message (STATUS "Build time ${CMAKE_BUILD_TIME}")

add_compile_options("-march=native")

add_compile_definitions(PROJECT_NAME="${PROJECT_NAME}"
                        PROJECT_VERSION="${PROJECT_VERSION}"
                        PROJECT_DESCRIPTION="${PROJECT_DESCRIPTION}"
                        SYSTEM_ID="${CMAKE_CXX_PLATFORM_ID}"
                        COMPILER_ID="${CMAKE_CXX_COMPILER_ID}"
                        COMPILER_VERSION="${CMAKE_CXX_COMPILER_VERSION}"
                        BUILD_TIME="${CMAKE_BUILD_TIME}")

include (GNUInstallDirs) 

# by default, install into $HOME/.local (not /usr/local), so that no root access (and sudo!!) is needed
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  if((CMAKE_SYSTEM_NAME STREQUAL "Windows") AND (NOT CMAKE_CROSSCOMPILING))
    set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/${PROJECT_NAME}" CACHE PATH "Default install path" FORCE)
  else()
    set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/.local" CACHE PATH "Default install path" FORCE)
  endif()
endif()

message (STATUS "Install prefix =  ${CMAKE_INSTALL_PREFIX}") 

# use install rpath, no need to use LD_LIBRARY_PATH 
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_LIBDIR}")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH OFF)

find_package(Eigen3 3.4 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})
message(STATUS "Eigen3 include ${EIGEN3_INCLUDE_DIR}")

find_package(HDF5 REQUIRED COMPONENTS C)
include_directories(${HDF5_INCLUDE_DIRS})
message(STATUS "HDF5 include: ${HDF5_INCLUDE_DIRS}")
message(STATUS "HDF5 libs: ${HDF5_LIBRARIES}")

include_directories(src)

include(Externals)
include(HelperFunctions)

file(GLOB gendedx_files CONFIGURE_DEPENDS
    ${PROJECT_SOURCE_DIR}/dedx/*.cpp)

add_library(dedx SHARED
    src/dedx.cpp
    src/straggling.cpp
    ${gendedx_files}
)
target_sources(dedx
    PUBLIC FILE_SET HEADERS
        BASE_DIRS src
        FILES 
            src/corteo.h
            src/dedx.h
)

add_dedx_library(zbl)
add_dedx_library(lj)
add_dedx_library(krc)
add_dedx_library(moliere)

set (LIB_TARGET ${PROJECT_NAME_LOWERCASE})
set (EXE_TARGET ${PROJECT_NAME_LOWERCASE}_exe)

add_library(${LIB_TARGET} SHARED
    ${PTABLE_DATA_CPP}
    src/ion.cpp   
    src/target.cpp
    src/flight_path_calc.cpp
    src/mccore.cpp   
    src/ion_beam.cpp
    src/event_stream.cpp
    src/h5serialize.cpp
    src/tally.cpp
    src/parse_json.cpp    
    src/mcdriver.cpp
    src/dedx_calc.cpp
)
target_sources(${LIB_TARGET}
    PUBLIC FILE_SET HEADERS 
        BASE_DIRS src
        FILES
            src/corteo_xs.h
            src/xs.h
            src/periodic_table.h
            src/ion.h 
            src/cascade_queue.h
            src/random_vars.h
            src/target.h 
            src/flight_path_calc.h
            src/dedx_calc.h
            src/mccore.h 
            src/geometry.h
            src/ion_beam.h 
            src/arrays.h
            src/event_stream.h 
            src/tally.h 
            src/mcdriver.h 
)

set(PROJECT_LIBS
    "-lm"
    ${HDF5_LIBRARIES}
    dedx
    xs_zbl xs_lj xs_krc xs_moliere
)

target_link_libraries(${LIB_TARGET} ${PROJECT_LIBS})

FetchContent_MakeAvailable(external_json)
FetchContent_MakeAvailable(external_highfive)
target_include_directories(${LIB_TARGET} PRIVATE
    ${external_json_SOURCE_DIR}/single_include/ 
    ${external_highfive_SOURCE_DIR}/include)

add_executable(${EXE_TARGET}
    src/main.cpp
)
# adjust exe name to be the same as the library's
set_property(TARGET ${EXE_TARGET} PROPERTY OUTPUT_NAME ${PROJECT_NAME_LOWERCASE})

FetchContent_MakeAvailable(external_cxxopts)
target_include_directories(${EXE_TARGET} PRIVATE
    ${external_cxxopts_SOURCE_DIR}/include
    extern/progressbar/include)

target_link_libraries(${EXE_TARGET}
    ${LIB_TARGET}
)

# add_subdirectory(doc)

set(installable_libs
    dedx
    xs_zbl xs_lj xs_krc xs_moliere
    ${LIB_TARGET} )

install(TARGETS ${installable_libs} ${EXE_TARGET}
    LIBRARY
    FILE_SET HEADERS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME_LOWERCASE}
    RUNTIME
)

if (BUILD_GUI)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

find_package(QT NAMES Qt6 Qt5 REQUIRED COMPONENTS Widgets)
find_package(Qt${QT_VERSION_MAJOR} REQUIRED COMPONENTS Widgets)
message(STATUS "Found Qt${QT_VERSION_MAJOR}")

find_package(Qwt REQUIRED)
message(STATUS "Found Qwt ${QWT_VERSION_STRING} at ${QWT_INCLUDE_DIR}")

FetchContent_MakeAvailable(external_qmatplotwidget)

set(GUI_SOURCES
    src/gui/main.cpp
    src/gui/mainui.cpp
    src/gui/opentrim.qrc
    src/gui/simulationoptionsview.cpp
    src/gui/periodictablewidget.cpp
    src/gui/materialsdefview.cpp
    src/gui/optionsmodel.cpp
    src/gui/floatlineedit.cpp
    src/gui/mydatawidgetmapper.cpp
    src/gui/regionsview.cpp
    src/gui/welcomeview.cpp
    src/gui/runview.cpp
    src/gui/mcdriverobj.cpp
    src/gui/resultsview.cpp
    src/gui/simcontrolwidget.cpp
    ${external_qmatplotwidget_SOURCE_DIR}/src/QMatPlotWidget.cpp
    src/gui/jsedit/jsedit.cpp
)

set (GUI_TARGET ${PROJECT_NAME_LOWERCASE}-gui)

if(${QT_VERSION_MAJOR} GREATER_EQUAL 6)
    qt_add_executable(${GUI_TARGET}
        MANUAL_FINALIZATION
        ${GUI_SOURCES}
    )
else()
    add_executable(${GUI_TARGET}
        ${GUI_SOURCES}
    )
endif()

FetchContent_MakeAvailable(external_json)
target_include_directories(${GUI_TARGET} PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${external_json_SOURCE_DIR}/single_include/
    ${external_qmatplotwidget_SOURCE_DIR}/src
)

target_link_libraries(${GUI_TARGET} PRIVATE
    Qt${QT_VERSION_MAJOR}::Widgets
    Qwt::Qwt
    ${LIB_TARGET}
)

# Qt for iOS sets MACOSX_BUNDLE_GUI_IDENTIFIER automatically since Qt 6.1.
# If you are developing for iOS or macOS you should consider setting an
# explicit, fixed bundle identifier manually though.
if(${QT_VERSION} VERSION_LESS 6.1.0)
  set(BUNDLE_ID_OPTION MACOSX_BUNDLE_GUI_IDENTIFIER com.example.ionsui)
endif()

set_target_properties(${GUI_TARGET} PROPERTIES
    ${BUNDLE_ID_OPTION}
    MACOSX_BUNDLE_BUNDLE_VERSION ${PROJECT_VERSION}
    MACOSX_BUNDLE_SHORT_VERSION_STRING ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
    MACOSX_BUNDLE TRUE
    WIN32_EXECUTABLE TRUE
)

install(TARGETS ${GUI_TARGET}
    BUNDLE DESTINATION .
    LIBRARY
    RUNTIME
)

message(STATUS "Install ${GUI_TARGET} at ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}")

if(QT_VERSION_MAJOR EQUAL 6)
    qt_finalize_executable(${GUI_TARGET})
endif()


endif () ## BUILD GUI





