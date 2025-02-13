include(FetchContent)

set(FETCHCONTENT_BASE_DIR ${CMAKE_BINARY_DIR}/external)
set(FETCHCONTENT_QUIET FALSE)
set(FETCHCONTENT_TRY_FIND_PACKAGE_MODE NEVER)

FetchContent_Declare(external_cxxopts
   GIT_REPOSITORY https://github.com/jarro2783/cxxopts.git
   GIT_TAG v3.1.1
   GIT_SUBMODULES_RECURSE FALSE
   GIT_SHALLOW TRUE
   EXCLUDE_FROM_ALL
   SOURCE_SUBDIR dummy
)
FetchContent_MakeAvailable(external_cxxopts)

FetchContent_Declare(external_highfive
   GIT_REPOSITORY git@github.com:BlueBrain/HighFive.git
   GIT_TAG v2.9.0
   GIT_SUBMODULES_RECURSE FALSE
   GIT_SHALLOW TRUE
   EXCLUDE_FROM_ALL
   SOURCE_SUBDIR dummy
)
FetchContent_MakeAvailable(external_highfive)

FetchContent_Declare(external_isotope
   GIT_REPOSITORY git@github.com:Gregstrq/Isotope-data.git
   GIT_TAG main
   GIT_SUBMODULES_RECURSE FALSE
   GIT_SHALLOW TRUE
   EXCLUDE_FROM_ALL
   SOURCE_SUBDIR dummy
)
FetchContent_MakeAvailable(external_isotope)

FetchContent_Declare(external_json
   GIT_REPOSITORY git@github.com:nlohmann/json.git
   GIT_TAG v3.11.3
   GIT_SUBMODULES_RECURSE FALSE
   GIT_SHALLOW TRUE
   EXCLUDE_FROM_ALL
   SOURCE_SUBDIR dummy
)
FetchContent_MakeAvailable(external_json)

FetchContent_Declare(external_periodic
   GIT_REPOSITORY git@github.com:Bowserinator/Periodic-Table-JSON.git
   GIT_TAG v.4.0.0
   GIT_SUBMODULES_RECURSE FALSE
   GIT_SHALLOW TRUE
   EXCLUDE_FROM_ALL
   SOURCE_SUBDIR dummy
)
FetchContent_MakeAvailable(external_periodic)

if (BUILD_GUI)
   FetchContent_Declare(external_qmatplotwidget
      GIT_REPOSITORY https://gitlab.com/qdaq/qmatplotwidget.git
      GIT_TAG master
      GIT_SUBMODULES_RECURSE FALSE
      GIT_SHALLOW TRUE
      EXCLUDE_FROM_ALL
      SOURCE_SUBDIR dummy
   )
   FetchContent_MakeAvailable(external_qmatplotwidget)
endif()
