# add_executable(gendedx
#     src/gendedx.cpp
#     src/periodic_table.h
#     ${PTABLE_DATA_CPP})

add_executable(genptable EXCLUDE_FROM_ALL
    source/src/genptable.cpp
)

FetchContent_MakeAvailable(external_isotope)
FetchContent_MakeAvailable(external_periodic)

set(ISOTOPE_CSV ${external_isotope_SOURCE_DIR}/isotopes_data.csv)
set(PTABLE_CSV ${external_periodic_SOURCE_DIR}/PeriodicTableCSV.csv)
set(PTABLE_DATA_CPP periodic_table.cpp)

add_custom_command(OUTPUT ${PTABLE_DATA_CPP}
    COMMAND genptable ${ISOTOPE_CSV} ${PTABLE_CSV} ${PTABLE_DATA_CPP}
    DEPENDS genptable
)

add_executable(gencorteo EXCLUDE_FROM_ALL
    source/src/gencorteo.cpp
)

FetchContent_MakeAvailable(external_cxxopts)

target_include_directories(gencorteo PRIVATE
    ${external_cxxopts_SOURCE_DIR}/include)

function(add_dedx_library screeningName)
    add_custom_command(OUTPUT xs_${screeningName}_data.cpp
        COMMAND gencorteo -s ${screeningName} DEPENDS gencorteo
    )
    add_library(xs_${screeningName} SHARED
        ${CMAKE_CURRENT_BINARY_DIR}/xs_${screeningName}_data.cpp
    )
endfunction()
