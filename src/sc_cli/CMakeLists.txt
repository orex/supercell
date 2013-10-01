set(Boost_USE_STATIC_LIBS ${LINK_STATIC_LIBS})

find_package(Boost 1.46 REQUIRED COMPONENTS regex program_options filesystem system)

set(sc_cli_source ${sc_sources} main.cpp parse_d2o_input.cpp)
set(sc_cli_libraries ${sc_libraries} ${Boost_LIBRARIES})

include_directories(${sc_includes})

add_executable(supercell ${sc_cli_source})

target_link_libraries(supercell ${sc_cli_libraries})

install(TARGETS supercell RUNTIME DESTINATION bin)