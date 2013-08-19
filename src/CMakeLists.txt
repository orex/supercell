#set common includes
set(sc_includes ${sc_includes} ${CMAKE_CURRENT_SOURCE_DIR}/core)

set(sc_includes ${sc_includes} /usr/local/include/openbabel-2.0/)
set(sc_includes ${sc_includes} ${CMAKE_SOURCE_DIR}/libs/utils/)

#set common sources
file(GLOB core_sources ${CMAKE_CURRENT_SOURCE_DIR}/core/*.cpp)
set(sc_sources ${sc_sources} ${core_sources})

set(sc_sources ${sc_sources} ${CMAKE_SOURCE_DIR}/libs/utils/string_utils.cpp)
set(sc_sources ${sc_sources} ${CMAKE_SOURCE_DIR}/libs/utils/rnd_utils.cpp)

#set_common_libraries
set(sc_libraries ${sc_libraries} openbabel)

add_subdirectory(sc_cli)
#add_subdirectory(sc_qt)