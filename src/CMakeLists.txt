#set common includes
set(sc_includes ${sc_includes} ${CMAKE_CURRENT_SOURCE_DIR}/core)

set(sc_includes ${sc_includes} /usr/local/include/openbabel-2.0/)
set(sc_includes ${sc_includes} ${CMAKE_SOURCE_DIR}/libs/common_utils/)

#set common sources
file(GLOB core_sources ${CMAKE_CURRENT_SOURCE_DIR}/core/*.cpp)
set(sc_sources ${WRAP_MEMCPY_SRC} ${sc_sources} ${core_sources})

set(sc_sources ${sc_sources} ${CMAKE_SOURCE_DIR}/libs/common_utils/others/string_utils.cpp
                             ${CMAKE_SOURCE_DIR}/libs/common_utils/others/rnd_utils.cpp)

#set_common_libraries
add_library(openbabel UNKNOWN IMPORTED)
set_property(TARGET openbabel PROPERTY IMPORTED_LOCATION "${OB_LIB}")

set(sc_libraries ${sc_libraries} openbabel)

add_subdirectory(sc_cli)
#add_subdirectory(sc_qt)