#set common includes
set(sc_includes ${sc_includes} ${CMAKE_CURRENT_SOURCE_DIR}/core)

find_package(OpenBabel2 REQUIRED)

set(sc_includes ${sc_includes} ${OpenBabel2_INCLUDE_DIRS})
set(sc_includes ${sc_includes} ${CMAKE_SOURCE_DIR}/libs/common_utils/)

find_package(Eigen3 REQUIRED)
set(sc_includes ${sc_includes} ${EIGEN3_INCLUDE_DIR})

#set common sources
file(GLOB core_sources ${CMAKE_CURRENT_SOURCE_DIR}/core/*.cpp)
set(sc_sources ${WRAP_MEMCPY_SRC} ${sc_sources} ${core_sources})

set(sc_sources ${sc_sources} ${CMAKE_SOURCE_DIR}/libs/common_utils/others/string_utils.cpp
                             ${CMAKE_SOURCE_DIR}/libs/common_utils/others/rnd_utils.cpp
                             ${CMAKE_SOURCE_DIR}/libs/common_utils/cryst_tools/cryst_tools.cpp
                             ${CMAKE_SOURCE_DIR}/libs/common_utils/cryst_tools/comb_points.cpp)

set(sc_libraries ${sc_libraries} openbabel)

add_subdirectory(sc_cli)
#add_subdirectory(sc_qt)