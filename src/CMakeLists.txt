#set common includes
set(sc_includes
        ${CMAKE_CURRENT_SOURCE_DIR}/core
        ${CMAKE_CURRENT_SOURCE_DIR}/utils
        ${CMAKE_SOURCE_DIR}/3rd_party/gemmi/include
        ${CMAKE_SOURCE_DIR}/3rd_party/xxHash)

#set common sources
file(GLOB core_sources ${CMAKE_CURRENT_SOURCE_DIR}/core/*.cpp)
set(sc_sources ${core_sources}
        ${CMAKE_CURRENT_SOURCE_DIR}/utils/others/string_utils.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/utils/cryst_tools/cryst_tools.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/utils/cryst_tools/comb_points.cpp)

add_subdirectory(sc_cli)
add_subdirectory(tests)
