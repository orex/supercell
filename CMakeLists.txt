cmake_minimum_required (VERSION 2.8)

project(supercell)
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)

set(CMAKE_BUILD_TYPE Debug)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

find_package(Boost 1.49)

add_subdirectory(libs)
add_subdirectory(src)
add_subdirectory(doc)