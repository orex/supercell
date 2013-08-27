cmake_minimum_required (VERSION 2.8)

project(supercell)
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wno-sign-compare -Wextra " )

if(CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE Release)
endif(CMAKE_BUILD_TYPE STREQUAL "")

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(libs)
add_subdirectory(src)
add_subdirectory(doc)
