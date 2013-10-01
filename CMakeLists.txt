cmake_minimum_required (VERSION 2.8)

project(supercell)
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wno-sign-compare -Wextra " )

option(LINK_STATIC_LIBS "Build project with static libs" OFF)

if(LINK_STATIC_LIBS)
  set(OB_LIB_NAMES libopenbabel.a openbabel.lib openbabel)
else()
  set(OB_LIB_NAMES libopenbabel.so openbabel.dll openbabel)
endif()

if( ${CMAKE_COMPILER_IS_GNUCXX} AND ${LINK_STATIC_LIBS} )
  set(CMAKE_CXX_LINK_EXECUTABLE "${CMAKE_CXX_LINK_EXECUTABLE} -static-libstdc++ -static-libgcc")
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
endif()


find_library(OB_LIB openbabel NAMES ${OB_LIB_NAMES})

if(CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE Release)
endif(CMAKE_BUILD_TYPE STREQUAL "")

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(libs)
add_subdirectory(src)
add_subdirectory(doc)
