cmake_minimum_required (VERSION 2.8.6)

project(supercell)
include(FeatureSummary)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O3 -Wall -Wno-sign-compare -Wextra " )

option(LINK_STATIC_LIBS "Build project with static libs" OFF)

if(LINK_STATIC_LIBS)
  set(OB_LIB_NAMES libopenbabel.a openbabel.lib openbabel)
else()
  set(OB_LIB_NAMES libopenbabel.so openbabel.dll openbabel)
endif()

#GLIBC_2.4
option(GLIBC_24_COMPATIBLE "Build project compatible with GLIBC 2.4. Gold linker required." OFF)
if(GLIBC_24_COMPATIBLE AND UNIX)
  set(WRAP_MEMCPY_SRC ${CMAKE_SOURCE_DIR}/libs/common_utils/others/memcpy.c)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,--wrap=memcpy")  
endif()

if( ${CMAKE_COMPILER_IS_GNUCXX} )
  if( ${LINK_STATIC_LIBS} )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libgcc -static-libstdc++ -Wl,-no-gnu-unique")
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  endif()
endif()

find_library(OB_LIB openbabel NAMES ${OB_LIB_NAMES})

if(CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE Release)
endif(CMAKE_BUILD_TYPE STREQUAL "")

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(libs)
add_subdirectory(src)
add_subdirectory(doc)