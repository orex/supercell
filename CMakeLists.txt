cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

project(supercell)
include(FeatureSummary)
include(CTest)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(BUILD_PROGRAM "Build supercell program" ON)
option(BUILD_DOC "Build supercell manual and tutorial" ON)

if (BUILD_PROGRAM)
  if (CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE Release)
  endif (CMAKE_BUILD_TYPE STREQUAL "")

  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wno-sign-compare -Wextra ")

  option(LINK_STATIC_LIBS "Build project with static libs" OFF)

  if (${CMAKE_COMPILER_IS_GNUCXX})
    if (${LINK_STATIC_LIBS})
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libgcc -static-libstdc++")
      set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    endif ()
  endif ()

  include_directories(${CMAKE_CURRENT_SOURCE_DIR})
  add_subdirectory(src)

endif (BUILD_PROGRAM)

if(BUILD_DOC)
  add_subdirectory(doc)
endif(BUILD_DOC)
