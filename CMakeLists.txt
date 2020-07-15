cmake_minimum_required (VERSION 3.1)

project(supercell)
include(FeatureSummary)

option(BUILD_PROGRAM "Build supercell program" ON)
option(BUILD_DOC "Build supercell manual and tutorial" ON)

if (BUILD_PROGRAM)
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O2 -Wall -Wno-sign-compare -Wextra ")

  option(LINK_STATIC_LIBS "Build project with static libs" OFF)

  #GLIBC_2.4
  option(GLIBC_24_COMPATIBLE "Build project compatible with GLIBC 2.4. Gold linker required." OFF)
  if (GLIBC_24_COMPATIBLE AND UNIX)
    set(WRAP_MEMCPY_SRC ${CMAKE_SOURCE_DIR}/libs/common_utils/others/memcpy.c)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,--wrap=memcpy")
  endif ()

  if (${CMAKE_COMPILER_IS_GNUCXX})
    if (${LINK_STATIC_LIBS})
      if(NOT WIN32)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libgcc -static-libstdc++ -fuse-ld=gold -Wl,--no-gnu-unique")
      else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libgcc -static-libstdc++ ")
      endif()
      set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    endif ()
  endif ()

  if (CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE Release)
  endif (CMAKE_BUILD_TYPE STREQUAL "")

  include_directories(${CMAKE_CURRENT_SOURCE_DIR})
  add_subdirectory(libs)
  add_subdirectory(src)

endif (BUILD_PROGRAM)

if(BUILD_DOC)
  add_subdirectory(doc)
endif(BUILD_DOC)
