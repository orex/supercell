#set common includes
set(sc_includes ${sc_includes} ${CMAKE_CURRENT_SOURCE_DIR}/core)

find_package(OpenBabel2 REQUIRED 2.4.1)

add_definitions(-DOLD_OB_PERIODIC_TABLE)
foreach(OB_ID ${OpenBabel2_INCLUDE_DIRS})
  if(EXISTS ${OB_ID}/openbabel/elements.h)
    remove_definitions(-DOLD_OB_PERIODIC_TABLE)
  endif()
endforeach()

set(sc_includes ${sc_includes} ${OpenBabel2_INCLUDE_DIRS})
set(sc_includes ${sc_includes} ${CMAKE_SOURCE_DIR}/libs/common_utils/)

find_package(Eigen3 QUIET)
if(${EIGEN3_FOUND})
  message(STATUS "Eigen3 library FOUND by cmake")
endif()

# The following workaround only need for Ubuntu 14.04 with non-default version of cmake (travis-ci case).
if(NOT DEFINED EIGEN3_FOUND)
  if(EXISTS "/usr/include/eigen3/Eigen/Core")
    set(EIGEN3_INCLUDE_DIR "/usr/include/eigen3/")
    set(EIGEN3_FOUND TRUE)
    message(STATUS "Eigen FOUND directly")
  else()
    message(FATAL_ERROR "Eigen3 library not FOUND!")
  endif()
endif()
set(sc_includes ${sc_includes} ${EIGEN3_INCLUDE_DIR})

find_package(LibArchive)
if( ${LibArchive_FOUND} )
  set(sc_includes ${sc_includes} ${LibArchive_INCLUDE_DIRS})
  add_definitions(-DLIBARCHIVE_ENABLED)
  set(CMAKE_REQUIRED_INCLUDES ${LibArchive_INCLUDE_DIRS})
  if(LINK_STATIC_LIBS)
    find_package(BZip2 REQUIRED)
    find_package(LibLZMA REQUIRED)
    set(CMAKE_REQUIRED_LIBRARIES ${LibArchive_LIBRARIES} ${LIBLZMA_LIBRARIES} ${BZIP2_LIBRARIES})
    set(LibArchive_LIBRARIES ${LibArchive_LIBRARIES} ${LIBLZMA_LIBRARIES} ${BZIP2_LIBRARIES})
  else()
    set(CMAKE_REQUIRED_LIBRARIES ${LibArchive_LIBRARIES})
  endif()
  include(CheckFunctionExists) 
  CHECK_FUNCTION_EXISTS(archive_write_set_format_filter_by_ext LibArchive_format_filter_by_ext)
  if( ${LibArchive_format_filter_by_ext} )
    add_definitions( -DLIBARCHIVE_PATCH_DISABLE )
  endif() 
else()
  message(WARNING "libarchive is not found. Packed output will be disabled.")
endif()

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
