#set common includes
set(sc_includes ${sc_includes} ${CMAKE_CURRENT_SOURCE_DIR}/core)

find_package(PkgConfig)

find_package(Threads REQUIRED)
set(sc_includes ${sc_includes} ${CMAKE_SOURCE_DIR}/libs/common_utils/)

find_package(Eigen3 QUIET)
if(${EIGEN3_FOUND})
  message(STATUS "Eigen3 library FOUND by cmake")
endif()

# The following workaround only need for Ubuntu 14.04 with non-default version of cmake (travis-ci case).
if(NOT DEFINED EIGEN3_FOUND)
  find_path(EIGEN3_INCLUDE_DIR NAMES Eigen/Core
      HINTS
      ${CMAKE_INSTALL_PREFIX}/include
      PATH_SUFFIXES eigen3 
    )  
    if(EXISTS "${EIGEN3_INCLUDE_DIR}/Eigen/Core")
    set(EIGEN3_FOUND TRUE)
    message(STATUS "Eigen FOUND directly")
  else()
  endif()
endif()

if(NOT ${EIGEN3_FOUND})
  message(FATAL_ERROR "Eigen3 library not FOUND!")
endif()

find_package(TBB COMPONENTS tbb)
if( NOT ${TBB_FOUND} )
  if( ${PkgConfig_FOUND} )
    pkg_check_modules(TBB REQUIRED IMPORTED_TARGET intel-tbb)
    add_library(TBB::tbb ALIAS PkgConfig::TBB)
  else()
    find_package(TBB REQUIRED COMPONENTS tbb)
  endif()
elseif( LINK_STATIC_LIBS )
  get_property(pmt TARGET TBB::tbb PROPERTY IMPORTED_LOCATION_RELEASE)
  string(REPLACE "libtbb.so" "libtbb_static.a" pmt ${pmt})
  set_property(TARGET TBB::tbb PROPERTY IMPORTED_LOCATION_RELEASE ${pmt})
endif()


set(sc_includes ${sc_includes} ${EIGEN3_INCLUDE_DIR}
        ${CMAKE_SOURCE_DIR}/3rd_party/gemmi/include
        ${CMAKE_SOURCE_DIR}/3rd_party/xxHash)


find_package(LibArchive)
if( ${LibArchive_FOUND} )
  set(sc_includes ${sc_includes} ${LibArchive_INCLUDE_DIRS})
  add_definitions(-DLIBARCHIVE_ENABLED)
  set(CMAKE_REQUIRED_INCLUDES ${LibArchive_INCLUDE_DIRS})
  if(LINK_STATIC_LIBS)
    find_package(BZip2 REQUIRED)
    find_package(LibLZMA REQUIRED)
    find_package(ZLIB REQUIRED)
    set(LibArchive_LIBRARIES ${LibArchive_LIBRARIES} ${LIBLZMA_LIBRARIES} ${BZIP2_LIBRARIES} ${ZLIB_LIBRARIES})
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
set(sc_sources ${sc_sources} ${core_sources})

set(sc_sources ${sc_sources} ${CMAKE_SOURCE_DIR}/libs/common_utils/others/string_utils.cpp
                             ${CMAKE_SOURCE_DIR}/libs/common_utils/cryst_tools/cryst_tools.cpp
                             ${CMAKE_SOURCE_DIR}/libs/common_utils/cryst_tools/comb_points.cpp)

add_subdirectory(sc_cli)
