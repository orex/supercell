set(Boost_USE_STATIC_LIBS ${LINK_STATIC_LIBS})
find_package(Boost 1.63 REQUIRED COMPONENTS program_options filesystem system)
include_directories(${Boost_INCLUDE_DIRS})

find_package(PkgConfig)

find_package(Threads REQUIRED)

find_package(Eigen3 REQUIRED)
add_definitions     ( ${EIGEN3_DEFINITIONS} )
include_directories ( ${EIGEN3_INCLUDE_DIRS} )

find_package(TBB COMPONENTS tbb)
if( NOT ${TBB_FOUND} )
    set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${PROJECT_SOURCE_DIR}/3rd_party/cmake/")
    find_package(TBB COMPONENTS tbb)
    if( NOT ${TBB_FOUND} )
        find_package(PkgConfig REQUIRED)
        pkg_search_module(PKGTBB REQUIRED IMPORTED_TARGET intel-tbb tbb)
        add_library(TBB::tbb ALIAS PkgConfig::PKGTBB)
    else()
        find_package(TBB REQUIRED COMPONENTS tbb)
    endif()
endif()

if( TBB_FOUND AND LINK_STATIC_LIBS )
    get_property(pmt TARGET TBB::tbb PROPERTY IMPORTED_LOCATION_RELEASE)
    string(REGEX REPLACE "libtbb\.(so|dylib)" "libtbb_static.a" pmt ${pmt})
    set_property(TARGET TBB::tbb PROPERTY IMPORTED_LOCATION_RELEASE ${pmt})
endif()

if (NOT TBB_USE_GLIBCXX_VERSION AND UNIX AND NOT APPLE)
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        # using Clang
        string(REPLACE "." "0" TBB_USE_GLIBCXX_VERSION ${CMAKE_CXX_COMPILER_VERSION})
    endif()
endif()

if (TBB_USE_GLIBCXX_VERSION)
    add_definitions(-DTBB_USE_GLIBCXX_VERSION=${TBB_USE_GLIBCXX_VERSION})
endif()

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
else()
    message(WARNING "libarchive is not found. Packed output will be disabled.")
endif()

include_directories(${sc_includes})
add_executable(supercell ${sc_sources} main.cpp parse_d2o_input.cpp)
target_link_libraries(supercell Boost::program_options Boost::filesystem Boost::system
        ${LibArchive_LIBRARIES} TBB::tbb Threads::Threads )

install(TARGETS supercell RUNTIME DESTINATION bin)
