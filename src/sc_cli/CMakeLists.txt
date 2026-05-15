set(Boost_USE_STATIC_LIBS ${LINK_STATIC_LIBS})
find_package(Boost 1.71 REQUIRED COMPONENTS program_options)
include_directories(${Boost_INCLUDE_DIRS})

find_package(PkgConfig)

find_package(Threads REQUIRED)

find_package(Eigen3 REQUIRED)

find_package(TBB CONFIG COMPONENTS tbb)
if(NOT TBB_FOUND)
    find_package(PkgConfig REQUIRED)
    pkg_search_module(PKGTBB REQUIRED IMPORTED_TARGET tbb)
    add_library(TBB::tbb ALIAS PkgConfig::PKGTBB)
endif()

if(LINK_STATIC_LIBS AND TBB_FOUND)
    get_property(pmt TARGET TBB::tbb PROPERTY IMPORTED_LOCATION_RELEASE)
    if(pmt)
        string(REGEX REPLACE "libtbb\\.(so|dylib)" "libtbb_static.a" pmt ${pmt})
        set_property(TARGET TBB::tbb PROPERTY IMPORTED_LOCATION_RELEASE ${pmt})
    endif()
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
        # libarchive 3.6+ links zstd / lz4 transitively when those headers are
        # present at build time. Pull them into the static link list if found.
        find_library(_zstd_lib NAMES zstd)
        if(_zstd_lib)
            list(APPEND LibArchive_LIBRARIES ${_zstd_lib})
        endif()
        find_library(_lz4_lib NAMES lz4)
        if(_lz4_lib)
            list(APPEND LibArchive_LIBRARIES ${_lz4_lib})
        endif()
    endif()
else()
    message(WARNING "libarchive is not found. Packed output will be disabled.")
endif()

include_directories(${sc_includes})
add_executable(supercell ${sc_sources} main.cpp parse_d2o_input.cpp)
target_link_libraries(supercell Boost::program_options Eigen3::Eigen
        ${LibArchive_LIBRARIES} TBB::tbb Threads::Threads )

# Optional: pin libm exp() to its GLIBC_2.2.5 symbol so the binary loads on
# pre-2.29 glibc hosts (RHEL 8 / Ubuntu 18.04 / Debian 10 / Leap 15). x86_64
# Linux only; documented in glibc_compat.c.
option(LEGACY_GLIBC_EXP
       "Pin libm exp() to its GLIBC_2.2.5 symbol via a wrap shim (x86_64 Linux only)."
       OFF)
if(LEGACY_GLIBC_EXP)
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|amd64")
        target_sources(supercell PRIVATE glibc_compat.c)
        # GCC's --wrap is a final-link-stage rewrite that doesn't compose with
        # -flto: LTO discards __wrap_exp as "unused" at the IR level before
        # --wrap gets a chance to redirect calls to it. Compile this TU
        # without LTO so the wrapper survives as a regular ELF symbol.
        set_source_files_properties(glibc_compat.c PROPERTIES
            COMPILE_OPTIONS "-fno-lto")
        target_link_options(supercell PRIVATE "LINKER:--wrap=exp")
    else()
        message(WARNING "LEGACY_GLIBC_EXP only applies on x86_64 Linux; ignoring.")
    endif()
endif()

install(TARGETS supercell RUNTIME DESTINATION bin)
