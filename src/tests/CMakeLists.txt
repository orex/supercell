set(Boost_USE_STATIC_LIBS ${LINK_STATIC_LIBS})
find_package(Boost 1.71 REQUIRED COMPONENTS unit_test_framework)
find_package(Eigen3 REQUIRED)

if(NOT Boost_FOUND)
    return()
endif()

set(_src_root ${CMAKE_CURRENT_SOURCE_DIR}/..)
set(_utils    ${_src_root}/utils)
set(_core     ${_src_root}/core)
set(_3rdparty ${CMAKE_SOURCE_DIR}/3rd_party)

# Shared include set for all tests.
set(_test_includes
    ${_src_root}
    ${_core}
    ${_utils}
    ${_3rdparty}/gemmi/include
    ${_3rdparty}/xxHash)

# ---- test_random_sampler ----
add_executable(test_random_sampler
        ${_core}/common_types.cpp
        test_random_sampler.cpp)
target_include_directories(test_random_sampler PRIVATE ${_test_includes})
target_link_libraries(test_random_sampler Boost::unit_test_framework Eigen3::Eigen)
add_test(NAME TestRandomSampler COMMAND test_random_sampler)

# ---- test_string_utils ----
add_executable(test_string_utils
        ${_utils}/others/string_utils.cpp
        test_string_utils.cpp)
target_include_directories(test_string_utils PRIVATE ${_test_includes})
target_link_libraries(test_string_utils Boost::unit_test_framework)
add_test(NAME TestStringUtils COMMAND test_string_utils)

# ---- test_comb_points ----
add_executable(test_comb_points
        ${_utils}/cryst_tools/comb_points.cpp
        test_comb_points.cpp)
target_include_directories(test_comb_points PRIVATE ${_test_includes})
target_link_libraries(test_comb_points Boost::unit_test_framework)
add_test(NAME TestCombPoints COMMAND test_comb_points)

# ---- test_combinatorics ----
add_executable(test_combinatorics test_combinatorics.cpp)
target_include_directories(test_combinatorics PRIVATE ${_test_includes})
target_compile_definitions(test_combinatorics PRIVATE
        COMBINATORICS_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data_combinatorics")
target_link_libraries(test_combinatorics Boost::unit_test_framework)
add_test(NAME TestCombinatorics COMMAND test_combinatorics)

# ---- test_cryst_tools (gemmi-based, was OpenBabel) ----
add_executable(test_cryst_tools
        ${_core}/common_types.cpp
        ${_core}/cif_io.cpp
        ${_utils}/others/string_utils.cpp
        ${_utils}/cryst_tools/cryst_tools.cpp
        ${_utils}/cryst_tools/comb_points.cpp
        test_cryst_tools.cpp)
target_include_directories(test_cryst_tools PRIVATE ${_test_includes})
target_compile_definitions(test_cryst_tools PRIVATE
        CT_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data_cryst_tools")
target_link_libraries(test_cryst_tools Boost::unit_test_framework Eigen3::Eigen)
add_test(NAME TestCrystTools COMMAND test_cryst_tools)
