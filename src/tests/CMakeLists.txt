include_directories(${CMAKE_SOURCE_DIR}/others/)

find_package(Boost 1.46 REQUIRED COMPONENTS unit_test_framework)

find_package(Eigen3 REQUIRED)
add_definitions     ( ${EIGEN3_DEFINITIONS} )
include_directories ( ${EIGEN3_INCLUDE_DIRS} )

if(Boost_FOUND)
    add_executable(test_random_sampler ${CMAKE_SOURCE_DIR}/src/core/common_types.cpp test_random_sampler.cpp)
    target_link_libraries(test_random_sampler Boost::unit_test_framework)
    add_test(NAME TestRandomSampler COMMAND test_random_sampler)
endif()

#add_subdirectory(regression)