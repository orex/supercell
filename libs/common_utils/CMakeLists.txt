cmake_minimum_required (VERSION 2.8)

project(common_utils)

find_package(Boost 1.46 COMPONENTS unit_test_framework)

if(Boost_UNIT_TEST_FRAMEWORK_FOUND)
  enable_testing()
else()
  message(STATUS "No Boost::Test found. All tests have been disabled.")
endif()  

set(CMAKE_BUILD_TYPE Debug)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

if(Boost_UNIT_TEST_FRAMEWORK_FOUND)
  add_subdirectory(tests)
endif()    
