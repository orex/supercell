include_directories(${CMAKE_SOURCE_DIR}/others/)

add_executable(test_str test_str.cpp 
                        ${CMAKE_SOURCE_DIR}/others/string_utils.cpp)

target_link_libraries(test_str ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

add_test(NAME Test_string_utils_module COMMAND test_str)
