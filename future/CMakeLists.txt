cmake_minimum_required(VERSION 3.10)

project(SeAN)

set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -fprofile-arcs -ftest-coverage -fPIC -Wall -Wextra -Wshadow")
set(CMAKE_CXX_FLAGS "-O3 -march=native")

include_directories(include)

add_executable(sean src/main.cc)

if(CMAKE_BUILD_TYPE MATCHES DEBUG)
    enable_testing()

    add_executable(test_excited_state test/test_excited_state)
    add_executable(test_nucleus test/test_nucleus)

    add_test(NAME Execution COMMAND sean)
    add_test(NAME ExcitedState COMMAND test_excited_state)
    add_test(NAME Nucleus COMMAND test_nucleus)
    message(STATUS "Debug build activated.")
    message(STATUS "To run tests, execute\n\tmake test")
    message(STATUS "To view code coverage after testing, execute, e.g.,\ngcovr -r ${CMAKE_SOURCE_DIR} --exclude ${CMAKE_SOURCE_DIR}/include/TestUtilities.h\n   if 'gcovr' is installed.")
endif()