cmake_minimum_required (VERSION 3.9 FATAL_ERROR)
project (sean)

add_subdirectory("src/")
include_directories("include/")

add_executable(sean src/SeAN.cpp)
target_link_libraries(sean sean_lib)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED)
include(${ROOT_USE_FILE})
message(STATUS "The ROOT version is: ${ROOT_VERSION}")
target_link_libraries(sean ${ROOT_LIBRARIES})

set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra -Wconversion -Wsign-conversion")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ftree_vectorize -march=native")
set(CMAKE_CXX_FLAGS_RELWITHDEBUGINFO "-g ${CMAKE_CXX_FLAGS_RELEASE}")
