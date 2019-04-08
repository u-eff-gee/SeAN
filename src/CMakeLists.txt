cmake_minimum_required (VERSION 3.9 FATAL_ERROR)

include_directories("../include/")
add_library(sean_lib Absorption.cpp CrossSection.cpp Experiment.cpp InputReader.cpp Integrator.cpp PhononDensity.cpp Plotter.cpp Settings.cpp Target.cpp Writer.cpp)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS}$)
find_package(ROOT REQUIRED COMPONENTS MathMore)
include(${ROOT_USE_FILE})
target_link_libraries(sean_lib ${ROOT_LIBRARIES})

#set(FFTW3_DIR "/usr/local/lib/cmake/fftw3")
#message(STATUS "Looking for FFTW3Config.cmake in: ${FFTW3_DIR}. This path may have to be adapted manually in src/CMakeLists.txt")
find_package(FFTW3 REQUIRED)
target_link_libraries(sean_lib "${FFTW3_LIBRARY_DIRS}/libfftw3.so")

find_package(OpenMP REQUIRED)
target_link_libraries(sean_lib "${OpenMP_CXX}")
