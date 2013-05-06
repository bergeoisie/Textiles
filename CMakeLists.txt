cmake_minimum_required (VERSION 2.6)
project (textiles)
 
#find_package(armadillo)
find_package(Boost)

# configure a header file to pass some of the CMake settings
# to the source code
#configure_file (
#  "${PROJECT_SOURCE_DIR}/caentConfig.h.in"
#  "${PROJECT_BINARY_DIR}/caentConfig.h"
#  )
 
# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
include_directories("${PROJECT_BINARY_DIR}")

set(CMAKE_CXX_FLAGS "-g -Wall -std=c++11 -stdlib=libc++ -O3")

add_library(textileHelper textileHelper.cpp)
add_library(output output.cpp)

# add the executable
add_executable(FromSSE FromSSE.cpp)
add_executable(GMSquare gmsquare.cpp)

target_link_libraries(FromSSE output textileHelper)
target_link_libraries(GMSquare output textileHelper)