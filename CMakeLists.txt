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


set(CMAKE_CXX_FLAGS "-g -Wall -std=c++11 -O3")


# Compiler-specific C++11 activation.
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    execute_process(
        COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
    if (NOT (GCC_VERSION VERSION_GREATER 4.7 OR GCC_VERSION VERSION_EQUAL 4.7))
        message(FATAL_ERROR "${PROJECT_NAME} requires g++ 4.7 or greater.")
    endif ()
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
else ()
    message(FATAL_ERROR "Your C++ compiler does not support C++11.")
endif ()

add_library(textileHelper textileHelper.cpp)
add_library(output output.cpp)
#add_library(SSEGenerator SSEGenerator.cpp)

# add the executable
add_executable(FromSSE FromSSE.cpp)
add_executable(GMSquare gmsquare.cpp)
add_executable(GMSSEs GMSSEs.cpp)
add_executable(OneOneFour OneOneFour.cpp)
add_executable(OneOneTwo OneOneTwo.cpp)
add_executable(ThreeByThreeA ThreeByThreeA.cpp)
add_executable(ThreeByThreeASquared ThreeByThreeASquared.cpp)
add_executable(ThreeByThreeB ThreeByThreeB.cpp)
add_executable(ThreeByThreeBSquared ThreeByThreeBSquared.cpp)
add_executable(ThreeByThreeC ThreeByThreeC.cpp)
add_executable(Odense Odense.cpp)
add_executable(MoreOdense MoreOdense.cpp)
add_executable(ODFATesting ODFATesting.cpp)


target_link_libraries(FromSSE output textileHelper)
target_link_libraries(GMSquare output textileHelper)
target_link_libraries(GMSSEs output textileHelper)
target_link_libraries(OneOneFour output textileHelper)
target_link_libraries(OneOneTwo output textileHelper)
target_link_libraries(ThreeByThreeA output textileHelper)
target_link_libraries(ThreeByThreeB output textileHelper)
target_link_libraries(ThreeByThreeASquared output textileHelper)
target_link_libraries(ThreeByThreeBSquared output textileHelper)
target_link_libraries(ThreeByThreeC output textileHelper)
target_link_libraries(Odense output textileHelper)
target_link_libraries(MoreOdense output textileHelper)
target_link_libraries(ODFATesting output textileHelper)