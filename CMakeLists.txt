cmake_minimum_required(VERSION 3.0.0)
project(tqsim-cpp VERSION 0.1.0 LANGUAGES CXX)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O0 -fno-omit-frame-pointer")
#set(CMAKE_OSX_ARCHITECTURES "arm64")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")

include(CTest)
enable_testing()

add_executable(tqsim-cpp main.cpp
        src/circuit.cpp)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

add_subdirectory(src)

#enable_testing()
#include(GoogleTest)
add_subdirectory(test)