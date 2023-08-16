cmake_minimum_required(VERSION 3.19...3.27)
project(tqsim-cpp VERSION 0.1.0 LANGUAGES CXX)

# enable organization of targets into folders
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
	message(STATUS "Setting build type to 'Release' as none was specified.")
	set(CMAKE_BUILD_TYPE
	    Release
	    CACHE STRING "Choose the type of build." FORCE)
	# Set the possible values of build type for cmake-gui, ccmake
	set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel"
	             "RelWithDebInfo")
endif()

# Require C++ standard and disable standard extensions
set_property(GLOBAL PROPERTY CMAKE_CXX_STANDARD_REQUIRED ON)
set_property(GLOBAL PROPERTY CXX_EXTENSIONS OFF)

option(ENABLE_IPO "Enable Interprocedural Optimization, aka Link Time Optimization (LTO)" ON)
if(ENABLE_IPO)
	include(CheckIPOSupported)
	check_ipo_supported(RESULT ipo_supported OUTPUT ipo_output)
	# enable inter-procedural optimization if it is supported
	if(ipo_supported)
		set(CMAKE_INTERPROCEDURAL_OPTIMIZATION
		    TRUE
		    CACHE BOOL "Enable Interprocedural Optimization" FORCE)
	else()
		message(DEBUG "IPO is not supported: ${ipo_output}")
	endif()
endif()

# add main library code
add_subdirectory(src)

# add main executable
add_executable(tqsim main.cpp)
target_link_libraries(tqsim PUBLIC tqsim-cpp)

# add test code
option(BUILD_TQSIM_TESTS "Also build tests for the TQSIM project" ON)
if(BUILD_TQSIM_TESTS)
	enable_testing()
	include(GoogleTest)
	add_subdirectory(test)
endif()
