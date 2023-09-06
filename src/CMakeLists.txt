file(GLOB_RECURSE EXTERN_HEADER_FILES ${PROJECT_SOURCE_DIR}/extern/Eigen/*.hpp ${PROJECT_SOURCE_DIR}/extern/pcg/*.hpp)
file(GLOB_RECURSE HEADER_FILES ${PROJECT_SOURCE_DIR}/include/*.hpp)
file(GLOB_RECURSE SOURCE_FILES ${PROJECT_SOURCE_DIR}/src/*.cpp)

set(target_name tqsim-cpp)
add_library(${target_name} STATIC ${SOURCE_FILES} ${HEADER_FILES} ${EXTERN_HEADER_FILES})

# set include directories
target_include_directories(${target_name} PUBLIC ${PROJECT_SOURCE_DIR}/extern)
target_include_directories(${target_name} PUBLIC ${PROJECT_SOURCE_DIR}/include)

# add MQT::Core target
set(BUILD_MQT_CORE_TESTS
		OFF
		CACHE BOOL "Build MQT Core tests")
add_subdirectory("${PROJECT_SOURCE_DIR}/extern/mqt-core" "extern/mqt-core"
		EXCLUDE_FROM_ALL)

# link to the MQT::Core libraries
target_link_libraries(${PROJECT_NAME} PUBLIC MQT::CoreDD)

# set required C++ standard
target_compile_features(${target_name} PUBLIC cxx_std_17)

# always include debug symbols (avoids common problems with LTO)
target_compile_options(${target_name} PUBLIC -g)

# enable platform-dependent optimizations
target_compile_options(${target_name} PUBLIC -mtune=native)
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-march=native HAS_MARCH_NATIVE)
if(HAS_MARCH_NATIVE)
	target_compile_options(${target_name} PUBLIC -march=native)
endif()

# enable some more optimizations in release mode
target_compile_options(${target_name} PUBLIC $<$<CONFIG:RELEASE>:-fno-math-errno
                       -ffinite-math-only -fno-trapping-math>)

# enable some more options for better debugging
target_compile_options(
		${target_name} PUBLIC $<$<CONFIG:DEBUG>:-fno-omit-frame-pointer
		-fno-optimize-sibling-calls -fno-inline-functions>)
