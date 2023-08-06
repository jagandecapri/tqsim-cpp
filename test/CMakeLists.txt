include(FetchContent)
FetchContent_Declare(
        googletest
        # Specify the commit you depend on and update it regularly.
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.14.0 # release-1.10.0
        # URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../src)
# 'Google_Tests_run' is the target name
add_executable(Google_Tests_run test_circuit.cpp
        test_operator_generator.cpp
        test_basis_generator.cpp)
target_link_libraries(Google_Tests_run gtest gtest_main gmock gmock_main)
target_link_libraries(Google_Tests_run src)

include(GoogleTest)