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

# create an executable in which the tests will be stored
add_executable(tqsim-test
               test_circuit.cpp
               test_operator_generator.cpp
               test_basis_generator.cpp)

# link the Google test infrastructure and a default main function to the test executable.
target_link_libraries(tqsim-test PRIVATE tqsim-cpp gmock gtest_main)
# discover tests
gtest_discover_tests(
		tqsim-test
		WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
		PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}" DISCOVERY_TIMEOUT 60)
set_target_properties(tqsim-test PROPERTIES FOLDER tests)
