project(tqsim-cpp)

set(HEADER_FILES
        circuit.h
        operator_generator.h
)

set(SOURCE_FILES
        circuit.cpp
        operator_generator.cpp
)

add_library(src STATIC ${SOURCE_FILES} ${HEADER_FILES})