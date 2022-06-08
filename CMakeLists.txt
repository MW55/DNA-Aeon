cmake_minimum_required(VERSION 3.16.3)
project(arithmetic_modulator_error_correction)
set(ABSL_PROPAGATE_CXX_STD ON)
set(BUILD_TESTING OFF)
set(BUILD_TESTING_STATIC OFF)
set(BUILD_TESTING_SHARED OFF)

#Disabling Zippy extras
set(BUILD_SAMPLES OFF)
set(BUILD_TESTS OFF)
set(CREATE_DOCS OFF)

add_subdirectory(abseil-cpp)
add_subdirectory(Zippy)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DENABLE_DEBUG_MACRO -ggdb3")
set(CMAKE_CXX_STANDARD 20)
#set(CMAKE_CXX_STANDARD 11)
#-O4 -ffast-math
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -ansi -O4 -ffast-math -march=native")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -ansi -O4 -ffast-math -march=native")
#  -O4 -ffast-math -march=native -Q --help=target
# -pedantic
include_directories(cpp)

add_executable(arithmetic_modulator_error_correction
        cpp/src/ProbabilityEval.cpp
        cpp/include/ProbabilityEval.h
        cpp/single_include/nlohmann/json.hpp
        cpp/src/FreqTable.cpp
        cpp/include/FreqTable.h
        cpp/src/bitIO.cpp
        cpp/include/bitIO.h
        cpp/include/FastaParser.h
        cpp/src/commons.cpp
        cpp/include/commons.h
        cpp/single_include/argparse.hpp cpp/single_include/CRC.h cpp/single_include/robin_hood.h cpp/src/ACBase.cpp cpp/include/ACBase.h cpp/src/ACEncode.cpp cpp/include/ACEncode.h cpp/src/ACDecode.cpp cpp/include/ACDecode.h cpp/src/ECDecoding.cpp cpp/include/ECDecoding.h cpp/src/main.cpp cpp/include/debug_log.h)

target_link_libraries(arithmetic_modulator_error_correction absl::btree absl::flat_hash_set absl::flat_hash_map)

enable_testing ()
add_subdirectory (test)
