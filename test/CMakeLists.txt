ADD_EXECUTABLE(test_framework test_arithmetic_encoder.cc test_base.h
        ../cpp/src/ProbabilityEval.cpp
        ../cpp/include/ProbabilityEval.h
        ../cpp/single_include/nlohmann/json.hpp
        ../cpp/src/FreqTable.cpp
        ../cpp/include/FreqTable.h
        ../cpp/src/bitIO.cpp
        ../cpp/include/bitIO.h
        ../cpp/include/FastaParser.h
        ../cpp/include/commons.h
        ../cpp/src/commons.cpp
        ../cpp/single_include/argparse.hpp
        ../cpp/single_include/CRC.h
        ../cpp/src/ACBase.cpp
        ../cpp/include/ACBase.h
        ../cpp/src/ACEncode.cpp
        ../cpp/include/ACEncode.h
        ../cpp/src/ACDecode.cpp
        ../cpp/include/ACDecode.h
        ../cpp/src/ECDecoding.cpp
        ../cpp/include/ECDecoding.h
        ../cpp/include/debug_log.h
        ../cpp/single_include/robin_hood.h)
target_link_libraries(test_framework absl::btree absl::flat_hash_set absl::flat_hash_map)
ADD_TEST(test_framework test_arithmetic_encoder)
