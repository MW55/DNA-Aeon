ADD_EXECUTABLE(test_framework test_arithmetic_encoder.cc test_base.h
        ../src/ProbabilityEval.cpp
        ../include/ProbabilityEval.h
        ../single_include/nlohmann/json.hpp
        ../src/FreqTable.cpp
        ../include/FreqTable.h
        ../src/bitIO.cpp
        ../include/bitIO.h
        ../include/FastaParser.h
        ../include/commons.h
        ../src/commons.cpp
        ../single_include/argparse.hpp
        ../single_include/CRC.h
        ../src/ACBase.cpp
        ../include/ACBase.h
        ../src/ACEncode.cpp
        ../include/ACEncode.h
        ../src/ACDecode.cpp
        ../include/ACDecode.h
        ../src/ECDecoding.cpp
        ../include/ECDecoding.h
        ../include/debug_log.h
        ../single_include/robin_hood.h)
target_link_libraries(test_framework absl::btree absl::flat_hash_set absl::flat_hash_map)
ADD_TEST(test_framework test_arithmetic_encoder)
