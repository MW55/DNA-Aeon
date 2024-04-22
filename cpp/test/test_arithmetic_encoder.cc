#include <mutex>
#include <utility>
#include <filesystem>
#include <absl/synchronization/internal/thread_pool.h>
#include "test_base.h"
#include "include/commons.h"
#include "include/FastaParser.h"
#include "include/FreqTable.h"
#include "include/ACEncode.h"
#include "include/ECDecoding.h"
#include "../Zippy/library/Zippy/ZipArchive.hpp"

#define TEST_FILENAME "tmp_test_output.file"
#define MULTITHREADED
using namespace std;

list<tuple<string, vector<unsigned char>>> results = list<tuple<string, vector<unsigned char>>>();
std::mutex *res_lock;

void test_data_load() {
    robin_hood::unordered_set<string> codewords = parseFasta("./codewords/t2.fasta");
    robin_hood::unordered_map<string, vector<string>> motif = readConcScheme("./codewords/motifs.json");
    int codewordLen = getCodewordLen(codewords);
    size_t motifLen = motif.size();
    ALEPH_ASSERT_EQUAL(codewordLen, 10);
    ALEPH_ASSERT_EQUAL(motifLen, 31);
}

void test_simple_en_decode(tuple<uint8_t, string, string, string, string> tup, nlohmann::json config) {
    auto[sync, filename, data, pred1, pred2] = std::move(tup);
    results = list<tuple<string, vector<unsigned char>>>();
    res_lock = new std::mutex();
    robin_hood::unordered_set<string> codewords = parseFasta("./codewords/t2.fasta");
    robin_hood::unordered_map<string, vector<string>> motif = readConcScheme("./codewords/motifs.json");
    int codewordLen = getCodewordLen(codewords);
    auto propMap = ProbMap(codewordLen, false, codewords, motif);
    robin_hood::unordered_map<string, int32_t> freqDict = propMap.freqDict();
    robin_hood::unordered_map<string, char2double> pMap = propMap.createTransitionDict(freqDict);
    results = list<tuple<string, vector<unsigned char>>>();
    FreqTable freqs = FreqTable(motif, "", 0, false, codewordLen, pMap);
    freqs.calcFreqs();
    do_encode(data, sync, freqs, 0, filename, results, res_lock);
    ALEPH_ASSERT_EQUAL(results.size(), 1);
    ALEPH_ASSERT_STRING_EQUAL(static_cast<string>(get<0>(results.front())), filename);
    vector<unsigned char> msg = get<1>(results.front());
    string s(msg.begin(), msg.end());
    //array<int, 2> rate{3, 4};
    robin_hood::unordered_map<string, char2double> tMap = ProbMap(codewordLen, true, codewords,
                                                            motif).createTransitionDict(freqDict);
    //uint8_t threshold = 0;
    results.clear();
    config["decode"]["length"] = s.size();
    config["general"]["sync"] = sync;
    std::array<int, 2> e_rate{config["decode"]["metric"]["fano"]["rate"]["low"], config["decode"]["metric"]["fano"]["rate"]["high"]};
    do_decode(s, freqs, tMap, motif, config, codewordLen, results, res_lock, e_rate);
    ALEPH_ASSERT_EQUAL(results.size(), 1);
    msg = get<1>(results.front());
    string t(msg.begin(), msg.end());
    ALEPH_ASSERT_STRING_EQUAL(static_cast<string>(get<0>(results.front())), pred1);
    ALEPH_ASSERT_STRING_EQUAL(t, data);
    string s2 = s;
    s2.at(s2.size() - 1) = 'T';
    //s2.at(s2.size() - 3) = 'C';
    results.clear();
    config["decode"]["length"] = s2.size();
    do_decode(s2, freqs, tMap, motif, config, codewordLen, results, res_lock, e_rate);
    ALEPH_ASSERT_EQUAL(results.size(), 1);
    msg = get<1>(results.front());
    string t2(msg.begin(), msg.end());
    ALEPH_ASSERT_STRING_EQUAL(static_cast<string>(get<0>(results.front())), pred2);
    ALEPH_ASSERT_STRING_EQUAL(t2, data);
}

void test_load_zip_store_fasta_circle(const robin_hood::unordered_set<string> &input_data, uint8_t sync, nlohmann::json config) {
    //array<int, 2> rate{3, 4};
    robin_hood::unordered_set<string> codewords = parseFasta("./codewords/t2.fasta");
    robin_hood::unordered_map<string, vector<string>> motif = readConcScheme("./codewords/motifs.json");
    int codewordLen = getCodewordLen(codewords);
    size_t motifLen = motif.size();
    ALEPH_ASSERT_EQUAL(codewordLen, 10);
    ALEPH_ASSERT_EQUAL(motifLen, 31);
    res_lock = new std::mutex();

    auto propMap = ProbMap(codewordLen, false, codewords, motif);
    robin_hood::unordered_map<string, int32_t> freqDict = propMap.freqDict();
    robin_hood::unordered_map<string, char2double> pMap = propMap.createTransitionDict(freqDict);
    results = list<tuple<string, vector<unsigned char>>>();
    FreqTable freqs = FreqTable(motif, "", 0, false, codewordLen, pMap);
    freqs.calcFreqs();
    absl::synchronization_internal::ThreadPool pool(
            static_cast<int>(min(static_cast<unsigned long>(std::thread::hardware_concurrency()),
                                 input_data.size())));
    for (auto &data: input_data) {
#ifdef MULTITHREADED
        pool.Schedule([data, sync, &capture1 = freqs, &capture0 = results] {
            return do_encode(data, sync, capture1, 0, "ent", capture0, res_lock);
        });
#else
        do_encode(data, sync, freqs, 0, "ent", std::ref(results), res_lock);
#endif
    }
    while (results.size() != input_data.size()) {
        this_thread::sleep_for(chrono::milliseconds(500));
    }
    ALEPH_ASSERT_EQUAL(results.size(), input_data.size());

    //test fasta output:
    remove(TEST_FILENAME);
    write_to_fasta(TEST_FILENAME, results);
    // re-read the fasta file and compare the results:
    robin_hood::unordered_set<string> content = parseFasta(TEST_FILENAME);
    for (auto &ent: results) {
        string str(get<1>(ent).begin(), get<1>(ent).end());
        ALEPH_ASSERT_THROW(content.contains(str));
    }
    //test zip output
    remove(TEST_FILENAME);
    string filename = TEST_FILENAME;
    filename += ".zip";
    write_to_zip(filename, true, true, true, results);

    ALEPH_ASSERT_THROW(std::filesystem::exists(filename));
    results.clear();
    // read the zip file and compare the results:
    Zippy::ZipArchive inarch(TEST_FILENAME);
    std::vector<std::string> ents = inarch.GetEntryNames();
    for (auto &ent: ents) {
        string inp = inarch.GetEntry(ent).GetDataAsString();
        //int64_t messageLength = static_cast<int64_t>(inp.size());
        config["decode"]["length"] = static_cast<int64_t>(inp.size());
        config["general"]["sync"] = sync;
#ifdef MULTITHREADED
        pool.Schedule(
                [inp, &capture = freqs, &capture4 = pMap, &capture5 = motif, &capture6 = config, &capture7 = codewordLen, &capture1 = results] {
                    return do_decode(inp, capture, capture4, capture5, capture6, capture7, capture1,
                                     res_lock, e_rate);
                });
#else
        do_decode(inp, freqs, pMap, motif, config, codewordLen, std::ref(results), res_lock);
#endif
        // make sure all entries in the zip exist in the fasta file
        ALEPH_ASSERT_THROW(content.contains(ent));
    }
    ALEPH_ASSERT_EQUAL(results.size(), ents.size());
    for (auto &elem: results) {
        string str(get<1>(elem).begin(), get<1>(elem).end());
        ALEPH_ASSERT_THROW(input_data.contains(str));
    }
    remove(TEST_FILENAME);
}

int main(int, char **) {
    robin_hood::unordered_set<string> input_data = {"HelloWorld"s, "test1235"s, "long test 123"s, "\x00""abcdefghij"s,
                                              "test""\x00""123"s, "1234567891234567",
                                              "12345678912345671234567891234567", "12345678",
                                              "DÃ\u00151v;\bu<)@cD%t\u001Ba\u0019(y\u007Fx9¬‹\u001D3›D”ÒÉ"s};
    vector<tuple<uint8_t, string, string, string, string>> t = {{5, "filename1", "data",      "-10.907847", "-4.659920"},
                                                                {3, "filename1", "test123",   "-20.828143", "-20.815485"},
                                                                {4, "filename1", "HalloWelt", "-26.054357", "-26.048063"}};
    
    auto config = R"(
    {
        "general": {
            "sync": 4,
            "as_fasta": false,
            "codebook": {
                "words": "./codewords/codebook_test.fasta",
                "motifs": "./codewords/codebook_test.json"
            },
            "threads": 4,
            "zip": {
                "most_common_only": true,
                "decodable_only": true
            }
        },
        "encode":{
            "input": "data\\test.txt",
            "output": "data\\encoded.txt",
            "min_length": 0,
            "same_length": false,
            "update_config": false
        },
        "decode":{
            "input": "data\\encoded.txt",
            "output": "data\\decoded.txt",
            "length": 0,
            "threshold": {
                "loop": 1,
                "finish": 0,
                "checkpoint": 3
            },
            "metric": {
                "fano": {
                    "rate": {"low": 3, "high": 4},
                    "error_probability": 0.05
                },
                "penalties": {
                    "crc": 0.1,
                    "no_hit": 8
                }
            },
            "queue": {
                "size": 200000,
                "runs": 5,
                "reduce": 0.5
            }
        }
    }
    )"_json;
    
    test_data_load();
    for (auto &i: t) {
        test_simple_en_decode(i, config);
    }
    test_load_zip_store_fasta_circle(input_data, 5, config);
    this_thread::sleep_for(chrono::milliseconds(1000));
}
