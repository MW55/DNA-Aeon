//
// Created by michael on 14.12.21.
//
//#include "iostream"
#include "string"
//#include "../../Zippy/library/Zippy.hpp"
#include "../../libraries/Zippy/library/Zippy.hpp"
#include "absl/container/flat_hash_set.h"
#include "absl/container/flat_hash_map.h"
#include "single_include/nlohmann/json.hpp"
#include "single_include/robin_hood.h"
#include "debug_log.h"

using namespace std;

uint64_t getMostFrequentElement(vector<uint64_t> &sizes) {
    if (sizes.empty())
        return -1;
    unordered_map<uint64_t, int> freq_count;
    for (const auto &item : sizes)
        freq_count[item]++;
    auto most_freq_int = std::max_element(freq_count.begin(), freq_count.end(),
                             [] (const auto &x, const auto &y) {return x.second < y.second;});
    return most_freq_int->first;
}


void write_to_zip(string outputFile, bool add_index, bool most_common_only, bool decodable_only, list<tuple<string, vector<unsigned char>>> &results) {
/*
 Writes the results to a zip file.
 @param outputFile: The name of the zip file.
 @param add_index: If true, the index is added to the file name.
 @param results: The results to write. Each result is a tuple of the name of the file and the data.
 */
    vector<uint64_t> sizes;
    std::vector<std::string> ents;
    Zippy::ZipArchive outarch;
    if (!outputFile.ends_with(".zip")) {
        outputFile += ".zip";
    }
    outarch.Create(outputFile);
    uint32_t i = 0;
    for (const tuple<string, vector<unsigned char>> &tuple1: results) {
        string filename = get<0>(tuple1);
        // optionally do not add packets that were not decodeable
        if (filename != "1000.000000" || decodable_only == false){
            if (add_index)
                filename = to_string(i).append("_").append(filename);
            else
                filename = filename.append("_enc");
            vector<unsigned char> res = get<1>(tuple1);
            // Add filename to a map of sizes, so that deviating packets can be removed ToDo make optional!
            outarch.AddEntry(filename, res);
            sizes.push_back(res.size());
            ents.push_back(filename);
            i++;
        }
    }
    //Have to save it before being able to get the sizes.
    outarch.Save();
    if(most_common_only) {
    // optionally removing all entries that differ in size from the most frequent element
        uint_fast64_t mostCommonSize = getMostFrequentElement(sizes);
        for (auto &ent : ents){
            Zippy::ZipEntry b = outarch.GetEntry(ent);
            if (b.UncompressedSize() != mostCommonSize){
                outarch.DeleteEntry(ent);
            }
        }
    }
    outarch.Save();
    outarch.Close();
    DEBUG("Wrote to zip file: " << outputFile);
}

robin_hood::unordered_map<string, vector<string>> readConcScheme(const string &concScheme) {
    ifstream ifs(concScheme);
    nlohmann::json jsonFile = nlohmann::json::parse(ifs);
    robin_hood::unordered_map<string, vector<string>> motif;
    for (const auto &[key, val]: jsonFile["motif"].items()) {
        auto vec = jsonFile["motif"][key].get<vector<string>>();
        motif[key] = vec;
    }
    return motif;
}

/**
 * @brief Get the Codeword Len object
 * returns the size of the first element of the set
 * @param codewords 
 * @return int 
 */

int getCodewordLen(robin_hood::unordered_set<string> &codewords) {
    return codewords.begin()->length();
}
