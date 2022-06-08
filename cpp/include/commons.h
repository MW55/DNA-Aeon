//
// Created by michael on 14.12.21.
//

#ifndef ARITHMETIC_MODULATOR_ERROR_CORRECTION_COMMONS_H
#define ARITHMETIC_MODULATOR_ERROR_CORRECTION_COMMONS_H

#include "absl/container/flat_hash_set.h"
#include "absl/container/flat_hash_map.h"
#include "single_include/robin_hood.h"

using namespace std;

extern void write_to_zip(string outputFile, bool add_index, bool most_common_only, bool decodable_only, list <tuple<string, vector<unsigned char>>> &results);

extern robin_hood::unordered_map<string, vector<string>> readConcScheme(const string &concScheme);

extern int getCodewordLen(robin_hood::unordered_set<string> &codewords);

#endif //ARITHMETIC_MODULATOR_ERROR_CORRECTION_COMMONS_H
