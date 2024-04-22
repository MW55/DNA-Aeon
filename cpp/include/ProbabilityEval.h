//
// Created by wintermute on 9/8/21.
//


#ifndef ARITHMETIC_MODULATOR_PY_PROBABILITYEVAL_H
#define ARITHMETIC_MODULATOR_PY_PROBABILITYEVAL_H

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "single_include/robin_hood.h"
#include "string"
#include "vector"
#include "optional"
#include <mutex>

using namespace std;

typedef robin_hood::unordered_map<char, double> char2double;
typedef robin_hood::unordered_map<string, char2double> string2char2double;
typedef robin_hood::unordered_map<string, double> string2double;
typedef robin_hood::unordered_map<string, int32_t> string2int32;
typedef robin_hood::unordered_map<string, vector<string>> string2stringVec;

class ProbabilityEval {

protected:
    bool normalized;
    string2char2double updatedDict;
    void loop();
    void update_cb();
    void addBase(char newBase);
    void removeMotifs();

public:
    ProbabilityEval(string &seq, string2stringVec& motif, int start,
                    bool normalized,int codewordLen, robin_hood::unordered_map<string, char2double> *probDict);

    int codewordLen;
    string2stringVec motif;
    string seq;
    int length;
    int start;
    robin_hood::unordered_map<string, char2double> *probDict;
    char2double nextProbsSingleLetter();
    string2double nextProbsSequence();
    double calcCumProbs();

};

class ProbMap {
    private:
        int codewordLen;
        bool normalized;
        robin_hood::unordered_set<string> codewords;
        string2stringVec motif;
    robin_hood::unordered_map<string, char2double> calcFirstBase(string2int32 &fd_map) const;
        static int32_t checkFd(string ele, char base, const string2int32& fd_map);

    public:
        ProbMap(int cwL, bool norm, robin_hood::unordered_set<string>& cwds, string2stringVec& mot);
        robin_hood::unordered_map<string, char2double> createTransitionDict(string2int32 &freq_dict);
        string2int32 freqDict();
};

#endif //ARITHMETIC_MODULATOR_PY_PROBABILITYEVAL_H
