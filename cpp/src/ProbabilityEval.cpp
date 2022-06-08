//
// Created by wintermute on 9/8/21.
//

#include "include/ProbabilityEval.h"

using namespace std;

extern std::mutex *res_lock;

ProbabilityEval::ProbabilityEval(string &seq, string2stringVec &mot, int srt,
                                 bool norm, int cwL, robin_hood::unordered_map<string, char2double> *pDict) :
        normalized(norm),
        codewordLen(cwL),
        motif(mot),
        seq(seq),
        length(static_cast<int>(seq.length())),
        start(srt),
        probDict(pDict) {}


void ProbabilityEval::loop() {
    removeMotifs();
    start += codewordLen;
}

void ProbabilityEval::update_cb() {
    while ((length - start) >= codewordLen) {
        ProbabilityEval::loop();
    }
}

char2double ProbabilityEval::nextProbsSingleLetter() {
    update_cb();
    if (length % codewordLen == 0) {
        char2double res = (updatedDict.contains("")) ? updatedDict[""] : (*probDict)[""];
        return res;
    } else {
        string seqpart = seq.substr(start);
        if (updatedDict.contains(seqpart)){
            return updatedDict[seqpart];
        } else if (probDict->contains(seqpart)) {
            return (*probDict)[seqpart];
        } else {
            char2double res = {{'A', 0.0},
                               {'C', 0.0},
                               {'G', 0.0},
                               {'T', 0.0}
            };
            return res;
        }
    }
}

string2double ProbabilityEval::nextProbsSequence() {
    update_cb();
    string seqpart = seq.substr(start);
    if (updatedDict.count(seqpart)){
        string2double res;
        for (const auto &[key, value]: updatedDict[seqpart]) {
            res[seq += key] = value;
        }
        return res;
    } else {
        string2double res;
        for (const auto &[key, value]: (*probDict)[seqpart]) {
            res[seq += key] = value;
        }
        return res;
    }
}

void ProbabilityEval::addBase(char newBase) {
    seq += newBase;
    length++;
    update_cb();
}

void ProbabilityEval::removeMotifs() {
    updatedDict.clear();
    int end = start+codewordLen;
    string substring = seq.substr(0, end);
    for (const auto & [key, value] : motif) {
        if (substring.ends_with(key)) {
            for (const auto &ele: motif[key]) {
                uint32_t ele_len = ele.length() - 1;
                string pre_minus_one = ele.substr(0, ele_len);
                res_lock->lock();
                updatedDict[pre_minus_one] = (*probDict)[pre_minus_one];
                res_lock->unlock();
                updatedDict[pre_minus_one][ele[ele_len]] = 0;
                updatedDict[ele] = {{'A', 0},
                                    {'T', 0},
                                    {'C', 0},
                                    {'G', 0}};
            }
        }
    }
}

double ProbabilityEval::calcCumProbs() {
    start = 0;
    double prob = (*probDict)[""][seq[0]];
    for (string::size_type i = 1; i < seq.size(); i++){
        if (i%codewordLen == 0){
            loop();
            if (updatedDict.count("")) {
                prob *= updatedDict[""][seq[i]];
            } else {
                prob *= (*probDict)[""][seq[i]];
            }
        } else {
            string seqpart = seq.substr(start, i - start);
            if (updatedDict.count(seqpart)){
                prob *= updatedDict[seqpart][seq[i]];
            } else {
                prob *= (*probDict)[seqpart][seq[i]];
            }
        }
        if (!static_cast<bool>(prob)) {
            return 0.0;
        }
    }
    return prob;
}


ProbMap::ProbMap(int cwL, bool norm,robin_hood::unordered_set<string>& cwds, string2stringVec& mot):
    codewordLen(cwL),
    normalized(norm),
    codewords(cwds),
    motif(mot)
{}

robin_hood::unordered_map<string, char2double> ProbMap::createTransitionDict(string2int32 &freq_dict) {
    robin_hood::unordered_map<string, char2double> fd_map = calcFirstBase(freq_dict);
    for (const auto &ele: freq_dict) {
        if ((int)ele.first.length() < codewordLen) {
            int32_t A = checkFd(ele.first, 'A', freq_dict);
            int32_t T = checkFd(ele.first, 'T', freq_dict);
            int32_t G = checkFd(ele.first, 'G', freq_dict);
            int32_t C = checkFd(ele.first, 'C', freq_dict);
            if (normalized) {
                double csum = ele.second;//A + T + G + C;
                if (static_cast<bool>(csum)) {
                    fd_map[ele.first] = {{'A', 1.0 * A / csum},
                                         {'T', 1.0 * T / csum},
                                         {'C', 1.0 * C / csum},
                                         {'G', 1.0 * G / csum}};
                } else {
                    fd_map[ele.first] = {{'A', 0.0},
                                         {'T', 0.0},
                                         {'C', 0.0},
                                         {'G', 0.0}};
                }
            } else {
                fd_map[ele.first] = {{'A', A},
                                     {'T', T},
                                     {'C', C},
                                     {'G', G}};
            }
        }
    }
    return fd_map;
}

robin_hood::unordered_map<string, char2double> ProbMap::calcFirstBase(string2int32 &fd_map) const {
    robin_hood::unordered_map<string, char2double> pdict;
    if (normalized){
        double sum_first_base = fd_map["A"] + fd_map["T"] + fd_map["G"] + fd_map["C"];
        pdict[""] = {{'A', fd_map["A"] / sum_first_base},
                     {'T', fd_map["T"] / sum_first_base},
                     {'C', fd_map["C"] / sum_first_base},
                     {'G', fd_map["G"] / sum_first_base}};
    } else {
        pdict[""] = {{'A', fd_map["A"]},
                     {'T', fd_map["T"]},
                     {'C', fd_map["C"]},
                     {'G', fd_map["G"]}};
    }
    return pdict;
}

string2int32 ProbMap::freqDict() {
    string2int32 fq_map;
    for (int i = 0; i < codewordLen; i++) {
        for (const auto &ele: codewords) {
            fq_map[ele.substr(0, i + 1)] += 1;
        }
    }
    return fq_map;
}

int32_t ProbMap::checkFd(string ele, char base, const string2int32 &fd_map) {
    ele += base;
    auto res = fd_map.find(ele);
    if (res == fd_map.end()) {
        return 0;
    } else {
        return res->second;
    }
}
