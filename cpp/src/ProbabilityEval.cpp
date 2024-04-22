//
// Created by wintermute on 9/8/21.
//

#include "include/ProbabilityEval.h"

using namespace std;

extern std::mutex *res_lock;

/**
 * @brief Construct a new Probability Eval:: Probability Eval object
 * 
 * @param seq 
 * @param mot 
 * @param srt 
 * @param norm 
 * @param cwL 
 * @param pDict 
 */

ProbabilityEval::ProbabilityEval(string &seq, string2stringVec &mot, int srt,
                                 bool norm, int cwL, robin_hood::unordered_map<string, char2double> *pDict) :
        normalized(norm),
        codewordLen(cwL),
        motif(mot),
        seq(seq),
        length(static_cast<int>(seq.length())), //it's UB to cast size_t to int
        start(srt),
        probDict(pDict) {}

/**
 * @brief the loop calls removeMotifs
 * and increments the start variable by codewordLen
 */

void ProbabilityEval::loop() {
    removeMotifs();
    start += codewordLen;
}

/**
 * @brief update_cb calls loop
 * until the length of the sequence is greater than or equal to codewordLen
 * (start is increased by codewordLen each time loop is called)
 */
void ProbabilityEval::update_cb() {
    while ((length - start) >= codewordLen) {
        ProbabilityEval::loop();
    }
}

/**
 * @brief this function returns the probability of the next base (A, T, C, G)
 * 1) it calls update_cb 
 * 2) if length is a multiple of codewordlen, it returns the probability of the empty string ?
 * 3) otherwise, it returns the probability of the next base in the sequence
 * 4) if updatedDict contains the next base, it returns the probability of the next base
 * 5) if probDict contains the next base, it returns the probability of the next base
 * 6) otherwise, it returns a char2double with all values set to 0
 * 
 * @return char2double (unordered_map<char, double>)
 */

char2double ProbabilityEval:: nextProbsSingleLetter() {
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
/**
 * @brief this functions returns a string2double (unordered_map<string, double>)
 * 1) it calls update_cb
 * 2) it creates a string seqpart that is the substring of seq from start to the end
 * 3) if updatedDict contains seqpart, it returns updatedDict[seqpart]
 * 
 * @return string2double 
 */
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

/**
 * @brief this function 
 * add newBase to the sequence (concatenates newBase to the end of seq)
 * increments length by 1
 * calls update_cb
 * 
 * @param newBase 
 */
void ProbabilityEval::addBase(char newBase) {
    seq += newBase;
    length++;
    update_cb();
}

/**
 * @brief the removeMotifs function removes motifs from the updatedDict
 * 
 */

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

/**
 * @brief 
 * 
 * @return double 
 */

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
/**
 * @brief Construct a new Prob Map:: Prob Map object
 * 
 * @param cwL 
 * @param norm 
 * @param cwds 
 * @param mot 
 */
ProbMap::ProbMap(int cwL, bool norm,robin_hood::unordered_set<string>& cwds, string2stringVec& mot):
    codewordLen(cwL),
    normalized(norm),
    codewords(cwds),
    motif(mot)
{}

/**
 * @brief this function creates a transition dictionary
 * 1) you call calcFirstBase with the freq_dict
 * 2) you iterate through the freq_dict
 * 3) if the length of the element is less than codewordLen
 * 4) you calculate the frequency of each base in the element
 * 5) if normalized is true, you check if ele.second is not 0 (int32_t)
 * 6) fd_map[ele.first] is set to the probability of each base in the element
 * 7) if normalized is false, fd_map[ele.first] is set to the frequency of each base in the element
 * 
 * @param freq_dict 
 * @return robin_hood::unordered_map<string, char2double> 
 */

string2char2double ProbMap::createTransitionDict(string2int32 &freq_dict) {
    string2char2double fd_map = calcFirstBase(freq_dict);
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

/**
 * @brief a string2int32 is given and you generate a string2char2double
 * 1) if normalized is true, you calculate the sum of each base
 * 2) you create a string2char2double pdict 
 * 3) you add the empty string to pdict with the probability of each base 
 * @param fd_map 
 * @return robin_hood::unordered_map<string, char2double> 
 */

string2char2double ProbMap::calcFirstBase(string2int32 &fd_map) const 
{
    string2char2double pdict;
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

/**
 * @brief this function returns a string2int32 (unordered_map<string, int32_t>)
 * 1) it creates a string2int32 fq_map
 * 2) it iterate through codewordLen (?)
 * 3) for each element in codewords
 * 4) you count the substring of the element from 0 to i+1 of each element in codewords
 * 
 * @return string2int32 
 */

string2int32 ProbMap::freqDict() {
    string2int32 fq_map;
    for (int i = 0; i < codewordLen; i++) {
        for (const auto &ele: codewords) {
            fq_map[ele.substr(0, i + 1)] += 1;
        }
    }
    return fq_map;
}

/**
 * @brief the checkFd function checks if the element is in the fd_map
 * 1) it concatenates the element and the base
 * 2) it searches for the element in the fd_map (umap<string, int32_t>)
 * 3) if the element is not found, it returns 0
 * 4) if found, it returns the int32_t value of the element
 * 
 * @param ele 
 * @param base 
 * @param fd_map 
 * @return int32_t 
 */

int32_t ProbMap::checkFd(string ele, char base, const string2int32 &fd_map) {
    ele += base;
    auto res = fd_map.find(ele);
    if (res == fd_map.end()) {
        return 0;
    } else {
        return res->second;
    }
}
