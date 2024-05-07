//
// Created by wintermute on 9/29/21.
//

#ifndef ARITHMETIC_MODULATOR_ERROR_CORRECTION_ECDECODING_H
#define ARITHMETIC_MODULATOR_ERROR_CORRECTION_ECDECODING_H

#include <mutex>
#include <memory>

#include "include/commons.h"
#include "include/FreqTable.h"
#include "include/ProbabilityEval.h"
#include "include/ACDecode.h"
#include "include/bitIO.h"
#include "absl/container/btree_set.h"
#include "absl/container/internal/btree_container.h"
#include "absl/container/flat_hash_map.h"
#include "single_include/robin_hood.h"

using namespace std;

struct SeqEntry {
    double metric;
    string seq;
    uint32_t pos;
    Deflate ac;
    FreqTable freq;
    string lastCrc;

    bool operator<(const SeqEntry &comp) const {
        return metric < comp.metric;
    }

    SeqEntry(const SeqEntry &other) = default;

    SeqEntry(double metric, string &seq, int pos, Deflate acDeflate, FreqTable &freqTable, string lCrc = "") :
            metric(metric),
            seq(seq),
            pos(pos),
            ac(move(acDeflate)),
            freq(freqTable),
            lastCrc(std::move(lCrc)) { /*DEBUG("SeqEntry constructor");*/ };
};

struct compare_ptr {
    bool operator() (const std::unique_ptr<SeqEntry> &lhs, const std::unique_ptr<SeqEntry> &rhs) const {
       return lhs->metric < rhs->metric;
    }
};

class ECdecoding {

private:
    uint32_t messageLen;
    bool withCwProbs;
    string readSeq;
    string outp;
    bool endFailed;
    FreqTable frequencyMap;
    absl::btree_multiset<std::unique_ptr<SeqEntry>, compare_ptr> queue;
    robin_hood::unordered_map<string, std::pair<int, SeqEntry>> crcCheckpoints;
    vector<SeqEntry> candidates;
    robin_hood::unordered_map<string, char2double> &probMap;
    bool seqFlag;
    double errorProb;
    array<int, 2> rate;
    array<double, 2> fanos;
    int queueCounter;
    nlohmann::json config;

    void checkpointCheck(const SeqEntry &sequence);

    static double codewordFunc(double probs);

    static char defaultReturn(string &basicString, uint32_t index);

    void updateState(SeqEntry &sequence,
                     bool hit, char2double &nextProbs, char base, int increNum);

    void write(char symbol);

    void metricProbs(SeqEntry &sequence,
                     array<double, 2> &fanoMetrics, char2double &nextProbs, array<char, 4> &bases);

    void fanoCheck(SeqEntry &sequence, char2double &nextProbs, char base);

    void checkDels(SeqEntry &sequence);

    void checkIns(SeqEntry &sequence, char base, char2double &nextProbs);

    void queueInsert(SeqEntry &sequence);

    void queueCheck(array<double, 2> &fanoMetrics, array<char, 4> &bases);

    char2double calcNextProbs(SeqEntry &bestSequence);

    void mainLoop(array<double, 2> &fanoMetrics, array<char, 4> &bases, int itCount);

    SeqEntry checkCandidate(SeqEntry &can, unsigned long threshold);

    static array<double, 2> getFanos(double errorProb, array<int, 2> &rate);

    long double calcFano(double &probsBase, long double &probsCurrSeq, bool hit) const;


public:
    ECdecoding(string inp, FreqTable &freqs, robin_hood::unordered_map<string, char2double> &probMap,
               bool withCWProbs, nlohmann::json &config);

    SeqEntry decode(int codewordLen, 
                    robin_hood::unordered_map<string, vector<string>> &conc, 
                    nlohmann::json &config,
                    std::array<int, 2> &e_rate);

};

extern void do_decode(const string &inp, FreqTable &freqs,
                      robin_hood::unordered_map<string, char2double> &tMap, robin_hood::unordered_map<string, vector<string>> &motif,
                      nlohmann::json &config, int &codewordLen, list<tuple<string, vector<unsigned char>>> &results,
                      std::mutex *res_lock,
                      std::array<int, 2> &e_rate);

#endif //ARITHMETIC_MODULATOR_ERROR_CORRECTION_ECDECODING_H
