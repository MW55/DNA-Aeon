//
// Created by wintermute on 9/13/21.
//

#include "include/FreqTable.h"

using namespace std;

FreqTable::FreqTable(string2stringVec& motif, string seq,
                     int start, bool normalized, int codewordLen,
                     robin_hood::unordered_map<string, char2double>& pDict):
                     ProbabilityEval(seq, motif, start, normalized,
                                     codewordLen, &pDict),

    wordEnd(0),
    start(start),
    codewordLen(codewordLen)
    {}


void FreqTable::setTotal(double csum=0.0) {
    (static_cast<bool>(csum) ? this -> total = csum : this -> total = accumulate(frequencies.begin(), frequencies.end(), 0.0));
}

double FreqTable::get(int symbol) {
    checkSym(symbol);
    return frequencies[symbol];
}

uint32_t FreqTable::getSymLimit() const {
    return frequencies.size();
}

void FreqTable::calcNewFreqs(char newSym) {
    addBase(newSym);
    wordEnd++;
    if (wordEnd == codewordLen){
        this -> seq  = "";
        wordEnd = 0;
    }
    calcFreqs();
}

void FreqTable::calcFreqs() {
    vector<double> freqs(256, 0);
    this -> frequencies = freqs;
    if (updatedDict.contains(seq)) {
        setFreqs(updatedDict);
    } else {
        setFreqs(*probDict);
    }
    initCum();
}

void FreqTable::setFreqs(robin_hood::unordered_map<string, char2double> &usedMap){
    vector<double> baseFreqs;
    baseFreqs.reserve(4);
    for (const auto &[key, val] : usedMap[seq]){
        frequencies[int(key)] = val;
    }
    for (const auto &entry : usedMap[seq])
        baseFreqs.push_back(entry.second);
    setTotal(accumulate(baseFreqs.begin(), baseFreqs.end(), 0.0));
}

uint64_t FreqTable::getTotal() const {
    return total;
}

uint64_t FreqTable::getLow(int symbol) {
    checkSym(symbol);
    if (cumulative.empty())
        initCum();
    return cumulative[symbol];
}

uint64_t FreqTable::getHigh(int symbol) {
    checkSym(symbol);
    if (cumulative.empty())
        initCum();
    return cumulative[symbol+1];
}

void FreqTable::checkSym(int symbol) const {
    if (0 <= symbol && symbol < (int)frequencies.size())
        return;
    throw out_of_range("Symbol is out of range.");
}

void FreqTable::initCum() {
    vector<double> cumul;
    cumul.reserve(257);
    double sum = 0;
    cumul.push_back(0);
    for (const auto &freq : frequencies) {
        sum += freq;
        cumul.push_back(sum);
    }
    cumulative = cumul;
}
