//
// Created by wintermute on 9/13/21.
//

#ifndef ARITHMETIC_MODULATOR_ERROR_CORRECTION_FREQTABLE_H
#define ARITHMETIC_MODULATOR_ERROR_CORRECTION_FREQTABLE_H

#include "ProbabilityEval.h"
#include <single_include/nlohmann/json.hpp>
#include <fstream>
#include <iostream>

using namespace std;
class FreqTable: public ProbabilityEval {
    private:
        uint16_t wordEnd;

        void setTotal(double csum);
        void initCum();
        void checkSym(int symbol) const;
        void setFreqs(robin_hood::unordered_map<string, char2double> &usedMap);


    public:
        FreqTable(string2stringVec& motif, string seq,
                           int start, bool normalized, int codewordLen, robin_hood::unordered_map<string, char2double>& pDict);
        vector<double> cumulative;
        int total{};
        int start;
        vector<double> frequencies;
        uint8_t codewordLen;

        [[nodiscard]] uint32_t getSymLimit() const;
        double get(int symbol);
        void calcNewFreqs(char newSym);
        void calcFreqs();
        [[nodiscard]] uint64_t getTotal() const;
        [[nodiscard]] uint64_t getLow(int symbol);
        [[nodiscard]] uint64_t getHigh(int symbol);
};



#endif //ARITHMETIC_MODULATOR_ERROR_CORRECTION_FREQTABLE_H
