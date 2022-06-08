//
// Created by wintermute on 9/22/21.
//

#ifndef ARITHMETIC_MODULATOR_ERROR_CORRECTION_ARITHMETICCODE_H
#define ARITHMETIC_MODULATOR_ERROR_CORRECTION_ARITHMETICCODE_H

#include <cstdint>
#include "FreqTable.h"

using namespace std;

class acBase {
protected:
    uint64_t stateBits;
    uint64_t range;
    uint64_t mask;
    uint64_t low;
    uint64_t high;

    virtual void shift() {};

    virtual void underflow() {};

public:
    virtual ~acBase() = default;

    explicit acBase(int nBits);

    void update(FreqTable &freqs, int symbol);

    uint64_t hr;
    uint64_t qr;
    uint64_t mrange;
    uint64_t max;

};

#endif //ARITHMETIC_MODULATOR_ERROR_CORRECTION_ARITHMETICCODE_H
