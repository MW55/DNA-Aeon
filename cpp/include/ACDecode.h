//
// Created by wintermute on 9/27/21.
//

#ifndef ARITHMETIC_MODULATOR_ERROR_CORRECTION_ACDECODE_H
#define ARITHMETIC_MODULATOR_ERROR_CORRECTION_ACDECODE_H

#include "bitIO.h"
#include "FreqTable.h"
#include "ACBase.h"

using namespace std;

class Deflate : public acBase {
private:
    int nBits;
    uint16_t code;
    vector<bool> buffer;

    void underflow() override;
    //void writeToBuffer(int Bit);


public:
    Deflate(int nBits, BitOutStream bitout);
    int nUnderflow;
    BitOutStream bitout;


    void write(FreqTable &freqs, int sym);

    void finish();

    void shift() override;
};

[[maybe_unused]] void deflating(FreqTable &freqs, istream &inStream, BitOutStream out);


#endif //ARITHMETIC_MODULATOR_ERROR_CORRECTION_ACDECODE_H
