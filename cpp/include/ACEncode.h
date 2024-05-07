//
// Created by wintermute on 9/24/21.
//

#ifndef ARITHMETIC_MODULATOR_ERROR_CORRECTION_ACENCODE_H
#define ARITHMETIC_MODULATOR_ERROR_CORRECTION_ACENCODE_H

#include "ACBase.h"
#include "bitIO.h"
#include "FreqTable.h"
#include <mutex>
#include <sstream>

using namespace std;

class Inflate : public acBase {
private:
    int nBits;

    void underflow() override;

public:
    explicit Inflate(int nBits, BitInStream &bitin);

    BitInStream &bitin;
    bool readingDone;
    uint64_t code;

    void shift() override;
    int read(FreqTable &freqs);

    int readCB();

    int readCBEOF(BitInStream &input);
};


// you need T to have the put method implemented (using concepts would be better?)
template<typename T>
void inflating(FreqTable &freqs, BitInStream &bitin, T &out, int minLen = 0) {
    Inflate inf(16, bitin);
    bool flag = false;
    int len = 0;
    while (true) {
        int sym = inf.read(freqs);
        if (len >= minLen - 1) {
            if (inf.readingDone) {
                if(bitin.streamDone && (bitin.remBits <= 0)){
                    out.put(char(sym));
                    break;
                }
            }
        }
        out.put(char(sym));
        len++;
    }
}

extern void do_encode(const string &data, uint8_t sync, FreqTable freqs, uint16_t minLen, const string &filename,
                      list<tuple<string, vector<unsigned char>>> &results, std::mutex *res_lock);

#endif //ARITHMETIC_MODULATOR_ERROR_CORRECTION_ACENCODE_H
