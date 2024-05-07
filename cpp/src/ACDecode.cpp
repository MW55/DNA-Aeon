//
// Created by wintermute on 9/27/21.
//

#include "include/ACDecode.h"
#include <utility>

using namespace std;

Deflate::Deflate(int nBits, BitOutStream bitout):
        acBase(nBits),
        nBits(nBits),
        code(0),
        buffer(),
        nUnderflow(0),
        bitout(move(bitout))
{}


//---
void Deflate::write(FreqTable &freqs, int sym) {
    try {
        update(freqs, sym);
    } catch (...){
        bitout.failedSync = true;
    }
}

/**
 * @brief this function is called in the write function
 * we just check that bitout.failedSync didn't throw an exception which would set it to 1
 * if it's not set, we call finish on the bitout stream
 * 
 */
void Deflate::finish() {
    if (bitout.failedSync)
        return;
    bitout.done = true;
    bitout.endStream();
}

/**
 * @brief this shifts function called in the update function of ACBase.cpp
 * 1) bit is the most significant bit of low
 * 2) write bit to the bitout stream
 * 3) for loop from 0 to nUnderflow
 *     write bit^1 to the bitout stream (XOR^1 is the same as NOT)
 * 
 */

void Deflate::shift() {
    int bit = static_cast<int>(low >> (nBits - 1));
    bitout.write(bit);
    for (int i = 0; i < nUnderflow; ++i) {
        bitout.write(bit^1);
    }
    nUnderflow = 0;
}

/**
 * @brief this underflow function called in the update function of ACBase.cpp
 * 1) increment nUnderflow
 * 
 */
void Deflate::underflow() {
    nUnderflow++;
}

/**
 * @brief the function deflating 
 * 1) creates a Deflate object with 16 bits and the output stream
 * 2) test if sym is not 0, then calculate the new frequencies
 *    3) sym is the next character in the input stream
 *   4) if sym is -1 (means EOF ?), finish the Deflate object ()
 * 5) write the frequencies and the symbol to the Deflate object
 * @param freqs 
 * @param inStream 
 * @param out 
 */

[[maybe_unused]] void deflating(FreqTable &freqs, istream& inStream, BitOutStream out){
    Deflate def = Deflate(16, std::move(out));
    int sym = 0;
    while (true){
        if (sym)
            freqs.calcNewFreqs(static_cast<char>(sym));
        sym = inStream.get();
        if(sym == -1){
            def.finish();
            break;
        }
        def.write(freqs, sym);
    }
}
