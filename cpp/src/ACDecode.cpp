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

void Deflate::write(FreqTable &freqs, int sym) {
    try {
        update(freqs, sym);
    } catch (...){
        bitout.failedSync = true;
    }
}

void Deflate::finish() {
    if (bitout.failedSync)
        return;
    bitout.done = true;
    bitout.endStream();
}

void Deflate::shift() {
    int bit = static_cast<int>(low >> (nBits - 1));
    bitout.write(bit);
    for (int i = 0; i < nUnderflow; ++i) {
        bitout.write(bit^1);
    }
    nUnderflow = 0;
}

void Deflate::underflow() {
    nUnderflow++;
}

[[maybe_unused]] void deflating(FreqTable &freqs, istream& inStream, BitOutStream out){
    Deflate def = Deflate(16, move(out));
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
