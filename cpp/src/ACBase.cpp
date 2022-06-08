//
// Created by wintermute on 9/22/21.
// AC logic is based on Project Nayuki "Reference arithmetic coding", https://github.com/nayuki/Reference-arithmetic-coding
//

#include "include/ACBase.h"

acBase::acBase(int nBits):
    stateBits(nBits),
    range(1ULL << stateBits),
    mask(range - 1),
    low(0),
    high(mask),
    hr(range >> 1),
    qr(hr >> 1),
    mrange(qr + 2),
    max(mrange)
{}

void acBase::update(FreqTable &freqs, int symbol) {
    uint64_t h = high;
    uint64_t l = low;
    uint64_t ran = h - l + 1;
    uint64_t tot = freqs.getTotal();
    uint64_t symLow = freqs.getLow(symbol);
    uint64_t symHigh = freqs.getHigh(symbol);
    assert(symHigh != symLow);

    uint64_t nl =  l + symLow * ran / tot;
    uint64_t nh =  l + symHigh * ran / tot - 1;

    low = nl;
    high = nh;

    while (((low ^ high) & hr) == 0) {
        shift();
        low = ((low << 1) & mask);
        high = ((high << 1) & mask) | 1;
    }

    while ((low &~ high & qr) != 0) {
        underflow();
        low = (low << 1) ^ hr;
        high = ((high ^ hr) << 1) | hr | 1;
    }

}
