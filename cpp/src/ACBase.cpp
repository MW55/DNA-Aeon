//
// Created by wintermute on 9/22/21.
// AC logic is based on Project Nayuki "Reference arithmetic coding", https://github.com/nayuki/Reference-arithmetic-coding
//

#include "include/ACBase.h"

/**
 * @brief AC stands for arithmetic coding
 * Construct a new ac Base::ac Base object
 * stateBits is the number of bits in the state
 * range is 2^stateBits
 * mask is range - 1 (all bits set to 1)
 * low is 0
 * high is mask
 * hr is range / 2
 * qr is range / 4
 * mrange is qr + 2 
 * max is mrange
 * @param nBits 
 */

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

/**
 * @brief this function is called from the write function in ACDecode.cpp
 * 1)
 * 2) tot is total of freqs 
 * 3) symlow => checkSymbol if(cumulative.empty) initCumu
 *      checkSym if (symbol is [0, frequencies.size()])
 * 4) symhigh is same but offset of +1
 * 
 * 5) while loop1 (low ^ high) & hr == 0
 *      means either hr is 0 or every bits are equals ?
 *      shift (virtual function) so it calls shift of ACDecode.cpp
 *      update ?
 * 6) while loop2 (low &~ high & qr) != 0
 *      means qr = 1, low = 1, high = 0 to loop
 *      underflow (virtual function) so it calls underflow of ACDecode.cpp
 * 7)
 * @param freqs 
 * @param symbol 
 */
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
