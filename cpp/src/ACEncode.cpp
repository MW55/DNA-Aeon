//
// Created by wintermute on 9/24/21.
// AC logic is based on Project Nayuki "Reference arithmetic coding", https://github.com/nayuki/Reference-arithmetic-coding
//

#include "include/ACEncode.h"
#include "include/debug_log.h"

using namespace std;

Inflate::Inflate(int nBits, BitInStream &bitin):
    acBase(nBits),
    nBits(nBits),
    bitin(bitin),
    readingDone(false),
    code(0)
    {
    for (int i=0; i<nBits; i++) {
        int tempB = readCB();
        code = code << 1 | tempB;
    }
    }

int Inflate::read(FreqTable &freqs) {
    uint64_t total = freqs.getTotal();
    range = high - low + 1;
    assert(code >= low);
    uint64_t offset = code - low;
    assert(((offset + 1) * total - 1) / range >= 0);
    unsigned long val = ((offset + 1) * total - 1) / range;
    assert(val*range / total <= offset);
    assert(val < total);
    int start = 0;
    int end = freqs.getSymLimit();

    while (end - start > 1){
        int mid = (start + end) >> 1;
        if (freqs.getLow(mid) > val)
            end = mid;
        else
            start = mid;
    }
    assert(start + 1 == end);

    int sym = start;
    assert((freqs.getLow(sym) * range / total) <= offset && offset < (freqs.getHigh(sym) * range / total));
    update(freqs, sym);
    freqs.calcNewFreqs((char)sym);
    return sym;
}

void Inflate::shift() {
    int tempB = readCB();
    code = ((code << 1) & mask) | tempB;
}

void Inflate::underflow() {
    int tempB = readCB();
    code = (code & hr) | ((code << 1) & (mask >> 1)) | tempB;
}

int Inflate::readCB() {
    int temp = bitin.read();
    if (temp == -1) {
        (nBits == 16 ? temp = 1 : temp = 0);
        nBits--;
        if (nBits == -1) {
            readingDone = true;
        }
    }
    return temp;
}

void do_encode(const string &data, uint8_t sync, FreqTable freqs, uint16_t minLen, const string &filename,
               list<tuple<string, vector<unsigned char>>> &results, std::mutex *res_lock) {
    stringstream inStream(data);
    assert(inStream.good());
    BitInStream bin(inStream, sync);
    stringstream out;
    inflating(freqs, bin, out, minLen);
    std::string out_str = out.str();
    std::vector<unsigned char> res(out_str.begin(), out_str.end());
    {
        std::unique_lock<std::mutex> uLock(*res_lock);
        results.emplace_back(filename, res);
    }
}