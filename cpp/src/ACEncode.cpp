//
// Created by wintermute on 9/24/21.
// AC logic is based on Project Nayuki "Reference arithmetic coding", https://github.com/nayuki/Reference-arithmetic-coding
//

#include "include/ACEncode.h"
#include "include/debug_log.h"

using namespace std;

/**
 * @brief Construct a new Inflate:: Inflate object
 * nBits is fixed to 16
 * 
 * we get tempB from readCB() (nb_bytes |)
 * code is shifted to the left and tempB is added to it
 * (does it work like this ?)
 * 
 * @param nBits 
 * @param bitin 
 */

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

/**
 * @brief 
 * 
 * @param freqs 
 * @return int 
 */

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

/**
 * @brief readCB()
 * this function reads the bit from the BitInStream object
 * if it's -1 (EOF) 
 *  - if nBits is 16, set temp to 1, otherwise set to 0
 *  if nBits == -1, set readingDone to true
 * 
 * return temp
 * @return int you either return the number of bits read or 0 or 1 if EOF (if we are %16) 
 *
 */

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

/**
 * @brief 
 * you open a stringstream
 * the BitInStream is created
 * create an output stringstream
 * 
 * - inflating method ?
 * - create a std::vector<uc>res(out_str.begin(), out_str.end());
 *  
 * 
 * @param data     is a string of data to be encoded (inarch.GetEntry(ent).GetDataAsString())
 * @param sync     from ["general"]["sync"] (Number of bytes between sync (crc) steps)
 * @param freqs     
 * @param minLen   from ["encode"]["minLen"] (Minimum length of encoded data)
 * @param filename is the name (inarch.GetEntry(ent))
 * @param results  a reference to results list
 * @param res_lock a lock to protect results list
 */

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
