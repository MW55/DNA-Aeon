//
// Created by wintermute on 9/20/21.
//

#include <utility>
#include "include/bitIO.h"
#include "include/debug_log.h"

using namespace std;

BitInStream::BitInStream(istream &inStream, int sync) :
        in(inStream),
        cByte(0),
        bCount(0),
        sync(sync),
        flag(true),
        streamDone(false),
        remBits(0) {}

int BitInStream::read() {
    if (cByte == -1 || streamDone) {
        return endStream();
    }
    if (!remBits) {
        if (sync != 0 && (bCount % sync == 0) && !flag) {
            unsigned char temp = CRC::Calculate(crcBuffer.data(), crcBuffer.size() * sizeof(char), CRC::CRC_8_LTE());
            cByte = temp;
            crcBuffer.clear();
            remBits = 8;
            flag = true;
        } else {
            int temp = in.get();
            if (temp == -1) {
                return endStream();
            } else {
                crcBuffer.push_back(static_cast<unsigned char>(temp));
                inpBuffer.push_back(static_cast<unsigned char>(temp));
            }
            cByte = temp;
            remBits = 8;
            bCount++;
            flag = false;
        }
    }
    //assert remBits > 0
    remBits--;
    return (cByte >> remBits) & 1;
}

int BitInStream::endStream() {
    if (sync == 0) {
        streamDone = true;
    }
    if (!streamDone && !remBits) {
        inpBuffer.push_back(static_cast<unsigned char>(bCount));
        cByte = CRC::Calculate(inpBuffer.data(), inpBuffer.size() * sizeof(char), CRC::CRC_8());
        // set cByte to 1 if it would be 0, as it would lead to decoding problems.
        if (!cByte)
            cByte = 1;
        remBits = 7;
        streamDone = true;
        return (cByte >> remBits) & 1;
    } else if (streamDone) {
        remBits--;
        if (remBits < 0)
            return -1;
        else {
            return (cByte >> remBits) & 1;
        }
    } else {
        throw logic_error("Leftover bits after encoding.");
    }
}

void BitInStream::close() {
    cByte = -1;
    remBits = 0;
}

BitOutStream::BitOutStream(DecodedData out, int sync) :
        out(move(out)),
        cByte(static_cast<unsigned char>(0)),
        sync(sync),
        crc(8),
        maxBuffer(crc / 8),
        numCrcBytes(0),
        failedSync(false),
        bCount(0),
        crcFlag(false),
        nBitsDone(0),
        done(false) {}

bool BitOutStream::zeroByteEnd(){
    //remove the trailing 1 added during encoding to not have a zero byte crc
    out.data.pop_back();
    if (CRC::Calculate(out.data.data(), out.data.size() * sizeof(char), CRC::CRC_8()) != 0){
        return false;
    } else {
        // remove the bytecount and set nBitsDone to -1, indicating that the decoding is done
        out.data.pop_back();
        nBitsDone = -1;
        return true;
    }
}

/**
 * @brief this function is
 * 1) we set done to true when finish() is called
 * 2) we check the sync and if it's not 0 and done is true
 * 3) we add the remaining bits to the crcBuffer 
 * 4) then is not empty
 * 
 * 5) we call insert on out.data vector (so it copies the crcBuffer to the end of the out.data vector)
 * 6)
 */

void BitOutStream::endStream(){
    if (sync != 0 && done) {
        for (const auto &item: wBuffer) {
            crcBuffer.push_back(item);
        }
        if (!crcBuffer.empty()) {
            cByte = crcBuffer.back();
            crcBuffer.pop_back();
            crcBuffer.push_back(static_cast<unsigned char>(bCount-1));
        } else {
            crcBuffer.push_back(static_cast<unsigned char>(bCount));
        }

        out.data.insert(out.data.end(), crcBuffer.begin(), crcBuffer.end());
        out.data.push_back(cByte);
        int res = CRC::Calculate(out.data.data(), out.data.size() * sizeof(char), CRC::CRC_8());
        if (res != 0) {
            /*
             * if res is 7 it could be that during encoding the last crc was 0
             * and we changed it to 1, as the trailing zeros are not written out by the AC decoder
             * we then have to check if the last crc was zero using the zeroByteEnd method.
             */
            if (res != 7 || !zeroByteEnd()){
                throw logic_error("Last CRC check failed.");
            }
        } else {
            int to_be_removed = maxBuffer + 1;
            out.data.erase(prev(out.data.end(), to_be_removed), out.data.end());
            nBitsDone = -1;
        }
    }
}

void BitOutStream::write(int bit) {
    if (bit != 0 && bit != 1) {
        throw logic_error("Arg has to be binary.");
    }
    cByte = (cByte << 1) | static_cast<unsigned char>(bit);
    nBitsDone++;
    if (nBitsDone > 8) {
        failedSync = true;
        return;
    }
    if (nBitsDone == 8) {
        failedSync = false;
        if (sync == 0) {
            fillBuffer();
            writeBytes();
        } else if (crcFlag && bCount % sync == 0) {
            writeBytes();
            crcCheck(done);
        } else {
            fillBuffer();
            crcFlag = true;
        }
    }
}

void BitOutStream::crcCheck(bool filled) {
    if (numCrcBytes == maxBuffer) {
        if (CRC::Calculate(crcBuffer.data(), crcBuffer.size() * sizeof(char), CRC::CRC_8_LTE()) != 0) {
            failedSync = true;
        } else if(done){
            crcFlag = false;
            crcBuffer.clear();
            wBuffer.clear();
            numCrcBytes = 0;
            nBitsDone = 0;
        } else if (filled) {
            crcBuffer.clear();
            wBuffer.clear();
            if (sync == 1) {
                numCrcBytes = 0;
                crcFlag = false;
            }
        } else {
            crcBuffer.clear();
            crcFlag = false;
            numCrcBytes = 0;
            nBitsDone = 0;
            fillBuffer();
        }
    } else {
        crcBuffer.push_back(static_cast<unsigned char>(cByte));
        numCrcBytes++;
        cByte = static_cast<unsigned char>(0);
        nBitsDone = 0;
        crcCheck(true);
    }
}

void BitOutStream::fillBuffer() {
    wBuffer.push_back(static_cast<unsigned char>(cByte));
    cByte = static_cast<unsigned char>(0);
    nBitsDone = 0;
    bCount++;
}

void BitOutStream::writeBytes() {
    if (!wBuffer.empty()) {
        for (auto &item: wBuffer) {
            out.put(static_cast<char>(item));
            crcBuffer.push_back(item);
        }
        wBuffer.clear();
    }
}

void BitOutStream::close() {
    while (nBitsDone)
        write(0);
}