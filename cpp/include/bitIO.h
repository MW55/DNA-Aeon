//
// Created by wintermute on 9/20/21.
//

#ifndef ARITHMETIC_MODULATOR_ERROR_CORRECTION_BITIO_H
#define ARITHMETIC_MODULATOR_ERROR_CORRECTION_BITIO_H

#define CRCPP_USE_CPP11
#include "single_include/CRC.h"
#include <istream>
#include <ostream>
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

struct DecodedData {
    vector<unsigned char> data;
    void put(char symbol){
        data.push_back((unsigned char)symbol);
    };
};

class BitInStream {
private:
    istream &in;
    int cByte;
    int bCount;
    int sync;
    bool zeroByte = false;
    vector<unsigned char> inpBuffer;
    vector<unsigned char> crcBuffer;

public:
    explicit BitInStream(istream &inStream, int sync);
    bool flag;
    bool streamDone;
    int remBits;
    int read();
    void close();
    int endStream();
};

//template <typename T>
class BitOutStream {
    private:
        DecodedData out;
        //ostream &out;
        int sync;
        vector<unsigned char> wBuffer;
        vector<unsigned char> crcBuffer;
        int crc;
        int maxBuffer;
        int numCrcBytes;
        unsigned char lastCrc = 0;
        unsigned char cByte;
        void fillBuffer();
        void writeBytes();
        void crcCheck(bool rec);
        bool zeroByteEnd();

public:
        BitOutStream(DecodedData out, int sync);
        bool failedSync;
        int bCount;
        bool crcFlag;
        int nBitsDone;
        bool done;
        void write(int bit);
        void close();
        void endStream();
        [[nodiscard]] const vector<unsigned char> *get_data() {
            return &out.data;
        }
};



#endif //ARITHMETIC_MODULATOR_ERROR_CORRECTION_BITIO_H
