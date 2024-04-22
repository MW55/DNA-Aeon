//
// Created by wintermute on 9/29/21.
//

#include "include/ECDecoding.h"
#include <utility>
#include "include/debug_log.h"

using namespace std;

/**
 * @brief Construct a new ECdecoding::ECdecoding object
 * 1) set the messageLen to the length of the input string or to the value of the config
 * 2) set the withCwProbs to the value of the withCWProbs parameter
 * 3) set the readSeq to the input string
 * 4) set the endFailed to true
 * 5) set the frequencyMap to the freqs parameter
 * 6) set the probMap to the pMap parameter
 * 7) set the seqFlag to false
 * 8) set the errorProb to the value of the error_probability parameter of the config
 * 9) set the rate to the values of the rate parameter of the config
 * 10) set the fanos to the result of the getFanos function with the errorProb and rate as parameters
 * 
 * 
 * @param inp           read-only input string
 * @param freqs         frequency table
 * @param pMap          map of probabilities
 * @param withCWProbs   a boolean flag
 * @param config        json object with the configuration
 */

ECdecoding::ECdecoding(string inp, FreqTable &freqs, robin_hood::unordered_map<string, char2double> &pMap, bool withCWProbs, nlohmann::json &config) :
        messageLen(config["decode"]["length"] <= 0 ? inp.size() : static_cast<uint32_t>(config["decode"]["length"])),
        withCwProbs(withCWProbs),
        readSeq(std::move(inp)),
        endFailed(true),
        frequencyMap(freqs),
        probMap(pMap),
        seqFlag(false),
        errorProb(config["decode"]["metric"]["fano"]["error_probability"]),
        rate({config["decode"]["metric"]["fano"]["rate"]["low"],config["decode"]["metric"]["fano"]["rate"]["high"]}),
        fanos(getFanos(errorProb, rate)),
        queueCounter(0),
        config(config) {}

double ECdecoding::codewordFunc(double probs) {
    return probs;
}

void ECdecoding::write(char symbol) {
    outp += symbol;
}

char ECdecoding::defaultReturn(string &basicString, uint32_t index) {
    if ((unsigned long) index >= basicString.size())
        return 'N';
    else
        return basicString[index];
}

/**
 * @brief 
 * 
 * @param sequence 
 * @param fanoMetrics 
 * @param nextProbs 
 * @param bases 
 */

void ECdecoding::metricProbs(SeqEntry &sequence, array<double, 2> &fanoMetrics,
                             char2double &nextProbs, array<char, 4> &bases) {

    for (auto &ele: bases) {
        if (static_cast<bool>(nextProbs[ele])) {
            fanoCheck(sequence, nextProbs, ele);
            checkIns(sequence, ele, nextProbs);
        }
    }
    checkDels(sequence);
}



void ECdecoding::updateState(SeqEntry &sequence,
                             bool hit, char2double &nextProbs, char base, int increNum) {
    sequence.ac.write(sequence.freq, (int) base);
    if (sequence.ac.bitout.failedSync) {
        checkpointCheck(sequence);
    } else {
        if (!sequence.ac.bitout.crcFlag && sequence.ac.bitout.bCount) {
            // crcFlag is off (just had a crc check) and the crc check was passed successful.
            if (!crcCheckpoints.count(sequence.seq)) {
                crcCheckpoints.emplace(sequence.seq, std::pair(0, sequence));
                sequence.lastCrc = sequence.seq;
            }
        }
        sequence.freq.calcNewFreqs(base);

        double penalty = 0;
        //int failCount = crcCheckpoints.at(sequence.lastCrc).first;
        if (crcCheckpoints.at(sequence.lastCrc).first)
            penalty = config["decode"]["metric"]["penalties"]["crc"]; //0.6
        if (!hit) {
            penalty += (sequence.seq.size() / static_cast<int>(config["decode"]["metric"]["penalties"]["no_hit"]));
        }

        sequence.metric -= (fanos[hit] + codewordFunc(nextProbs[base]) - penalty);
        sequence.seq += base;
        sequence.pos += increNum;

        // Put the first bases in the crc check hashmap, so that subtrees can get metric penalties
        // if the subtree fails multiple times at the same crc checkpoint. This should
        // reduce the probability of getting stuck in a local optima.
        if (sequence.pos == 1 && !crcCheckpoints.contains(sequence.seq)) {
            crcCheckpoints.emplace(sequence.seq, std::pair(0, sequence));
            sequence.lastCrc = base;
        }
        queueInsert(sequence);
    }
}


void ECdecoding::checkpointCheck(const SeqEntry &sequence) {
    // the threshold is more or less arbitrary, this can be adjusted.
    // potential problem: if the metric of the last checkpoint is too high,
    // and the queue is full, it will be removed quickly.
    // if indelcheck for last checkpoint is true, do not put it in the queue again
    if (!sequence.seq.empty()) {
        if (crcCheckpoints.contains(sequence.lastCrc)) {
            if (crcCheckpoints.at(sequence.lastCrc).first > config["decode"]["threshold"]["checkpoint"]) {
                SeqEntry lastCheckpoint = crcCheckpoints.at(sequence.lastCrc).second;
                // Do the checkpoint check for the last checkpoint too, so that we can move even
                // further behind if necessary.
                checkpointCheck(lastCheckpoint);
                //crcCheckpoints.erase(sequence.lastCrc);
            } else {
                crcCheckpoints.at(sequence.lastCrc).first++;
            }
        }
    }
}

/**
 * @brief 
 * fanoCheck is called for each base in the sequence. ?
 * 
 * 
 * @param sequence 
 * @param nextProbs 
 * @param base 
 */

void ECdecoding::fanoCheck(SeqEntry &sequence, char2double &nextProbs, char base) {
    //switched sequence.pos to sequence.seq.size(), as it would lead to wrong results if indels are present
    if (sequence.seq.size() > messageLen)
        return;
    if (static_cast<bool>(nextProbs[base])) {
        //Input sequence at pos x == candidate sequence at pos x
        if (defaultReturn(readSeq, sequence.pos) == base) {
            SeqEntry newSeq = sequence;
            updateState(newSeq, true, nextProbs, base, 1);
        } else {
            SeqEntry newSeq = sequence;
            updateState(newSeq, false, nextProbs, base, 1);
        }
    }
}

/**
 * @brief 
 * 
 * @param sequence 
 */

void ECdecoding::checkDels(SeqEntry &sequence) {
    for (auto &baseDel: array<char, 4>{'A', 'T', 'C', 'G'}) {
        SeqEntry newSeqDel = sequence;
        char2double nextProbs = calcNextProbs(newSeqDel);
        if (static_cast<bool>(nextProbs[baseDel])) {
            updateState(newSeqDel, false, nextProbs, baseDel, 0);
        }
    }
}

/**
 * @brief 
 * 
 * @param sequence 
 * @param base 
 * @param nextProbs 
 */

void ECdecoding::checkIns(SeqEntry &sequence, char base, char2double &nextProbs) {
    if (defaultReturn(readSeq, sequence.pos + 1) == base) {
        SeqEntry newSeqIns = sequence;
        updateState(newSeqIns, false, nextProbs, base, 2);
    }
}

void ECdecoding::queueInsert(SeqEntry &sequence) {
    queue.insert(std::make_unique<SeqEntry>(sequence));
}

/**
 * @brief QueueCheck consists of the following steps: 
 * use config["decode"]["threshold"]["loop"] as a threshold
 * 
 * 1) create a vector of SeqEntry objects
 * 2) we push back the first element of the queue in the vector
 * 3) we erase the first element of the queue
 * while counter smaller than the threshold and the queue is not empty
 * 4) we check that the sequence size is not equal to the messageLen
 * 5) we push back the first element of the queue in the vector
 * 6) we erase the first element of the queue
 * 7) we increment the counter
 * 
 * 8) for each seqEntry in the seqEvals vector
 * 9) we calculate the nextProbs
 * 10) we calculate the metricProbs
 * 11) if the queue size is bigger than the config["decode"]["queue"]["size"]
 * 12) we increment the queueCounter¨
 * 13) if the queueCounter is bigger than the config["decode"]["queue"]["runs"]
 * 14) we throw a logic_error
 * 15) we create an iterator to the first element of the queue
 * 16) we erase the elements from the queue from the reduce_to position to the end
 * 
 * @param fanoMetrics 
 * @param bases 
 */

void ECdecoding::queueCheck(array<double, 2> &fanoMetrics, array<char, 4> &bases) {
    vector<SeqEntry> seqEvals;
    auto itr = queue.begin();
    seqEvals.push_back(**itr);
    queue.erase(itr);
    int counter = 0;
    while (!queue.empty() && counter < config["decode"]["threshold"]["loop"]) {
        itr = queue.begin();
        if ((**itr).seq.size() == messageLen)
            break;
        seqEvals.push_back(**itr);
        queue.erase(itr);
        counter++;
    }

    for (auto &seqEntry: seqEvals) {
        char2double nextProbs = calcNextProbs(seqEntry);
        metricProbs(seqEntry, fanoMetrics, nextProbs, bases);
    }
    if (queue.size() > config["decode"]["queue"]["size"]) {
        WARN("reached max queue size.");
        queueCounter++;
        if (queueCounter > config["decode"]["queue"]["runs"]) { //8 5
            throw logic_error("maximum queue counter reached.");
        }
        auto delete_itr = queue.begin();
        int reduce_to =  static_cast<int>( static_cast<int>( config["decode"]["queue"]["size"])*static_cast<double>( config["decode"]["queue"]["reduce"]));
        advance(delete_itr, reduce_to);
        queue.erase(delete_itr, queue.end());
    }
}

char2double ECdecoding::calcNextProbs(SeqEntry &bestSequence) {
    if (withCwProbs) {
        char2double nextProbs = ProbabilityEval(bestSequence.seq,
                                                this->frequencyMap.motif, 0, true, this->frequencyMap.codewordLen,
                                                &probMap).nextProbsSingleLetter();
        return nextProbs;
    } else {
        char2double nextProbs = {{'A', 1.0},
                                 {'T', 1.0},
                                 {'C', 1.0},
                                 {'G', 1.0}
        };
        return nextProbs;
    }
}

// ToDo calculate the fanos for each position seperatly, using the occurence probs as errorProb.
array<double, 2> ECdecoding::getFanos(double errorProb, array<int, 2> &rate) {
    double bias = (static_cast<float>(rate[0]) / static_cast<float>(rate[1]));
    array<double, 2> fanos = {(log2(2 * errorProb) - bias), (log2(2 * (1 - errorProb)) - bias)};
    return fanos;
}

/**
 * @brief 
 * 1) generate the bases
 * 2) extract e_rate from the config ?
 * 3) get the fanoMetrics(errorProb, e_rate)
 * 
 * 4) create a ProbabilityEval object 
 * 5) get the nextProbsSingleLetter from the ProbabilityEval object
 * 6) create a DecodedData object
 * 7) create a BitOutStream object
 * 8) create a Deflate object with 16 bits and the BitOutStream object ?
 * 9) create a SeqEntry object with the probability, the sequence, the position, the Deflate object and the frequencyMap
 * 10) copy the sequence to the baseLine object
 * 11) put the baseLine object in the crcCheckpoint map
 * 12) calculate the metricProbs
 * 13) start the mainLoop
 * 14) return the best candidate ?
 * 15) return SeqEntry object
 * 
 * @param codewordLen length of the codeword
 * @param motif       map of motifs (key: motif, value: vector of motifs), which is used to calculate the probabilities
 * @param config      json object with the configuration
 * @return SeqEntry   decoded sequence ?
 */


SeqEntry ECdecoding::decode(int codewordLen, robin_hood::unordered_map<string, vector<string>> &motif, nlohmann::json &config) {
    array<char, 4> bases{'A', 'T', 'C', 'G'};
    int itCount = 0;
    array<int, 2> e_rate = {config["decode"]["metric"]["fano"]["rate"]["low"],config["decode"]["metric"]["fano"]["rate"]["high"]};
    DEBUG("Starting decoding with error rate: " << errorProb << " and rates: " << e_rate[0] << " " << e_rate[1]);
    array<double, 2> fanoMetrics = getFanos(errorProb, e_rate);
    DEBUG("Fano metrics: " << fanoMetrics[0] << " " << fanoMetrics[1]);
    string s;
    ProbabilityEval currSeq = ProbabilityEval(s, motif, 0, true, codewordLen, &probMap);
    char2double nextProbs = currSeq.nextProbsSingleLetter();
    DecodedData dec = DecodedData();
    BitOutStream bitOut = BitOutStream(dec, config["general"]["sync"]);
    Deflate ac = Deflate(16, bitOut);
    SeqEntry sequence = SeqEntry(0.0, currSeq.seq, 0, ac, frequencyMap);
    //put in the empty sequence in the crcCheckpoint map as a last-resort
    SeqEntry baseLine = sequence;
    crcCheckpoints.emplace(baseLine.seq, std::pair(0, baseLine));
    metricProbs(sequence, fanoMetrics, nextProbs, bases);
    try {
        mainLoop(fanoMetrics, bases, itCount);
    } catch (const logic_error &e) {
        WARN(e.what());
        SeqEntry res = **queue.begin();
        res.metric = 1000;
        return res;
    }
    SUCCESS("Finished first run, checking CRC...");
    auto itr = queue.begin();
    SeqEntry res = **itr;
    queue.erase(itr);
    while (endFailed) {
        try {
            endFailed = false;
            res.ac.finish();
            SUCCESS("CRC SUCCESSFUL in iteration: " << itCount << " metric " << res.metric);
            res = checkCandidate(res, config["decode"]["threshold"]["finish"]);
        } catch (const logic_error &e) {
            INFO("CRC FAILED in iteration: " << itCount);
            endFailed = true;
            try {
                mainLoop(fanoMetrics, bases, itCount);
            } catch (const logic_error &e) {
                WARN(e.what());
                res.metric = 1000;
                return res;
            }
            itr = queue.begin();
            res = **itr;
            queue.erase(itr);
            itCount++;
        }
    }
    return res;
}

/**
 * @brief 
 * 1) check if we are not at the end of the queue
 * we queueCheck(fanoMetrics, bases) => how does it update ?
 * if the queue is empty, throw an out_of_range exception
 * 
 * @param fanoMetrics   array of fano metrics
 * @param bases         the 4 bases (A, T, C, G)
 * @param itCount       iteration counter
 */

void ECdecoding::mainLoop(array<double, 2> &fanoMetrics, array<char, 4> &bases, int itCount) {
    // Stopping when the best candidate has the desired length is suboptimal.
    while ((**queue.begin()).seq.size() != messageLen) {
        queueCheck(fanoMetrics, bases);
        itCount++;
        if (queue.empty()) {
            throw out_of_range("Queue is empty.");
        }
    }
}

/**
 * @brief the checkCandidate function
 * 1) push back the candidate in the candidates vector (from ECdecoding class)
 * 2) if the size of the candidates vector is bigger than the threshold
 * 3) get the best candidate by using the min_element function
 * 4) return the best candidate
 * 5) else set the endFailed to true and return the candidate
 * 
 * @param can 
 * @param threshold 
 * @return SeqEntry 
 */

SeqEntry ECdecoding::checkCandidate(SeqEntry &can, unsigned long threshold) {
    candidates.push_back(can);
    if (candidates.size() > threshold) {
        SeqEntry bestCan = *min_element(candidates.begin(), candidates.end());
        DEBUG("Final candidate: " << bestCan.seq);
        return bestCan;
    } else {
        endFailed = true;
        return can;
    }
}

/**
 * @brief this function is called from the main function
 * 1) try to create an ECdecoding
 * 2) call the decode function
 * 3) the dec.metric double is converted to a string
 * 4) the data is extracted from the Deflate object in the SeqEntry object
 * 5) then we aquire the lock for the results (unique_lock of std::mutex)
 * 6) we push back a list of tuples of (string, vector<uc>)
 * 
 * @param inp           read-only input string
 * @param freqs         frequency table
 * @param tMap          map of probabilities
 * @param motif         map of motifs
 * @param config        json object with the configuration
 * @param codewordLen   length of the codeword
 * @param results       list of tuples with the results
 * @param res_lock      mutex for the results
 */

void do_decode(const string &inp,
               FreqTable &freqs,
               robin_hood::unordered_map<string, char2double> &tMap, 
               robin_hood::unordered_map<string, vector<string>> &motif,
               nlohmann::json &config,
               int &codewordLen, 
               list<tuple<string, vector<unsigned char>>> &results,
               std::mutex *res_lock) {
    string metric_str = "ERR";
    vector<unsigned char> data = {};

    try {
        ECdecoding ecDec = ECdecoding(inp, freqs, tMap, true, config);
        SeqEntry dec = ecDec.decode(codewordLen, motif, config);
        metric_str = to_string(dec.metric);
        data = *dec.ac.bitout.get_data();
    }
    catch (...) {
        WARN("No candidate found.");
    }
    {
        //we emplace back even if data is empty ?
        //std::unique_lock<std::mutex> uLock(*res_lock);
        std::scoped_lock lock(*res_lock); //C++17
        results.emplace_back(metric_str, data); //shouldn't I use std::move of std::tuple<string, vector<unsigned char>> ?
        //results.emplace_back(std::move(std::make_tuple(metric_str, data))); ?
    }
}
