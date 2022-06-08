//
// Created by wintermute on 9/29/21.
//

#include "include/ECDecoding.h"
#include <utility>
#include "include/debug_log.h"

using namespace std;

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

void ECdecoding::checkDels(SeqEntry &sequence) {
    for (auto &baseDel: array<char, 4>{'A', 'T', 'C', 'G'}) {
        SeqEntry newSeqDel = sequence;
        char2double nextProbs = calcNextProbs(newSeqDel);
        if (static_cast<bool>(nextProbs[baseDel])) {
            updateState(newSeqDel, false, nextProbs, baseDel, 0);
        }
    }
}

void ECdecoding::checkIns(SeqEntry &sequence, char base, char2double &nextProbs) {
    if (defaultReturn(readSeq, sequence.pos + 1) == base) {
        SeqEntry newSeqIns = sequence;
        updateState(newSeqIns, false, nextProbs, base, 2);
    }
}

void ECdecoding::queueInsert(SeqEntry &sequence) {
    queue.insert(std::make_unique<SeqEntry>(sequence));
}

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


SeqEntry ECdecoding::decode(int codewordLen, robin_hood::unordered_map<string, vector<string>> &motif, nlohmann::json &config) {
    array<char, 4> bases{'A', 'T', 'C', 'G'};
    int itCount = 0;
    array<int, 2> e_rate = {config["decode"]["metric"]["fano"]["rate"]["low"],config["decode"]["metric"]["fano"]["rate"]["high"]};
    array<double, 2> fanoMetrics = getFanos(errorProb, e_rate);
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

void do_decode(const string &inp, FreqTable &freqs, robin_hood::unordered_map<string, char2double> &tMap, robin_hood::unordered_map<string, vector<string>> &motif,
               nlohmann::json &config, int &codewordLen, list<tuple<string, vector<unsigned char>>> &results, std::mutex *res_lock) {
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
        std::unique_lock<std::mutex> uLock(*res_lock);
        results.emplace_back(metric_str, data);
    }
}
