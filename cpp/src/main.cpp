//
// Created by wintermute on 10/12/21.
//
#include <single_include/argparse.hpp>
#include <absl/synchronization/internal/thread_pool.h>
#include "../include/commons.h"
#include "string"
#include "iostream"
#include "include/debug_log.h"
#include "../include/bitIO.h"
#include "../include/ACEncode.h"
#include "../include/ECDecoding.h"
#include "include/FastaParser.h"
#include "../../libraries/Zippy/library/Zippy/ZipArchive.hpp"


using namespace std;

list<tuple<string, vector<unsigned char>>> results = list<tuple<string, vector<unsigned char>>>();
std::mutex *res_lock;

/**
 * @brief the function is recursively called to encode the input file
 * if sameLength is true, the function will encode the input file with the same length
 * then recall
 * if update_config is true, the function will update the config file with the new length
 * with do_encode, the function will update the results
 * then we clear the results and call encode_zip again
 * 
 * 
 * @param config 
 * @param freqs       
 * @param minLen     (when sameLength is true, the minLen is the maxLen of the previous encoding)
 * @param sameLength (it's a bool from the config file)
 * @param configPath 
 */

void encode_zip(nlohmann::json config, const FreqTable &freqs, uint16_t minLen, bool sameLength, string configPath) {
    int maxLen = 0;
    Zippy::ZipArchive inarch(config["encode"]["input"]);
    std::vector<std::string> ents = inarch.GetEntryNames();
    /*for (auto &ent: ents){
        std::string data = inarch.GetEntry(ent).GetDataAsString();
        std::cout << ent << std::endl;
        std::cout << data << std::endl;
    }*/

    //static_cast of an unsigned long to an int can be UB if value is greater than 2^31
    absl::synchronization_internal::ThreadPool pool(static_cast<int>(min(static_cast<unsigned long>(config["general"]["threads"]),ents.size())));
    for (auto &ent: ents) {
        string data = inarch.GetEntry(ent).GetDataAsString();
        pool.Schedule(std::bind(do_encode, data, config["general"]["sync"], freqs, minLen, ent, std::ref(results), std::ref(res_lock)));
    }
    while (results.size() != ents.size()) {
        this_thread::sleep_for(chrono::milliseconds(500));
    }
    //either store as zip or as fasta...
    {
        std::unique_lock<std::mutex> uLock(*res_lock);
        DEBUG("results size 2 : " << results.size());
        if (sameLength) {
            DEBUG("SAME LENGTH");
            for (auto &res: results) {
                int size = get<1>(res).size();
                if (size > maxLen)
                {
                    //DEBUG("Size: " << size); outputs 2054
                    maxLen = size;
                }
            }
            DEBUG("Max length: " << maxLen);
        } else if (config["general"]["as_fasta"]) {
            DEBUG("Writing to fasta");
            write_to_fasta(config["encode"]["output"], results);
        } else {
            DEBUG("Writing to zip");
            write_to_zip(config["encode"]["output"], true, config["general"]["zip"]["most_common_only"], config["general"]["zip"]["decodable_only"], results);
        }
        if (config["encode"]["same_length"] && config["encode"]["update_config"]) {
            config[nlohmann::json::json_pointer("/decode/length")] = minLen;
            std::ofstream confout(configPath);
            confout << config;
            confout.flush();
            confout.close();
        }
    }
    if (sameLength) {
        DEBUG("SAME LENGTH_2");
        results.clear();
        encode_zip(config, freqs, maxLen, false, configPath);
    }
}

nlohmann::json parseValidateConfig(string configPath, bool encode) {
    ifstream cf(configPath);
    nlohmann::json config = nlohmann::json::parse(cf, nullptr, true, true);
    
    try {
        //sync required
        static_cast<void>(config.at(nlohmann::json::json_pointer("/general/sync")));
        //as_fasta default false
        config[nlohmann::json::json_pointer("/general/as_fasta")] = config.value<bool>(nlohmann::json::json_pointer("/general/as_fasta"), false);
        // codewords default codewords\codebook_test.fasta
        config[nlohmann::json::json_pointer("/general/codebook/words")] = config.value<string>(nlohmann::json::json_pointer("/general/codebook/words"), R"(./codewords/codebook_test.fasta)");
        // motifs default codewords\codebook_test.json
        config[nlohmann::json::json_pointer("/general/codebook/motifs")] = config.value<string>(nlohmann::json::json_pointer("/general/codebook/motifs"), R"(./codewords/codebook_test.json)");
        // threads default std::thread::hardware_concurrency()
        config[nlohmann::json::json_pointer("/general/threads")] = config.value<int>(nlohmann::json::json_pointer("/general/threads"), static_cast<int>(std::thread::hardware_concurrency()));
        // zip options
        config[nlohmann::json::json_pointer("/general/zip/most_common_only")] = config.value<bool>(nlohmann::json::json_pointer("/general/zip/most_common_only"), true);
        config[nlohmann::json::json_pointer("/general/zip/decodable_only")] = config.value<bool>(nlohmann::json::json_pointer("/general/zip/decodable_only"), true);
        
        if(encode) {
            //input & output required
            static_cast<void>(config.at(nlohmann::json::json_pointer("/encode/input")));
            static_cast<void>(config.at(nlohmann::json::json_pointer("/encode/output")));
            //min_length default 0
            config[nlohmann::json::json_pointer("/encode/min_length")] = config.value<int>(nlohmann::json::json_pointer("/encode/min_length"), 0);
            //same_length default false
            config[nlohmann::json::json_pointer("/encode/same_length")] = config.value<bool>(nlohmann::json::json_pointer("/encode/same_length"), false);
            //update_config default true (only if not zip or samelength)
            config[nlohmann::json::json_pointer("/encode/update_config")] = config.value<bool>(nlohmann::json::json_pointer("/encode/update_config"), true);
        } else {
            //input & output required
            static_cast<void>(config.at(nlohmann::json::json_pointer("/decode/input")));
            static_cast<void>(config.at(nlohmann::json::json_pointer("/decode/output")));
            //length default 0 (which will be changed to file.size())
            config[nlohmann::json::json_pointer("/decode/length")] = config.value<int>(nlohmann::json::json_pointer("/decode/length"), 0);
            // num additional sequences to check in each loop default 1
            config[nlohmann::json::json_pointer("/decode/threshold/loop")] = config.value<int>(nlohmann::json::json_pointer("/decode/threshold/loop"), 1);
            // num additional sequences to check at the end default 0
            config[nlohmann::json::json_pointer("/decode/threshold/finish")] = config.value<int>(nlohmann::json::json_pointer("/decode/threshold/finish"), 0);
            // crc checkpoint checks threshold default 3
            config[nlohmann::json::json_pointer("/decode/threshold/checkpoint")] = config.value<int>(nlohmann::json::json_pointer("/decode/threshold/checkpoint"), 3);
            // rate
            //low default to sync 
            config[nlohmann::json::json_pointer("/decode/metric/fano/rate/low")] = config.value<int>(nlohmann::json::json_pointer("/decode/metric/fano/rate/low"), static_cast<int>(config["general"]["sync"]));
            //high default to sync +1
            config[nlohmann::json::json_pointer("/decode/metric/fano/rate/high")] = config.value<int>(nlohmann::json::json_pointer("/decode/metric/fano/rate/high"), static_cast<int>(config["general"]["sync"]) + 1);
            //errorProb default 0.05
            config[nlohmann::json::json_pointer("/decode/metric/fano/error_probability")] = config.value<double>(nlohmann::json::json_pointer("/decode/metric/fano/error_probability"), 0.05);
            //panalties
            //crc failed default 0.1
            config[nlohmann::json::json_pointer("/decode/metric/penalties/crc")] = config.value<double>(nlohmann::json::json_pointer("/decode/metric/penalties/crc"), 0.1);
            //no hit divisor default 8
            config[nlohmann::json::json_pointer("/decode/metric/penalties/no_hit")] = config.value<int>(nlohmann::json::json_pointer("/decode/metric/penalties/no_hit"), 8);
            //queue
            //max size default: 200000
            config[nlohmann::json::json_pointer("/decode/queue/size")] = config.value<int>(nlohmann::json::json_pointer("/decode/queue/size"), 200000);
            //runs default 5
            config[nlohmann::json::json_pointer("/decode/queue/runs")] = config.value<int>(nlohmann::json::json_pointer("/decode/queue/runs"), 5);
            //reduce default 0.5
            config[nlohmann::json::json_pointer("/decode/queue/reduce")] = config.value<double>(nlohmann::json::json_pointer("/decode/queue/reduce"), 0.5);   
        }
    } catch (nlohmann::json::exception &err){
        cerr << err.what() << endl;
        exit(EXIT_FAILURE);
    }
    return config;
}

/**
 * @brief main function
 * 1) parse command line arguments
 * 2) check we are either encoding or decoding 
 * 3) parse and validate config file
 * 4) read codewords and motifs
 * 5) create probability map
 * 6) create frequency table
 * 7) transistion dict ?
 * 8) if encoding:
 *    - if input is a zip file, encode each file in the zip
 *    - if input is a file, encode the file
 * 9) if decoding:
 *   - if input is a zip file, decode each file in the zip 
 *   - if input is a fasta file, decode each sequence in the fasta file
 *   - if input is a file, decode the file
 * 
 * @param argc 
 * @param argv 
 * @return int 
 */

int main(int argc, char *argv[]) {
    res_lock = new std::mutex(); //why dynamic allocation?
    argparse::ArgumentParser program("arithmetic_error_correction", "1.0.0");
    program.add_argument("-e", "--encode").default_value(false).implicit_value(true).help("encode a file");
    program.add_argument("-d", "--decode").default_value(false).implicit_value(true).help("decode a file");
    program.add_argument("Config").required().help("config file with program parameters");
    try {
        program.parse_args(argc, argv);
    }
    catch (const runtime_error &err) {
        cerr << err.what() << endl;
        cerr << program << endl;
        exit(EXIT_FAILURE);
    }
    // ensure only en- OR decoding is selected
    bool encode = program.get<bool>("--encode");
    bool decode = program.get<bool>("--decode");
    if (!(encode ^ decode) || program.get<bool>("-h")) {
        cerr << program << endl;
        exit(EXIT_FAILURE);
    }

    string configPath = program.get<string>("Config");
    nlohmann::json config = parseValidateConfig(configPath, encode);
    
    robin_hood::unordered_set<string> codewords = parseFasta(config["general"]["codebook"]["words"]); //hp4_gc40-60.fasta
    robin_hood::unordered_map<string, vector<string>> motif = readConcScheme(config["general"]["codebook"]["motifs"]); //hp4_gc40-60.json"
    
    int codewordLen = getCodewordLen(codewords);
    auto propMap = ProbMap(codewordLen, false, codewords, motif); 
    string2int32 freqDict = propMap.freqDict();
    string2char2double pMap = propMap.createTransitionDict(freqDict);
    results = list<tuple<string, vector<unsigned char>>>();
    FreqTable freqs = FreqTable(motif, "", 0, false, codewordLen, pMap);
    freqs.calcFreqs();
    if (encode) {
        DEBUG("Encoding");
        if (static_cast<string>((std::string)config["encode"]["input"]).ends_with(".zip")) {
            DEBUG("Encoding zip file");
            encode_zip(config, freqs, config["encode"]["min_length"], config["encode"]["same_length"], configPath);
        } else {
            DEBUG("Encoding file (not zip)");
            ifstream inStream((std::string)config["encode"]["input"], ios::binary);
            assert(inStream.good());
            BitInStream bin(inStream, config["general"]["sync"]); //sync =2 for now
            ofstream out((std::string)config["encode"]["output"], ios::binary);
            inflating(freqs, bin, out, config["encode"]["min_length"]); //min_length = 0 for now
            if(config["encode"]["update_config"]){ //true in config.json
                out.flush();
                out.close();
                ifstream outin((std::string)config["encode"]["output"]);
                stringstream buffer;
                buffer << outin.rdbuf();
                string inp = buffer.str();
                int size = inp.size();
                config[nlohmann::json::json_pointer("/decode/length")] = size;
                std::ofstream confout(configPath);
                confout << config;
            }
        }
    } else {
        DEBUG("Decoding");
        robin_hood::unordered_map<string, char2double> tMap = ProbMap(codewordLen, true, codewords,
                                                                motif).createTransitionDict(freqDict);
        array<int, 2> e_rate = {config["decode"]["metric"]["fano"]["rate"]["low"],config["decode"]["metric"]["fano"]["rate"]["high"]};
        if (static_cast<string>(config["decode"]["input"]).ends_with(".zip")) {
            DEBUG("Decoding zip file");
            Zippy::ZipArchive inarch(static_cast<string>(config["decode"]["input"]));
            std::vector<std::string> ents = inarch.GetEntryNames();
            absl::synchronization_internal::ThreadPool pool(config["general"]["threads"]);
            //auto& done = pool.getCompletionFuture();
            std::vector<std::string> inpts;
            for (auto &ent: ents) {
                inpts.push_back(inarch.GetEntry(ent).GetDataAsString());
            }
            for (auto &inp: inpts) {
                //use lambda rather than std::bind
                pool.Schedule(std::bind(do_decode,
                                std::ref(inp), 
                                std::ref(freqs), 
                                std::ref(tMap), 
                                std::ref(motif), 
                                std::ref(config), 
                                std::ref(codewordLen), 
                                std::ref(results), 
                                std::ref(res_lock),
                                std::ref(e_rate)));
            }
            while (results.size() != ents.size()) {
                this_thread::sleep_for(chrono::milliseconds(500));
            }
            //done.wait();
            {
                //std::unique_lock<std::mutex> uLock(*res_lock);
                std::scoped_lock lock(*res_lock);
                write_to_zip(config["decode"]["output"], true, static_cast<bool>(config["general"]["zip"]["most_common_only"]), static_cast<bool>(config["general"]["zip"]["decodable_only"]),results);
            }
        } else if (config["general"]["as_fasta"]) {
            DEBUG("Decoding fasta file");
            // since we are returning a set, we wont have duplicates
            robin_hood::unordered_set<string> dna_lines = parseFasta(config["decode"]["input"]);
            absl::synchronization_internal::ThreadPool pool(config["general"]["threads"]);
            DEBUG(dna_lines.size());
            
            //DEBUG("Starting decoding with error rate: " << errorProb << " and rates: " << e_rate[0] << " " << e_rate[1]);
            //DEBUG("Fano metrics: " << fanoMetrics[0] << " " << fanoMetrics[1]);
            for (const string &line: dna_lines) {    
                pool.Schedule(std::bind(do_decode,  
                              std::ref(line), 
                              std::ref(freqs),
                              std::ref(tMap), 
                              std::ref(motif),
                              std::ref(config),
                              std::ref(codewordLen),
                              std::ref(results),
                              std::ref(res_lock),
                              std::ref(e_rate)));
            }
            //TODO THIS CONDITION SHOULD BE IMPROVED (right now we can not make "results" a set and thus might have duplicates)
            // ideally we would not want any (semi) busy waiting and instead be notified when the threadpool is done
            while (results.size() < dna_lines.size()) {
                this_thread::sleep_for(chrono::milliseconds(1000));
            }
            {
                DEBUG("Writing to zip");
                std::unique_lock<std::mutex> uLock(*res_lock);
                write_to_zip(config["decode"]["output"], true, static_cast<bool>(config["general"]["zip"]["most_common_only"]), static_cast<bool>(config["general"]["zip"]["decodable_only"]),results);
            }
        } else {
            DEBUG("Decoding file");
            ifstream iStream((std::string)config["decode"]["input"]);
            stringstream buffer;
            buffer << iStream.rdbuf();
            string inp = buffer.str();
            thread_local ECdecoding ecDec = ECdecoding(inp, freqs, tMap, true, config);
            SeqEntry dec = ecDec.decode(codewordLen, motif, config, e_rate);
            ofstream out((std::string)config["decode"]["output"], ios::out | ios::binary);
            const vector<unsigned char> str = *dec.ac.bitout.get_data();
            std::string seq(reinterpret_cast<const char *>(str.data()), str.size());
            out << seq;
            out.flush();
            out.close();
        }
    }
    std::flush(std::cout);
    return 0;
}
