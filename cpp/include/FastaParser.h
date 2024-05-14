//
// Created by michael on 22.11.21.
//

#ifndef ARITHMETIC_MODULATOR_ERROR_CORRECTION_FASTAPARSER_H
#define ARITHMETIC_MODULATOR_ERROR_CORRECTION_FASTAPARSER_H

#include <iostream>
#include <fstream>
#include "include/debug_log.h"

/**
 * @brief 
 * 
 * you open a file as a output stream and write the sequences to the file
 * you write that way :
 * - for each string, you have a vector of unsigned char
 * - you transform the vector of unsigned char to a "string"
 * - you write the string part of the tuple to the file (>8416_RU10)
 * - you write the vector of unsigned char to the file (casted to a string)
 * 
 * - once done you flush the file and close it
 * 
 * @param file_name 
 * @param list_of_sequences 
 * @return int 
 */

int write_to_fasta(const std::string &file_name, list<tuple<std::string, vector<unsigned char>>> &list_of_sequences) 
{
    std::ofstream fasta_file;
    fasta_file.open(file_name);
    for (tuple<std::string, vector<unsigned char>> &sequence: list_of_sequences) {
        vector<unsigned char> tmp = std::get<1>(sequence);
        std::string seq(reinterpret_cast<const char *>(tmp.data()), tmp.size());
        fasta_file << ">" << std::get<0>(sequence) << "\n" << seq << "\n";
    }
    fasta_file.flush();
    fasta_file.close();
    return 0;
}


robin_hood::unordered_set<std::string> parseFasta(const std::string &fileName) {
    /*
     * Parses a fasta file and returns a set of all sequences.
     */
    robin_hood::unordered_set<std::string> sequences = robin_hood::unordered_set<std::string>();
    std::ifstream input(fileName);
    if (!input.good()) {
        ERROR("Error opening: " << fileName << ". ");
        throw std::runtime_error("Error opening: " + fileName + " .");
    }
    std::string line, seq_id, DNA_sequence;

    while (std::getline(input, line)) {
        // line may be empty so you *must* ignore blank lines
        // or you have a crash waiting to happen with line[0]
        if (line.empty())
            continue;

        if (line[0] == '>') {
            // output previous line before overwriting seq_id
            // but ONLY if seq_id actually contains something
            if (!seq_id.empty())
                sequences.insert(DNA_sequence);
            seq_id = line.substr(1);
            DNA_sequence.clear();
        } else if (line[0] == ';') {
            //skip comment lines
            continue;
        } else {
            DNA_sequence += line;
        }
    }

    // output final entry ONLY if seq_id actually contains something
    if (!seq_id.empty())
        sequences.insert(DNA_sequence);
    return sequences;
}

#endif //ARITHMETIC_MODULATOR_ERROR_CORRECTION_FASTAPARSER_H
