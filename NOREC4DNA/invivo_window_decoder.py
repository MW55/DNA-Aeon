import copy
import os
import random
from io import BytesIO

from norec4dna import RU10Decoder, get_error_correction_decode
from norec4dna.helper.quaternary2Bin import tranlate_quat_to_byte
from norec4dna.rules.DNARules import DNARules
from norec4dna.rules.FastDNARules import FastDNARules

"""
This code reads a single strand of DNA and decodes it:
- A sliding window of size _PACKET_SEQ_LENGTH_ is used check if the window is a valid packet 
  (adheres to the rules and has a valid Reed-Solomon code)
  - If the window is a valid packet it gets added to the RU10 decoder and the window is shifted by _PACKET_SEQ_LENGTH_
  - else: the window will be moved by one base and the process will be repeated 
- To reduce the risk of having a single wrong packet propagating to the whole result, 
  the decoding process will be repeated after shuffling the parsed packets.  

To use this script for other experiments you need to change the following variables:
additionally you should ensure the correct ruleset and boundary is used 
"""
REED_SOLOMON_PARITY_LENGTH = 3
_error_correction = get_error_correction_decode("reedsolomon", REED_SOLOMON_PARITY_LENGTH)
RULES_DROP_LIMIT = 2.0
RULES = FastDNARules()
NUM_CHUNKS = 70
PACKET_SEQ_LENGTH = 300
USE_HEADER_CHUNK = False
ID_LEN_FORMAT = "H"
CRC_LEN_FORMAT = "I"
# this file should include a SINGLE fasta-entry:
INPUT_FILE = "assembly.fasta"
# OPTIONAL - only required for creating a match-html file for whitebox analysis:
# this file should include all original sequences encoded into DNA (in any order)
GROUND_TRUTH = ".INFILES/merged_file.fasta"


def create_decoder():
    """
    create a fresh decoder instance
    """
    decoder = RU10Decoder(INPUT_FILE, use_headerchunk=USE_HEADER_CHUNK, error_correction=_error_correction,
                          static_number_of_chunks=NUM_CHUNKS)
    decoder.number_of_chunks = NUM_CHUNKS
    decoder.read_all_before_decode = True
    return decoder


def load_fasta(fasta_file):
    """
    Loads fasta file and returns a dictionary of sequences
    """
    fasta_dict = {}
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                seq_name = line.strip().split()[0][1:]
                fasta_dict[seq_name] = ''
            else:
                fasta_dict[seq_name] += line.strip()
    return fasta_dict


def window_parse_packets():
    """
    This function uses a sliding window to find (correct) packets in the input sequnce
    """
    decoder = create_decoder()
    fasta = load_fasta(INPUT_FILE)
    window_start = 0
    packets = []
    correct = 0
    correct_seqs = []
    for k in fasta.keys():
        while (window_start + PACKET_SEQ_LENGTH) <= len(fasta[k]):
            line = fasta.get(k)[window_start:(window_start + PACKET_SEQ_LENGTH)]
            # ensure the packet adheres to the rules
            rule_err = RULES.apply_all_rules(line)
            if rule_err < RULES_DROP_LIMIT:
                new_pack = decoder.parse_raw_packet(BytesIO(tranlate_quat_to_byte(line)).read(),
                                                    crc_len_format=CRC_LEN_FORMAT,
                                                    packet_len_format="",
                                                    number_of_chunks_len_format="",
                                                    id_len_format=ID_LEN_FORMAT)
                # TODO: one could convert the RS-repaired data back to dna and check if it adheres to all rules - this
                # would be less restrictive since packets that violate constraints might still be repairable
                if new_pack is not None and new_pack != "CORRUPT":
                    # packet correct: add to correct_seqs and move the window by PACKET_SEQ_LENGTH
                    correct += 1
                    packets.append(new_pack)
                    correct_seqs.append(line)
                    window_start += PACKET_SEQ_LENGTH
                else:
                    # RS was unable to repair the packet: move the window by one base
                    window_start += 1
            else:
                # packet does not adhere to the rules: move the window by one base
                window_start += 1
    # (optional) write all correctly parsed sequences to a file
    # this will allow using the usual means of decoding the input (e.g. ConfigWorker)
    with open("correct_seqs.fasta", "w") as f:
        for i, seq in enumerate(correct_seqs):
            f.write(f">{i}\n{seq}\n")
    print(f"Correct sequences: {correct}")
    return packets


def main():
    packets = window_parse_packets()
    for i in range(1000):
        # start with a clean decoder
        _decoder = create_decoder()
        try:
            # shuffle the packets to reduce impact of single wrong packet
            random.shuffle(packets)
            for index, pack in enumerate(packets):
                if index != i:
                    _decoder.input_new_packet(pack)
            # if the decoder solved the matrix, the result should be correct
            # optionally, one could sanity check the result and reject the result if it is not correct
            # such a sanity check could be a check-sum at the ent of the file, or a check of the content
            # (e.g. only ASCII characters...) same for the filename stored in the header chunk (if available).
            if _decoder.solve(partial=True):
                _decoder.saveDecodedFile(null_is_terminator=True, print_to_output=True,
                                         return_file_name=True, partial_decoding=True)
                print(f"Finished after {i} runs.")
                return
        except Exception:
            # The decoder will raise various exceptions if the decoding fails or yields an invalid result
            pass


def merge_files(folder):
    # write content of each file in the folder into a single file
    with open(folder + "/merged_file.fasta", 'w') as outfile:
        for filename in os.listdir(folder):
            if filename.endswith(".txt"):
                with open(folder + "/" + filename, 'r') as infile:
                    outfile.write(f"\n>{filename}\n")
                    for line in infile:
                        outfile.write(line)


def create_match_html():
    correct_seqs = 0
    ground_truth = load_fasta(GROUND_TRUTH)
    ground_truth = set(x for x in ground_truth.values())
    in_vivoed = load_fasta(INPUT_FILE)
    in_vivoed_seq = list(in_vivoed.values())[0]
    vivo_html = copy.deepcopy(in_vivoed_seq)
    for k in copy.deepcopy(ground_truth):
        if k in in_vivoed_seq:
            correct_seqs += 1
            vivo_html = vivo_html.replace(k, f"<br/><font color='green'>{k}</font><br/>")
            ground_truth.remove(k)
    print("Correct sequences:", correct_seqs)
    with open("match.html", "w") as f:
        f.write('<font color="red">' + vivo_html.replace("<br/><br/>", "<br/>") + '</font>')
    print(ground_truth)
    return "match.html"


if __name__ == "__main__":
    # main()
    create_match_html()
