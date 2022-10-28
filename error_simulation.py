import math
import os
import sys
import pandas as pd
from numpy import count_nonzero
import numpy as np
import random
import subprocess
from copy import deepcopy
import json
from datetime import datetime
from contextlib import closing
from zipfile import ZipFile
import shutil
import warnings
warnings.filterwarnings("ignore")

NOREC4DNA_BASE_PATH = "/home/wintermute/projects/dna_aeon_review_clean/DNA-Aeon/NOREC4DNA"
DNA_AEON_PATH = "/home/wintermute/projects/dna_aeon_review_clean/DNA-Aeon"
INPUT_DATA = "/home/wintermute/projects/dna_aeon_review_clean/DNA-Aeon/data/D"
FILENAME = "D"
CONFIG = "/home/wintermute/projects/dna_aeon_review_clean/DNA-Aeon/config.json"
ENCODED_FILE = "/home/wintermute/projects/dna_aeon_review_clean/DNA-Aeon/data/encoded.fasta"


def modify_seq(original, pos_sub, pos_ins, pos_del):
    modified = deepcopy(original)
    if pos_sub:
        modified = substitutions(modified, pos_sub)
    if pos_ins:
        modified = insertions(modified, pos_ins)
    if pos_del:
        modified = deletions(modified, pos_del)
    return modified  # , len(pos_sub), len(pos_ins), len(pos_del)


def modify_seqs(seqs, results, num_subs, num_dels, num_ins):
    enc_data_len = sum([len(seq) for seq in seqs])
    one_seq_len = len(seqs[0])  # assuming all seqs have the same length
    all_pos_subs = np.random.choice(enc_data_len, num_subs, replace=False)
    all_pos_ins = np.random.choice(enc_data_len, num_ins, replace=False)
    all_pos_dels = np.random.choice(enc_data_len, num_dels, replace=False)
    mult = 1
    modified_seqs = []
    for seq in seqs:
        curr_low = one_seq_len * (mult - 1)
        curr_high = one_seq_len * mult
        pos_sub = [num % one_seq_len for num in all_pos_subs if (curr_high >= num > curr_low)]
        pos_ins = [num % one_seq_len for num in all_pos_ins if (curr_high >= num > curr_low)]
        pos_del = [num % one_seq_len for num in all_pos_dels if (curr_high >= num > curr_low)]
        modified_seqs.append(modify_seq(seq, pos_sub, pos_ins, pos_del))
        mult += 1
    return modified_seqs


def calculate_num_chunks(file, chunk_size):
    """
    calculate the number of chunks in the file
    :param file:
    :param chunk_size:
    :return:
    """
    file_size = os.path.getsize(file)
    return int(math.ceil(1.0 * file_size / chunk_size))


def mess_data_all_seqs(data, num_errors=0):
    """
    mutate DNA data by adding :num_errors: number of errors
    :param data:
    :param num_errors:
    :return:
    """
    data_length = len(data)
    error_pos = random.sample(range(data_length), num_errors)  # ensure that it will be a different position each time
    data = substitutions(data, error_pos)
    return data


def parse_fasta(filename):
    """
    parse a FASTA file into a list of sequences
    :param filename:
    :return:
    """
    global ac_seqlen

    with open(filename, 'r') as f:
        lines = f.readlines()
    sequences = []
    seq = ''
    for line in lines:
        if line[0] == '>':
            if seq != '':
                sequences.append(seq)
            seq = ''
        else:
            seq += line.strip()
    sequences.append(seq)
    ac_seqlen = len(sequences[0])
    print("Seqlen: " + str(ac_seqlen))
    return sequences


def encode_mutate_decode(file, encoder_function, decoder_function, code_name, num_errors, repeat=1, pre_encoded=False):
    """
    encode a file using a given algorithm, apply mess_data and try to decode it
    results will be stored in a csv/yaml file

    :param code_name: name of the evaluated code
    :param file:
    :param encoder_function: function that take the file as input and return a list of sequences
    :param decoder_function: function that take the list of sequences as input and returns information about the success of the decoding
    :param num_errors:
    :param repeat:
    :return:
    """
    results = []
    # open a file and read to variable
    with open(file, 'rb') as f:
        ground_truth = f.read()
    for i in range(repeat):
        print(i)
        print(datetime.now())
        
        if not pre_encoded:
            encoded_seq = encoder_function(file)
        
        encoded_seq = parse_fasta(ENCODED_FILE)

        res = dict()
        res["encoded_bases"] = len(encoded_seq[0].strip()) * len(
            encoded_seq)
        res["encoded_bytes"] = (res["encoded_bases"] * 2 / 8)
        print("Encoded bases: " + str(res["encoded_bases"]))
        
        res["num_subs"] = int(num_errors["substitutions"]*res["encoded_bases"])
        res["num_dels"] = int(num_errors["deletions"]*res["encoded_bases"])
        res["num_ins"] = int(num_errors["insertions"]*res["encoded_bases"])
        print("subs: " + str(res["num_subs"]) + " dels: " + str(res["num_dels"]) + " ins: " + str(res["num_ins"]))
        res["num_errors"] = sum([res["num_subs"], res["num_dels"], res["num_ins"]])
        res["code"] = code_name
        res["information_bytes"] = os.path.getsize(file)

        res["encoded_bases"] = len(encoded_seq[0].strip()) * len(
            encoded_seq)
        res["encoded_bytes"] = (res["encoded_bases"] * 2 / 8)
        
        global base_err_ratio
        base_err_ratio = res["num_errors"]/res["encoded_bases"]

        mutated_seqs = modify_seqs(encoded_seq, res, res["num_subs"], res["num_dels"],
                                   res["num_ins"])
        success = decoder_function(mutated_seqs, ground_truth)
        res["decoding_success"] = success
        print("Number of Errors: " + str(res["num_errors"]))
        print("Decoding success: " + str(success))
        # results.append(success)
        results.append(res)
    return results


def substitutions(original, pos, base=None):
    modified = deepcopy(original)
    for ele in pos:
        if not base:
            base = random.choice(list({"A", "T", "G", "C"}.difference(original[ele])))
        modified = modified[:ele] + base + modified[ele + 1:]
    return modified


def insertions(original, pos, base=None):
    modified = deepcopy(original)
    shift = 0
    pos.sort()
    for ele in pos:
        if not base:
            base = random.choice(list({"A", "T", "G", "C"}))
        modified = modified[:ele + shift] + base + modified[ele + shift:]
        shift += 1
    return modified


def deletions(original, pos):
    modified = deepcopy(original)
    shift = 0
    pos.sort()
    for ele in pos:
        modified = modified[:ele - shift] + modified[ele - shift + 1:]
        shift += 1
    return modified


def encode_dna_aeon(file):
    py_command = ("python3 " + DNA_AEON_PATH + "/encode.py -c " + CONFIG)
    
    process = subprocess.Popen(py_command.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    return


def decode_dna_aeon(sequences, validation_data):
    
    # write sequences to fasta:
    c = 0
    with open(DNA_AEON_PATH + "/data/mut_encoded.fasta", 'w') as f:
        for sequence in sequences:
            f.write(">" + str(c) + "\n")
            f.write(sequence + "\n")

    
    if os.path.exists(DNA_AEON_PATH + "/decoded.txt.zip"):
        os.remove(DNA_AEON_PATH + "/decoded.txt.zip")

    if os.path.exists(DNA_AEON_PATH + "/data/results/" + FILENAME):
        os.remove(DNA_AEON_PATH + "/data/results/" + FILENAME)


    filename = INPUT_DATA.split("/")[-1]
    py_command = ("python3 " + DNA_AEON_PATH + "/decode.py -c " + CONFIG)
    process = subprocess.Popen(py_command.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    
    if os.path.exists(DNA_AEON_PATH + "/decoded.txt.zip"):
        with closing(ZipFile(DNA_AEON_PATH + "/decoded.txt.zip")) as archive:
            count = len(archive.infolist())
            print("Number of files in zip: " + str(count))

    # validate
    if not os.path.exists(DNA_AEON_PATH + "/data/results/" + FILENAME):
        try:
            shutil.copy(DNA_AEON_PATH + "/decoded.txt.zip", DNA_AEON_PATH + "/data/failed/")
        except:
            print("failed copy.")
        return False
    
    try:
        messcheck = np.fromfile(DNA_AEON_PATH + "/data/results/" + FILENAME, dtype=np.uint8)
        validation_data = np.fromstring(validation_data, dtype=np.uint8)
        validation_data = np.append(validation_data, [0] * (len(messcheck) - len(validation_data)), axis=0)
        badbytes = count_nonzero(validation_data - messcheck)
        print(badbytes)
    except ValueError:
        badbytes = 1
    if badbytes:
        try:
            shutil.copy(DNA_AEON_PATH + "/data/results/" + FILENAME, DNA_AEON_PATH + "/data/failed/")
        except:
            pass
        try:
            shutil.copy(DNA_AEON_PATH + "/decoded.txt.zip", DNA_AEON_PATH + "/data/failed/")
        except:
            pass
    return badbytes == 0

if __name__ == '__main__':
    for i in [0.1]:
        (srate, drate, irate) = i * np.array([0.0238, 0.0082, 0.0039])
        print("multiplier: " + str(i))
        num_errors = {"substitutions": srate, "insertions": irate, "deletions": drate}
        code_name = "DNA-Aeon"
        filename = INPUT_DATA.split("/")[-1] 
        if os.path.exists(DNA_AEON_PATH + "/data/results/" + FILENAME):
            os.remove(DNA_AEON_PATH + "/data/results/" + FILENAME)
        res_list = encode_mutate_decode(INPUT_DATA, encode_dna_aeon,
                                        decode_dna_aeon, code_name, num_errors, repeat=100,
                                        pre_encoded=False)
        results = pd.DataFrame(res_list)
        results.to_csv(code_name + "_s" + str(res_list[0]["num_subs"]) + "_i" + str(res_list[0]["num_ins"]) + "_d" + str(res_list[0]["num_dels"]) + ".csv")
