from math import floor, log
from bitarray import bitarray
import bitarray.util
import argparse


def encode(filepath, codebook):
    inf_rate = floor(log(len(codebook), 2)) / len(codebook[0])
    print("Information rate: " + str(inf_rate))
    bits_per_codeword = floor(inf_rate*len(codebook[0]))
    with open(filepath, "rb") as f:
        bits = bitarray.bitarray()
        bits.fromfile(f)
        chunks = [bits[i:i+int(bits_per_codeword)] for i in range(0, len(bits), int(bits_per_codeword))]
        encoded_data = list()
        for chunk in chunks:
            encoded_data.append(codebook[bitarray.util.ba2int(chunk)])
        return encoded_data


def decode(filepath, codebook):
    inf_rate = floor(log(len(codebook), 2)) / len(codebook[0])
    print("Information rate: " + str(inf_rate))
    bits_per_codeword = floor(inf_rate * len(codebook[0]))
    bit_arr = bitarray.bitarray()
    seq_indices = dict(zip(codebook, range(0, len(codebook))))
    with open(filepath, 'r') as f:
        indata = ""
        infile = f.readlines()
        for line in infile[1::2]:
            indata += ("".join(line.strip()))
        encoded_data = [indata[i:i + len(codebook[0])] for i in range(0, len(indata), len(codebook[0]))]
    for i in range(len(encoded_data)):
        index = seq_indices[encoded_data[i]]
        if i == len(encoded_data)-1:
            bits = bitarray.util.int2ba(index)
            mod_eigth = (len(bit_arr) + len(bits)) % 8
            zeros_needed = mod_eigth if mod_eigth == 0 else (8 - mod_eigth)
            bits = bitarray.bitarray("0"*zeros_needed) + bits
        else:
            bits = bitarray.util.int2ba(index, bits_per_codeword)
        bit_arr.extend(bits)
    return bit_arr


def read_codebook(codebook_path):
    with open(codebook_path, "r") as f_:
        codebook = []
        lines = f_.readlines()
        for line in lines[1::2]:
            codebook.append(line.strip())
    return sorted(codebook, key=str.upper)


def write(output, path, encoding, seq_len = 0):
    if encoding:
        if seq_len != 0:
            c = 0
            with open(path, "w") as f:
                for i in range(0, len(output), seq_len):
                    seq = "".join(output[i:i + seq_len])
                    f.write(">" + str(c) + "\n")
                    f.write(seq + "\n")
                    c += 1
        else:
            with open(path, "w") as f:
                out = "".join(output)
                f.write(">0\n")
                f.write(out)
    else:
        with open(path, "wb") as f:
            output.tofile(f)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', '-in', help='Path to the input file to encode/decode', required=True, type=str)
    parser.add_argument('--output', '-out', help='Output filepath', required=True, type=str)
    parser.add_argument('--codebook', '-c', help='Path to the code book (FASTA)', required=True, type=str)
    parser.add_argument('--seq_len', '-l', help='sequence length of the output FASTA file (only for encoding), multiple of the codeword length, default: one long DNA string', required=False, type=int, default=0)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--encode', '-e', help='encode file', dest='encoding', action='store_true')
    group.add_argument('--decode', '-d', help='decode file', dest='decoding', action='store_true')
    args = parser.parse_args()
    input_path = args.input
    output_path = args.output
    codebook_path = args.codebook
    len_multi = args.seq_len
    if args.encoding:
        encoding = True
    elif args.decoding:
        encoding = False
    else:
        raise ValueError("Either encoding or decoding has to be set.")

    print("encoding: " + str(encoding))

    codebook = read_codebook(codebook_path)

    if encoding:
        encoded_data = encode(input_path, codebook)
        write(encoded_data, output_path, encoding, len_multi)
    else:
        decoded_data = decode(input_path, codebook)
        write(decoded_data, output_path, encoding, len_multi)

