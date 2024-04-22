from collections import defaultdict
import json
from subprocess import Popen, PIPE, STDOUT
import argparse


def read_fasta(seq, file_name):
    with open(file_name, "r") as file:
        for line in file:
            if not str.startswith(line, ">"):
                seq.append(line.strip().upper())
    return seq

def string_splitting(seq):
    str_map = dict()
    str_map['motif'] = defaultdict(list)
    for entry in seq:
        for i in range(1, len(entry)):
            start = entry[0:i]
            end = entry[i:len(entry)]
            str_map['motif'][start].append(end)
    return str_map

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate a codebook and the concatenation scheme using ContrainedKaos.')
    parser.add_argument('--Path', '-p', dest='path', type=str, action='store',
                        help='path to the ContrainedKaos executable.', default="./ConstrainedKaos/ConstrainedKaos.jar")
    parser.add_argument('--GC', '-GC', dest='gc', type=float, action='store',
                        help='GC content as float, only to be used if a fixed GC content is needed.')
    parser.add_argument('--GClow', '-GCl', dest='gcl', type=float, action='store',
                        help='Minimal GC content of each code word. Use in conjunction with --GChigh')
    parser.add_argument('--GChigh', '-GCh', dest='gch', type=float, action='store',
                        help='Maximal GC content of each code word. Use in conjunction with --GClow.')
    parser.add_argument('--Homopolymer', '-hp', dest='hp', type=int, action='store',
                        help='Homopolymer length constraint, -hp 4 does not allow homopolymers of length 4 or greater.')
    parser.add_argument('--Motifs', '-m', dest='motifs', type=str, action='store',
                        help='Path to a FASTA file containing undesired motifs.')
    parser.add_argument('--Output', '-o', dest='outp', type=str, action='store',
                        help='Destination to save the codewords.', required=True)
    parser.add_argument('--Length', '-l', dest='cw_length', type=int, action='store',
                        help='Length of the output code words, over 10 requires a large amount of memory.', required=True)
    args = parser.parse_args()

    if args.gc and (args.gcl or args.gch):
        parser.error('Only use either --GC (hard constraint) or the interval --GCLow and --GCHigh (soft constraint).')

    seq = []
    if args.motifs:
        seq = read_fasta(seq, args.motifs)
    if args.hp:
        seq.append("A"*args.hp)
        seq.append("G"*args.hp)
        seq.append("T"*args.hp)
        seq.append("C"*args.hp)

    str_map = string_splitting(seq)

    with open(args.outp[:-5] + 'json', "w") as outfile:
        json.dump(str_map, outfile)

    ck_command = ['java', '-jar', args.path, '-length', str(args.cw_length), '-output', str(args.outp)]
    if args.gc:
        ck_command.append('-gc')
        ck_command.append(str(args.gc))
    if args.gcl:
        ck_command.append('-gcStart')
        ck_command.append(str(args.gcl))
    if args.gch:
        ck_command.append('-gcEnd')
        ck_command.append(str(args.gch))
    if args.hp:
        ck_command.append('-hp')
        ck_command.append(str(args.hp))
    if args.motifs:
        ck_command.append('-input')
        ck_command.append(str(args.motifs))


    lines = Popen(ck_command, stdout=PIPE, stderr=STDOUT)
    for line in lines.stdout:
        print(line.decode())
