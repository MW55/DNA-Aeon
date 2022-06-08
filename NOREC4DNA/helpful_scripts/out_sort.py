import os
import bisect

directory = "norec4dna_new/"
sort_list = []

for file in os.listdir(directory):
    with open(directory + file, "r") as in_file:
        lines = in_file.readlines()
    for i in range(0, len(lines), 2):
        error_prob, seed = lines[i][1:].split(".")[0].split("_")
        dna_str = lines[i + 1]
        bisect.insort_left(sort_list, (error_prob, seed, dna_str))
with open("out.fasta", "w") as out:
    for error_prob, seed, dna_str in sort_list:
        out.write(">" + str(error_prob) + "_" + str(seed) + "\n" + dna_str)  # + "\n")
