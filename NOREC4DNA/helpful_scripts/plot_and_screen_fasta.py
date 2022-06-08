import galois
import matplotlib
import matplotlib.pyplot as plt
from norec4dna.rules.FastDNARules import FastDNARules

from numpy import zeros
import numpy as np
import pandas as pd
import seaborn as sns

from plot.plot_error_prob_all import plot_error_prob_for_all

matplotlib.rcParams["figure.dpi"] = 800


# GF = galois.GF_factory(2, 1)


def gff(x, header_len=20, repair_len=16):
    return (2 * x) / max(0, (header_len + x + repair_len))


def gf_poly_div_mod(poly_dividend, poly_divisor):
    return galois.Poly.divmod(poly_dividend, poly_divisor)


def create_poly_from_array(bool_array):
    return galois.Poly(bool_array, field=GF)


def plot_density_graph(max_len=500):
    plt.plot([gff(x) for x in range(1, max_len)])
    plt.xlabel("Length (nt)")
    plt.ylabel("Density (bits/nt)")
    plt.show()
    plt.close()


def chaos_plt(sequence):
    max_x = 5000
    max_y = 5000

    def midpoint(p, q):
        return 0.5 * (p[0] + q[0]), 0.5 * (p[1] + q[1])

    def select_pos_for_base(base):
        base = base.strip().upper()
        if base == "A":
            return 0, 0
        elif base == "T":
            return max_x, 0
        elif base == "G":
            return 0, max_y
        elif base == "C":
            return max_x, max_y
        else:
            raise RuntimeError("Illegal char: ", base)

    N = len(sequence)
    x = zeros(N)
    y = zeros(N)

    x[0], y[0] = select_pos_for_base(sequence[0])

    for idx, base in enumerate(sequence[1:]):
        k = select_pos_for_base(base)
        x[idx], y[idx] = midpoint(k, (x[idx - 1], y[idx - 1]))

    plt.scatter(x, y, s=1)
    plt.xticks([x for x in range(0, max_x + 500, 500)])
    plt.yticks([x for x in range(0, max_y + 500, 500)])
    plt.grid(True)


def create_fasta_from_fastq(fastq_filename: str):
    out_lines = []
    rules = FastDNARules()
    with open(fastq_filename, "r") as in_file:
        lines = in_file.readlines()[1:][::4]
    for line in lines:
        err = rules.apply_all_rules(line[18:182])
        if err < 1.0:
            out_lines.append(line[18:182])
    with open("R:/out.fasta", "w") as out_file:
        for line in out_lines:
            out_file.write(">todo\n")
            out_file.write(line + "\n")
    return lines


def fix_(in_file, out_file_str):
    correct = []
    rule = FastDNARules()
    with open(in_file, "r") as inf:
        lines = inf.readlines()
    for line in lines[1::2]:
        line = line.strip()
        if len(line) != 164:
            continue
        err_prob = rule.apply_all_rules(line)
        if err_prob < 1.0:
            correct.append((line, err_prob))
    with open(out_file_str, "w") as out_file:
        for line, err_prob in correct:
            out_file.write(f">%s\n" % err_prob)
            out_file.write(line.strip().replace("\n", "") + "\n")
    cleaned = find_dup_ids(out_file_str)
    with open(out_file_str, "w") as out_file:
        for line in cleaned:
            out_file.write(f">abc\n")
            out_file.write(line.strip().replace("\n", "") + "\n")


def screen(truth, to_check):
    in_set = set()
    wrong = []
    correct = []
    with open(truth, "r") as in_file:
        lines = in_file.readlines()
        for line in lines[1::2]:
            in_set.add(line.strip())
    with open(to_check, "r") as check_file:
        lines = check_file.readlines()
        for line in lines[1::2]:
            if not line.strip() in in_set:
                wrong.append(line.strip())
            else:
                correct.append(line.strip())
    return wrong, correct


def matching(wrongs, gt):
    res = []
    import difflib
    for wrong in wrongs:
        match = difflib.get_close_matches(wrong, gt)
        res.append((wrong, match))


def find_dup_ids(inf):
    ids_to_seq = dict()
    known_ids = set()

    dup_ids = set()
    with open(inf, "r") as in_file:
        lines = in_file.readlines()
        for line in lines[1::2]:
            _id = line[:8]
            if _id in known_ids:
                #print(f"Id %s already exists: " % _id)
                #print(line)
                #print(ids_to_seq[_id])
                dup_ids.add(_id)
            else:
                known_ids.add(_id)
                ids_to_seq[_id] = line
    # u_ids = known_ids - dup_ids
    #for _id in dup_ids:
    #    del ids_to_seq[_id]
    print(len(dup_ids))
    return ids_to_seq.values()


if __name__ == "__main__":
    fix_("DR1803S1_A.dereplicated.fasta", "DR1803S1_A.dereplicated_fixed.fasta")
    fix_("DR2403S2_A.dereplicated.fasta", "DR2403S2_A.dereplicated_fixed.fasta")
    fix_("DR2503S3_A.dereplicated.fasta", "DR2503S3_A.dereplicated_fixed.fasta")
    fix_("DR2603S4_A.dereplicated.fasta", "DR2603S4_A.dereplicated_fixed.fasta")

    #find_dup_ids("DR1803S1_A.dereplicated_fixed.fasta")

    wrong, correct = screen("RU10_output_102k.fasta", "DR1803S1_A.dereplicated_fixed.fasta")
    with open("test.fasta", "w") as out_f:
        for line in correct:
            out_f.write(">asdf\n")
            out_f.write(line.strip().replace("\n", "") + "\n")
    # print(wrong)
    # print(len(wrong))
    #
    # print(correct)
    # print(len(correct))
    #
    # print(matching(wrong[:20], correct))
    #
    # create_fasta_from_fastq("R:/Undetermined_S0_L001_R1_001.fastq")
    # GF = galois.GF_factory(2, 1)
    # x = galois.Poly([1, 1, 1, 1, 1, 1, 0, 0, 1], field=GF)
    # data = create_poly_from_array([1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1])
    # generator_poly = create_poly_from_array([1, 1, 0, 1, 1, 0])
    # res, rem = gf_poly_div_mod(data, generator_poly)
    # print(res)
    # print(rem)
    #
    # data -= rem
    # print("----")
    # res, rem = gf_poly_div_mod(data, generator_poly)
    # print(res)
    # print(rem)
    #
    # plot_density_graph(500)
    # # AACCACGCGATCGGCGACAGCGCTAGCTGTCGACGATAGTAGATAGAAAGTAGAGATGAGAGGGGGAATGATAGTAGTAGTAGAGATGACCCC"
    # with open("umr_logo_sw_scaled.png_RU10.fasta", "r") as in_file:  # .OUTFILES/output.fasta
    #     lines = in_file.readlines()
    # for line in lines[1::2]:
    #     chaos_plt(line.strip())
    # # chaos_plt("TCCTATCTTCACCCGTCGTGCGCCATGTAGGTCTATCTGAAGTCCGGCCTGACGTGAACGCACACTCCATCGAAGGCTAACGGCCGTGCTTTCTTACTAGCTGACTAACTGCCGTGAGAACGGTAGAGCTAGCGCACGGAAGTAAGACAGAACCTCACGATCACTTATCTTTCTGAAGTCATCACGCAATGTCTAAATTACAGCCGGCATCGCTAACGTTCGCAATCTCTTTCGGCATTAATAACACGCGTGCGACCGGACGCCAGGGCGCCCGTTCTTCCGATATGCAGGCCGTCTGCATAGC")
    # plt.show()
    #
    # # plot_error_prob_for_all(100)
    # ltt = {0: 14, 1: 15, 2: 14, 3: 12, 4: 14, 5: 15, 6: 17, 7: 15, 8: 14, 9: 24, 10: 19, 11: 18, 12: 19, 13: 14, 14: 15,
    #        15: 27, 16: 16, 17: 22, 18: 14, 19: 17, 20: 23, 21: 11, 22: 13, 23: 14, 24: 17, 25: 15, 26: 11, 27: 20,
    #        28: 17, 29: 19, 30: 18, 31: 16, 32: 21, 33: 17, 34: 20, 35: 16, 36: 21, 37: 17, 38: 15, 39: 18, 40: 11,
    #        41: 18, 42: 17, 43: 14, 44: 16, 45: 17, 46: 11, 47: 17, 48: 22, 49: 19, 50: 14, 51: 15, 52: 18, 53: 18,
    #        54: 18, 55: 14, 56: 13, 57: 11, 58: 20, 59: 14, 60: 19, 61: 8, 62: 22, 63: 23, 64: 8, 65: 13, 66: 15, 67: 15,
    #        68: 10, 69: 19, 70: 12, 71: 22, 72: 24, 73: 21, 74: 21, 75: 22, 76: 15, 77: 16, 78: 13, 79: 14, 80: 15,
    #        81: 20, 82: 17, 83: 15, 84: 12, 85: 18, 86: 10, 87: 18, 88: 25, 89: 15, 90: 10, 91: 14, 92: 14, 93: 11,
    #        94: 20, 95: 17, 96: 17, 97: 13, 98: 19, 99: 13, 100: 17, 101: 14, 102: 17, 103: 11, 104: 21, 105: 12,
    #        106: 15, 107: 18, 108: 17, 109: 20, 110: 16, 111: 11, 112: 20, 113: 21, 114: 15, 115: 16, 116: 16, 117: 16,
    #        118: 11, 119: 26, 120: 16, 121: 15, 122: 11, 123: 18, 124: 15, 125: 25, 126: 22, 127: 20, 128: 16, 129: 14,
    #        130: 15, 131: 22, 132: 18, 133: 15, 134: 15, 135: 16, 136: 19, 137: 17, 138: 13, 139: 17, 140: 23, 141: 12,
    #        142: 11, 143: 13, 144: 18, 145: 20, 146: 15, 147: 19, 148: 21, 149: 12, 150: 18, 151: 20, 152: 16, 153: 22,
    #        154: 16, 155: 15, 156: 18, 157: 18, 158: 15, 159: 19, 160: 20, 161: 17, 162: 24}
    # onlinet = {0: 109, 1: 112, 2: 106, 3: 130, 4: 90, 5: 130, 6: 133, 7: 130, 8: 134, 9: 158, 10: 124, 11: 95, 12: 130,
    #            13: 136, 14: 138, 15: 123, 16: 137, 17: 98, 18: 125, 19: 115, 20: 116, 21: 132, 22: 112, 23: 136,
    #            24: 122, 25: 135, 26: 128, 27: 144, 28: 116, 29: 124, 30: 149, 31: 98, 32: 134, 33: 157, 34: 140,
    #            35: 135, 36: 104, 37: 123, 38: 118, 39: 80, 40: 118, 41: 116, 42: 118, 43: 127, 44: 128, 45: 100,
    #            46: 116, 47: 144, 48: 132, 49: 92, 50: 135, 51: 96, 52: 92, 53: 144, 54: 116, 55: 113, 56: 101, 57: 134,
    #            58: 127, 59: 122, 60: 139, 61: 109, 62: 105, 63: 132, 64: 123, 65: 112, 66: 102, 67: 126, 68: 88, 69: 99,
    #            70: 122, 71: 161, 72: 112, 73: 152, 74: 74, 75: 126, 76: 150, 77: 118, 78: 121, 79: 129, 80: 110,
    #            81: 147, 82: 119, 83: 154, 84: 134, 85: 154, 86: 113, 87: 104, 88: 173, 89: 148, 90: 107, 91: 132,
    #            92: 120, 93: 159, 94: 121, 95: 134, 96: 136, 97: 130, 98: 151, 99: 118, 100: 133, 101: 109, 102: 132,
    #            103: 104, 104: 125, 105: 132, 106: 130, 107: 160, 108: 132, 109: 107, 110: 135, 111: 150, 112: 151,
    #            113: 124, 114: 92, 115: 122, 116: 146, 117: 96, 118: 150, 119: 95, 120: 134, 121: 139, 122: 142,
    #            123: 108, 124: 120, 125: 145, 126: 122, 127: 105, 128: 140, 129: 120, 130: 128, 131: 116, 132: 128,
    #            133: 150, 134: 93, 135: 96, 136: 145, 137: 122, 138: 117, 139: 110, 140: 123, 141: 107, 142: 96,
    #            143: 158, 144: 132, 145: 110, 146: 108, 147: 133, 148: 124, 149: 96, 150: 143, 151: 82, 152: 115,
    #            153: 128, 154: 82, 155: 152, 156: 142, 157: 121, 158: 147, 159: 125, 160: 146, 161: 132, 162: 110,
    #            163: 126, 164: 116, 165: 151, 166: 132, 167: 136, 168: 94, 169: 137, 170: 124, 171: 157, 172: 120,
    #            173: 102, 174: 92, 175: 138, 176: 103, 177: 121, 178: 126, 179: 122, 180: 102, 181: 111, 182: 114,
    #            183: 112, 184: 126, 185: 137, 186: 132, 187: 126, 188: 120, 189: 146, 190: 114, 191: 128, 192: 122,
    #            193: 118, 194: 134, 195: 127}
    # ru10t = {0: 85, 1: 102, 2: 94, 3: 106, 4: 85, 5: 69, 6: 95, 7: 83, 8: 92, 9: 98, 10: 90, 11: 112, 12: 89, 13: 92,
    #          14: 92, 15: 102, 16: 95, 17: 109, 18: 98, 19: 103, 20: 88, 21: 72, 22: 81, 23: 109, 24: 77, 25: 68, 26: 71,
    #          27: 98, 28: 87, 29: 68, 30: 107, 31: 78, 32: 91, 33: 73, 34: 91, 35: 101, 36: 85, 37: 83, 38: 74, 39: 102,
    #          40: 84, 41: 84, 42: 61, 43: 80, 44: 66, 45: 87, 46: 82, 47: 105, 48: 97, 49: 80, 50: 80, 51: 72, 52: 58,
    #          53: 76, 54: 82, 55: 63, 56: 76, 57: 73, 58: 71, 59: 60, 60: 73, 61: 95, 62: 87, 63: 92, 64: 86, 65: 91,
    #          66: 80, 67: 91, 68: 91, 69: 99, 70: 74, 71: 72, 72: 82, 73: 81, 74: 69, 75: 73, 76: 80, 77: 82, 78: 87,
    #          79: 79, 80: 101, 81: 86, 82: 88, 83: 83, 84: 75, 85: 89, 86: 61, 87: 72, 88: 87, 89: 77, 90: 79, 91: 78,
    #          92: 88, 93: 88, 94: 93, 95: 65, 96: 74, 97: 59, 98: 90, 99: 78, 100: 89, 101: 93, 102: 83, 103: 110,
    #          104: 106, 105: 83, 106: 95, 107: 92, 108: 80, 109: 90, 110: 72, 111: 82, 112: 93, 113: 92, 114: 81,
    #          115: 71, 116: 82, 117: 91, 118: 74, 119: 70, 120: 79, 121: 73, 122: 101, 123: 95, 124: 66, 125: 92,
    #          126: 83, 127: 70, 128: 73, 129: 83, 130: 85, 131: 88, 132: 88, 133: 78, 134: 84, 135: 89, 136: 64, 137: 85,
    #          138: 75, 139: 68, 140: 91, 141: 86, 142: 79, 143: 68, 144: 51, 145: 68, 146: 79, 147: 80, 148: 70, 149: 89,
    #          150: 71, 151: 94, 152: 74, 153: 88, 154: 82, 155: 88, 156: 62, 157: 91, 158: 102, 159: 73, 160: 86,
    #          161: 73, 162: 67}
    #
    # a = plt.plot(ltt.values(), color="b")
    # b = plt.plot(onlinet.values(), color="g")
    # c = plt.plot(ru10t.values(), color="r")
    # plt.legend(["LT", "Online", "RU10"])
    # plt.grid(True)
    # plt.ylabel("Number of packets containing the chunks")
    # plt.xlabel("Chunks")
    # plt.savefig("compare.pdf", bbox_inches="tight")
    # plt.show(block=False)
