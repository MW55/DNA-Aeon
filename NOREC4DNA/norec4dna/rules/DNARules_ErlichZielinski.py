import copy
import math

from norec4dna.rules.RuleParser import longestSequenceOfChar, gc_content


class DNARules_ErlichZielinski:
    def __init__(self, active_rules=None):
        self.active_rules = [
            DNARules_ErlichZielinski.homopolymers,
            DNARules_ErlichZielinski.gc_content,
            DNARules_ErlichZielinski.simple_motif_search,
            # DNARules_ErlichZielinski.windowed_gc_content,
        ]

    @staticmethod
    def homopolymers(data):
        return 1.0 if longestSequenceOfChar(data, "*")[1] > 3 else 0.0

    @staticmethod
    def gc_content(data):
        gcc = gc_content(data)
        return 1.0 if gcc > 55 or gcc < 45 else 0.0

    @staticmethod
    def windowed_gc_content(data, window_size=50):
        chunks = [gc_content(data[i:i + window_size]) for i in range(0, len(data), window_size)]
        m = min(chunks)
        ma = max(chunks)
        res = math.pow(ma - m, 2) / 100 * 0.05
        return min(1.0, res)

    @staticmethod
    def add_complementary(input_lst):
        def revert(in_chr):
            return {"G": "C", "C": "G", "T": "A", "A": "T"}[in_chr]

        tmp = copy.deepcopy(input_lst)

        for x in input_lst:
            tmp.append("".join([revert(chr) for chr in x]))
        return tmp

    @staticmethod
    def add_reverse(input_lst):
        tmp = copy.deepcopy(input_lst)
        for x in input_lst:
            tmp.append(x[::-1])
        return tmp

    @staticmethod
    def simple_motif_search(data):
        return 1.0 if any(
            [x in data for x in DNARules_ErlichZielinski.add_complementary(
                DNARules_ErlichZielinski.add_reverse(["ATAACTTCGTATAGCATACATTATACGAAGTTAT",
                                                      "ATAACTTCGTATAGCATACATTATACGAACGGTA",
                                                      "TACCGTTCGTATAGCATACATTATACGAAGTTAT",
                                                      "TACCGTTCGTATAGCATACATTATACGAACGGTA",
                                                      "TACCGTTCGTATATGGTATTATATACGAAGTTAT",
                                                      "TACCGTTCGTATATTCTATCTTATACGAAGTTAT",
                                                      "TACCGTTCGTATAGGATACTTTATACGAAGTTAT",
                                                      "TACCGTTCGTATATACTATACTATACGAAGTTAT",
                                                      "TACCGTTCGTATACTATAGCCTATACGAAGTTAT",
                                                      "ATAACTTCGTATATGGTATTATATACGAACGGTA",
                                                      "ATAACTTCGTATAGTATACCTTATACGAAGTTAT",
                                                      "ATAACTTCGTATAGTATACATTATACGAAGTTAT",
                                                      "ATAACTTCGTATAGTACACATTATACGAAGTTAT",
                                                      "GCATACAT",
                                                      "TGGTATTA",
                                                      "TTCTATCT",
                                                      "GGATACTT",
                                                      "TACTATAC",
                                                      "CTATAGCC",
                                                      "AGGTATGC",
                                                      "TTGTATGG",
                                                      "GGATAGTA",
                                                      "GTGTATTT",
                                                      "GGTTACGG",
                                                      "TTTTAGGT",
                                                      "GTATACCT",
                                                      "GTACACAT",
                                                      "GAAGAC",
                                                      "CTTCTG",
                                                      "GGTCTC",
                                                      "CCAGAG"]))]) else 0.0

    @staticmethod
    def apply_all_rules(packet):
        dna_data = packet.get_dna_struct(True)
        res_arr = [
            x(dna_data)
            for x in [
                DNARules_ErlichZielinski.homopolymers,
                DNARules_ErlichZielinski.gc_content,
                DNARules_ErlichZielinski.simple_motif_search,
                # DNARules_ErlichZielinski.windowed_gc_content,
            ]
        ]
        return sum(res_arr)


def apply_all_rules_with_data(packet):
    try:
        dna_data = packet.get_dna_struct(True)
    except:
        dna_data = packet
    res_arr = [
        x(dna_data)
        for x in [
            DNARules_ErlichZielinski.homopolymers,
            DNARules_ErlichZielinski.gc_content,
            DNARules_ErlichZielinski.simple_motif_search,
            # DNARules_ErlichZielinski.windowed_gc_content,
        ]
    ]
    return sum(res_arr), res_arr


if __name__ == "__main__":
    print(DNARules_ErlichZielinski.homopolymers("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"))
    print(DNARules_ErlichZielinski.gc_content("A" * 1000))
    print(DNARules_ErlichZielinski.windowed_gc_content(
        "ATTAGCGTATCCAATCAGCTGACACCAAAAAAAAAAAATTAGCGTATCCAATCAGCTGACACCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGCAAGCAAAAAAAAAAAAAAAAAAAAAAATTAGCGTATCCAATCAGCTGACACCAAAAAAAAAAAAAAAAAAAAAAAAAAGCAAGCAAAAAACCATGGTCGAAAAAAAAAAAAAAAAAAGCAAGCAAAAAAAAAAAAAAAAAAAAAAAGCAAGCAAAAAACCATGGTCG"))
