import copy
import math
from functools import partial

try:
    import pybloomfilter
except:
    import pybloom as pybloomfilter  # code for windows...

from norec4dna.rules.RuleParser import longestSequenceOfChar, microsatellite, length, strContainsIllegalChars, \
    charCountBiggerEqualThanX, gc_content, strContainsSub, strContainsSubRegex

try:
    from cdnarules import repeatRegion as rRegion
    from cdnarules import smallRepeatRegion as smallrRegion
except:
    print("C Module failed to load, falling back to slow mode")


    def rRegion(data, repeat_length=20):
        """

        :param data:
        :param repeat_length:
        :return:
        """
        for i in range(len(data) - repeat_length):
            subseq = data[i:i + repeat_length]
            if data[i + 1:].find(subseq) >= 0:
                return 1.0
        return 0.0


    def smallrRegion(data, repeat_length=9):
        count = 1
        for i in range(len(data) - repeat_length):
            subseq = data[i:i + repeat_length]
            if data[i + 1:].find(subseq) >= 0:
                count += 1
        return 1.0 if 1.0 * count * repeat_length / len(data) > 0.44 else count * repeat_length / len(data) * 0.5


def gc_error_calculation(gc_percentage):
    return (100 + (175 * gc_percentage) / 6 - (121 * gc_percentage ** 2) / 72 + (gc_percentage ** 3) / 36 - (
            gc_percentage ** 4) / 7200) / 100


def fs_gc_error_calculation(gc_percentage):
    return 1.0 if gc_percentage > 60 or gc_percentage < 40 else 0.0


def ts_gc_error_calculation(gc_percentage):
    return 1.0 if gc_percentage > 70 or gc_percentage < 30 else 0.0


def gc_strict_calculation(gc_percentage):
    return (
                   100 + 49970.8 * gc_percentage - 2582 * gc_percentage ** 2 + 41.6458 * gc_percentage ** 3 - 0.208229 * gc_percentage ** 4) / 100


def strict_homopolymers():
    return [0.0, 0.0, 0.2, 0.5, 0.8, 1.0]


def three_strict_homopolymers():
    return [0.0, 0.0, 0.0, 0.0, 1.0]


def lax_homopolymers():
    res = [0.0, 0.0]
    for x in strict_homopolymers():
        res.append(x)
    return res


class FastDNARules:
    def __init__(self, active_rules=None):
        self.nineteen_mers = pybloomfilter.BloomFilter(50000000, 0.0001)
        self.tmp_nineteen_mers = pybloomfilter.BloomFilter(50000, 0.0001)
        if active_rules is None:
            self.active_rules = [
                # FastDNARules.a_permutation,
                # FastDNARules.t_permutation,
                # FastDNARules.c_permutation,
                # FastDNARules.g_permutation,
                # FastDNARules.dinucleotid_runs,
                # FastDNARules.homopolymers,
                partial(FastDNARules.homopolymers, probs=three_strict_homopolymers()),
                # FastDNARules.overall_gc_content,
                # To change the GC error function:
                partial(FastDNARules.overall_gc_content, calc_func=fs_gc_error_calculation),
                # FastDNARules.windowed_gc_content,
                partial(FastDNARules.windowed_gc_content, calc_func=ts_gc_error_calculation),
                #  FastDNARules.long_strands,
                #  FastDNARules.illegal_symbols,
                # FastDNARules.trinucleotid_runs,
                # FastDNARules.random_permutations,
                FastDNARules.motif_search
                # FastDNARules.motif_regex_search,
                # FastDNARules.repeatRegion,
                # FastDNARules.smallRepeatRegion,
                # self.check_and_add_mers
            ]
        else:
            self.active_rules = active_rules

    def check_and_add_mers(self, data, length=19):
        chunks = [data[i:i + length] for i in range(0, len(data), length)]
        self.tmp_nineteen_mers.clear_all()
        res = 0.0
        for chunk in chunks:
            if chunk in self.nineteen_mers or chunk in self.tmp_nineteen_mers:
                res = 0.7
                break
            else:
                self.tmp_nineteen_mers.add(chunk)
        for chunk in chunks:
            self.nineteen_mers.add(chunk)
        return res

    @staticmethod
    def repeatRegion(data, repeat_length=20):
        return rRegion(data, repeat_length)

    @staticmethod
    def smallRepeatRegion(data, repeat_length=9):
        """
        :param data:
        :param repeat_length:
        :return:
        """
        return smallrRegion(data, repeat_length)

    # @staticmethod
    def apply_all_rules(self, packet):
        """
        Apply all rules to the DNA sequence of the given packet and returns the summed error probability.
        :param packet:
        :return:
        """
        try:
            dna_data = packet.get_dna_struct(True)
        except Exception as ex:
            dna_data = packet
        res_arr = [
            x(dna_data)
            for x in self.active_rules
        ]
        return sum(res_arr)

    # @staticmethod
    def apply_all_rules_with_data(self, packet):
        """

        :param packet:
        :return:
        """
        try:
            dna_data = packet.get_dna_struct(True)
        except:
            dna_data = packet
        res_arr = [
            x(dna_data)
            for x in self.active_rules
        ]
        return sum(res_arr), res_arr, packet

    @staticmethod
    def homopolymers(data, probs=None):
        """
        Calculates the dropchance based on the chance of homopolymers to mutate. Homopolymers are repeats of the same
        single unit and have a higher chance to mutate than regular polymers.
        :param data: The DNA sequence to check for homopolymers.
        :return: The dropchance based on the occurence of homopolymers.
        """
        if probs is None:
            probs = [0.0, 0.0, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.0]
        length = longestSequenceOfChar(data, '*')[1]
        return min(1.0, probs[length] if length < len(probs) else 1.0)

    @staticmethod
    def dinucleotid_runs(data):
        """
        Calculates the dropchance based on the chance of dinucleotid microsatellites to mutate. Microsatellites are
        repeats of two to six units that have a higher chance to mutate than regular polymers.
        The dropchances are additive in this case (@shouldDrop()).
        Function y = 0.00002x^2 - 0.0001x
        :param data: The DNA sequence to check for microsatellites with length 2.
        :return: The dropchance based on the occurence of microsatellites with length 2.
        """
        _len = microsatellite(data, 2)[0]
        return max(min(1.0, 0.00002 * _len ** 2 - 0.0001 * _len), 0.0)

    @staticmethod
    def trinucleotid_runs(data):
        """
        Calculates the dropchance based on the chance of trinucleotid microsatellites to mutate. Microsatellites are
        repeats of two to six units that have a higher chance to mutate than regular polymers.
        The dropchances are additive in this case (@shouldDrop()).
        Function y = 0.00002x^2 - 0.0001x
        :param data: The DNA sequence to check for microsatellites with length 3.
        :return: The dropchance based on the occurence of microsatellites with length 3.
        """
        _len = microsatellite(data, 3)[0]
        return max(min(1.0, 0.00002 * _len ** 2 - 0.0001 * _len), 0.0)

    @staticmethod
    def long_strands(data):
        """
        Calculates the dropchance based on the length of the sequence since longer strands have a higher chance to mutate.
        Function y = 0.00144363 * exp(0.0140971 * x), x > 117
        :param data: The sequence to check for the strand length.
        :return: The dropchance based on the strand length.
        """
        _len = length(data)
        if 117 <= _len:
            dropchance = 0.00144363 * math.exp(0.0140971 * _len)
        else:
            dropchance = 0.0
        return min(dropchance, 0.2)

    @staticmethod
    def random_permutations(data):
        """
        Calculates the dropchance for simulated random mutations in the DNA-Data.
        :param data: The sequence to simulate random mutations for.
        :return: The dropchance based on random mutations.
        """
        return 0.02

    @staticmethod
    def illegal_symbols(data):
        """
        Checks the DNA data for illegal symbols and returns a dropchance of 1.0 if the sequence contains them.
        :param data: The sequence to check for illegal symbols.
        :return: The dropchance based on the occurence of illegal symbols (0.0 or 1.0).
        """
        if strContainsIllegalChars(data, "ACGT"):
            return 1.0
        return 0.0

    @staticmethod
    def ch_cg_permutation(data, ch):
        """
        Calculates the dropchance based on the number of occurrences of a nucleotide, G and C in this case since their
        chance to mutate is generally higher.
        These errors are additiv.
        :param data: The sequence to check for the number of occurences of a nucleotide.
        :param ch: The nucleotide to check for.
        :return: The dropchance based on the number of occurences of the given nucleotide.
        """
        occ = charCountBiggerEqualThanX(data, ch)
        dropchance = math.floor(occ / 20) * 0.002
        return min(dropchance, 0.02)

    @staticmethod
    def ch_at_permutation(data, ch):
        """
        Calculates the dropchance based on the number of occurrences of a nucleotide, A and T in this case since their
        chance to mutate is generally lower.
        These errors are additiv.
        :param data: The sequence to check for the number of occurences of a nucleotide.
        :param ch: The nucleotide to check for.
        :return: The dropchance based on the number of occurences of the given nucleotide.
        """
        occ = charCountBiggerEqualThanX(data, ch)
        dropchance = math.floor(occ / 20) * 0.001
        return min(dropchance, 0.01)

    @staticmethod
    def a_permutation(data):
        """
        Calls @ch_at_permutation for 'A'
        :param data: Sequence
        :return: Dropchance
        """
        return FastDNARules.ch_at_permutation(data, 'A')

    @staticmethod
    def t_permutation(data):
        """
        Calls @ch_at_permutation for 'T'
        :param data: Sequence
        :return: Dropchance
        """
        return FastDNARules.ch_at_permutation(data, "T")

    @staticmethod
    def g_permutation(data):
        """
        Calls @ch_cg_permutation for 'G'
        :param data: Sequence
        :return: Dropchance
        """
        return FastDNARules.ch_cg_permutation(data, "G")

    @staticmethod
    def c_permutation(data):
        """
        Calls @ch_cg_permutation for 'C'
        :param data: Sequence
        :return: Dropchance
        """
        return FastDNARules.ch_cg_permutation(data, "C")

    @staticmethod
    def overall_gc_content(data, calc_func=gc_error_calculation):
        """
        Calculates the dropchance of the sequence based on the GC-content since GC-rich DNA is more stable.
        :param calc_func: function to calculate the error-prob based on the gc-content %
        :param data: The sequence to check for the GC-content.
        :return: The Dropchance based on the GC-content.
        """
        return max(0.0, min(1.0, calc_func(gc_content(data))))

    @staticmethod
    def windowed_gc_content(data, window_size=50, calc_func=gc_error_calculation):
        chunks = [FastDNARules.overall_gc_content(data[i:i + window_size], calc_func=calc_func) for i in
                  range(0, len(data), window_size)]
        return min(1.0, max(chunks))

    @staticmethod
    def motif_search(data):
        """
        Searches for undesired motifs with given error probabilities.
        :param data:
        :return:
        """
        _undes_motifs = [
            ("CTCGTAGACTGCGTACCA", 1.01),
            ("GACGATGAGTCCTGAGTA", 1.01),
            ("CTGTCTCTTATACACATCT", 1.01),
            ("TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG", 1.01),
            ("GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG", 1.01),
        ]

        undes_motifs = [
            # Promoter recognition motif (Euk).
            ("TATAAA", 0.01),
            # Promoter recognition motifs (Prok).
            ("TTGACA", 0.05),
            ("TGTATAATG", 0.05),
            # Polyadenylation signals (Euk).
            ("AATAAA", 0.01),
            ("TTGTGTGTTG", 0.01),
            # Lox sites.
            ("ATAACTTCGTATAGCATACATTATACGAAGTTAT", 1.01),
            ("ATAACTTCGTATAGCATACATTATACGAACGGTA", 1.01),
            ("TACCGTTCGTATAGCATACATTATACGAAGTTAT", 1.01),
            ("TACCGTTCGTATAGCATACATTATACGAACGGTA", 1.01),
            ("TACCGTTCGTATATGGTATTATATACGAAGTTAT", 1.01),
            ("TACCGTTCGTATATTCTATCTTATACGAAGTTAT", 1.01),
            ("TACCGTTCGTATAGGATACTTTATACGAAGTTAT", 1.01),
            ("TACCGTTCGTATATACTATACTATACGAAGTTAT", 1.01),
            ("TACCGTTCGTATACTATAGCCTATACGAAGTTAT", 1.01),
            ("ATAACTTCGTATATGGTATTATATACGAACGGTA", 1.01),
            ("ATAACTTCGTATAGTATACCTTATACGAAGTTAT", 1.01),
            # Lox site spacers not covered by the Lox sites.
            ("AGGTATGC", 1.01),
            ("TTGTATGG", 1.01),
            ("GGATAGTA", 1.01),
            ("GTGTATTT", 1.01),
            ("GGTTACGG", 1.01),
            ("TTTTAGGT", 1.01),
            ("GTACACAT", 1.01),
            # Restriction enzyme recognition motifs.
            # BpiI
            ("GAAGAC", 1.01),
            # inverse BpiI
            ("CTTCTG", 1.01),
            # BsaI
            ("GGTCTC", 1.01),
            # inverse BsaI
            ("CCAGAG", 1.01),

            ("CGTCTC", 0.01),
            ("GCGATG", 0.01),
            ("GCTCTTC", 0.01),
            # Oligo Adapters.
            ("CTCGTAGACTGCGTACCA", 0.01),
            ("GACGATGAGTCCTGAGTA", 0.01),
            # 5' extensions.
            ("GGTTCCACGTAAGCTTCC", 0.01),
            ("GCGATTACCCTGTACACC", 0.01),
            ("GCCAGTACATCAATTGCC", 0.01),
            # Twister Adapters:
            ("GAAGTGCCATTCCGCCTGACCT", 1.0),  # Twister 5' Adapter
            ("AGGCTAGGTGGAGGCTCAGTG", 1.0)  # Twister 3' Adapter
        ]
        # undes_motifs = FastDNARules.add_reverse_complementary(undes_motifs)
        dropchance = 0.0
        for motif in undes_motifs:
            if strContainsSub(data, motif[0]):
                dropchance += motif[1]
        return min(1.0, dropchance)

    @staticmethod
    def motif_regex_search(data):
        """
        Searches for undesired motifs with given error probabilities.
        :param data:
        :return:
        """
        undes_motifs = [
            # Promoter recognition motif (Euk).
            ("CANYYY", 0.01),
            ("ANCCAATCA", 0.01),
            ("KGGGCGGRRY", 0.01),
            ("KRGGCGKRRY", 0.01),
            # Promoter recognition motifs (Prok).
            ("AAAWWTWTTTTNNNAAA", 0.05),
            # Ribosomal binding site (Euk).
            ("RCCACCATGG", 0.05),
            # Ribosomal binding site (Prok).
            ("AGGAGGACAGCTAUG", 0.05),
            # Lox sites.
            ("ATAACTTCGTATAGTAYACATTATACGAAGTTAT", 0.01)
        ]
        dropchance = 0.0
        for motif in undes_motifs:
            if strContainsSubRegex(data, motif[0]):
                dropchance += motif[1]
        return min(1.0, dropchance)

    @staticmethod
    def add_complementary(input_lst):
        def revert(in_chr):
            return {"G": "C", "C": "G", "T": "A", "A": "T"}[in_chr]

        tmp = copy.deepcopy(input_lst)

        for x in input_lst:
            tmp.append("".join([revert(chr) for chr in x]))
        return tmp

    @staticmethod
    def add_reverse_complementary(input_lst):
        def revert(in_chr):
            return {"G": "C", "C": "G", "T": "A", "A": "T"}[in_chr]

        tmp = copy.deepcopy(input_lst)

        for x in input_lst:
            tmp.append(("".join([revert(chr[0]) for chr in x[0][::-1]]), x[1]))
        return tmp

    @staticmethod
    def add_reverse(input_lst):
        tmp = copy.deepcopy(input_lst)
        for x in input_lst:
            tmp.append(x[::-1])
        return tmp

    @staticmethod
    def simple_motif_search(data):
        return 1.0 if any([x in data for x in FastDNARules.add_complementary(["ATAACTTCGTATAGCATACATTATACGAAGTTAT",
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
                                                                              "CCAGAG"])]) else 0.0


if __name__ == "__main__":
    x = FastDNARules()
    x.motif_search("TTGACA")
    # print(x.add_reverse_complementary([("CTCGTAGACTGCGTACCA", 1.01)]))
    # print(x.check_and_add_mers("AAAAGAGAGAGAGAGAGAGCCCCCCCCCCCCCCCCCCCAACAGAGAGAGAGAGAGAG", 19))
    # print(x.check_and_add_mers("CCCCCCCCCCCCCCCCCCC", 19))
    print(x.homopolymers(
        "GGTCTCGCAAGTTACGTGTCTATTTAGCGCGGCATATCACAGCGGCGGTACGCATAACAGTTTACAGGGAAAGTAGATCATCAGGCGTGGCTAGGGAGCGCGTGTCCTCATTTGTTGAGGAGACGCTAAAGCACCCGGGTAGTAAATATCTGAACATGGGGGGG",
        probs=three_strict_homopolymers()))
    print(x.overall_gc_content(
        "GGTCTCGCAAGTTACGTGTCTATTTAGCGCGGCATATCACAGCGGCGGTACGCATAACAGTTTACAGGGAAAGTAGATCATCAGGCGTGGCTAGGGAGCGCGTGTCCTCATTTGTTGAGGAGACGCTAAAGCACCCGGGTAGTAAATATCTGAACATGGGGGGG"))
    # print(x.motif_search("ATGGTACGCAAGTCTACGAG"))
