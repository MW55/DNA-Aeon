import copy

from norec4dna.rules.RuleParser import shouldDropMax, shouldDrop, gc_content

try:
    from cdnarules import repeatRegion as rRegion
    from cdnarules import smallRepeatRegion as smallrRegion
except ImportError as ex:
    print("C Module failed to load, falling back to slow mode")
    from norec4dna.helper.fallback_code import r_region as rRegion
    from norec4dna.helper.fallback_code import small_r_region as smallrRegion


class DNARules:
    def __init__(self, active_rules=None):
        self.active_rules = []

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

    @staticmethod
    def homopolymers(data):
        """
        Calculates the dropchance based on the chance of homopolymers to mutate. Homopolymers are repeats of the same
        single unit and have a higher chance to mutate than regular polymers.
        :param data: The DNA sequence to check for homopolymers.
        :return: The dropchance based on the occurence of homopolymers.
        """
        return shouldDropMax(
            data,
            [
                ("longestSequenceOfChar(*,2)", 0.001),
                ("longestSequenceOfChar(*,3)", 0.005),
                ("longestSequenceOfChar(*,4)", 0.01),
                ("longestSequenceOfChar(*,5)", 0.05),
                ("longestSequenceOfChar(*,6)", 0.1),
                ("longestSequenceOfChar(*,7)", 0.5),
                ("longestSequenceOfChar(*,8)", 1.0),
            ],
        )

    @staticmethod
    def dinucleotid_runs(data):
        """
        Calculates the dropchance based on the chance of dinucleotid microsatellites to mutate. Microsatellites are
        repeats of two to six units that have a higher chance to mutate than regular polymers.
        The dropchances are additive in this case (@shouldDrop()).
        :param data: The DNA sequence to check for microsatellites with length 2.
        :return: The dropchance based on the occurence of microsatellites with length 2.
        """
        return shouldDrop(
            data,
            [
                ("microsatelliteLongerThanX(2,10)", 0.001),
                ("microsatelliteLongerThanX(2,15)", 0.002),
                ("microsatelliteLongerThanX(2,20)", 0.003),
                ("microsatelliteLongerThanX(2,25)", 0.004),
                ("microsatelliteLongerThanX(2,30)", 0.005),
                ("microsatelliteLongerThanX(2,35)", 0.006),
                ("microsatelliteLongerThanX(2,40)", 0.007),
                ("microsatelliteLongerThanX(2,45)", 0.008),
                ("microsatelliteLongerThanX(2,50)", 0.009),
                ("microsatelliteLongerThanX(2,55)", 0.010),
                ("microsatelliteLongerThanX(2,60)", 0.011),
                ("microsatelliteLongerThanX(2,65)", 0.012),
                ("microsatelliteLongerThanX(2,70)", 0.013),
                ("microsatelliteLongerThanX(2,75)", 0.014),
                ("microsatelliteLongerThanX(2,80)", 0.015),
                ("microsatelliteLongerThanX(2,85)", 0.016),
                ("microsatelliteLongerThanX(2,90)", 0.017),
                ("microsatelliteLongerThanX(2,95)", 0.018),
                ("microsatelliteLongerThanX(2,100)", 0.019),
            ],
        )

    @staticmethod
    def trinucleotid_runs(data):
        """
        Calculates the dropchance based on the chance of trinucleotid microsatellites to mutate. Microsatellites are
        repeats of two to six units that have a higher chance to mutate than regular polymers.
        The dropchances are additive in this case (@shouldDrop()).
        :param data: The DNA sequence to check for microsatellites with length 3.
        :return: The dropchance based on the occurence of microsatellites with length 3.
        """
        return shouldDrop(
            data,
            [
                ("microsatelliteLongerThanX(3,10)", 0.001),
                ("microsatelliteLongerThanX(3,15)", 0.002),
                ("microsatelliteLongerThanX(3,20)", 0.003),
                ("microsatelliteLongerThanX(3,25)", 0.004),
                ("microsatelliteLongerThanX(3,30)", 0.005),
                ("microsatelliteLongerThanX(3,35)", 0.006),
                ("microsatelliteLongerThanX(3,40)", 0.007),
                ("microsatelliteLongerThanX(3,45)", 0.008),
                ("microsatelliteLongerThanX(3,50)", 0.009),
                ("microsatelliteLongerThanX(3,55)", 0.010),
                ("microsatelliteLongerThanX(3,60)", 0.011),
                ("microsatelliteLongerThanX(3,65)", 0.012),
                ("microsatelliteLongerThanX(3,70)", 0.013),
                ("microsatelliteLongerThanX(3,75)", 0.014),
                ("microsatelliteLongerThanX(3,80)", 0.015),
                ("microsatelliteLongerThanX(3,85)", 0.016),
                ("microsatelliteLongerThanX(3,90)", 0.017),
                ("microsatelliteLongerThanX(3,95)", 0.018),
                ("microsatelliteLongerThanX(3,100)", 0.019),
            ],
        )

    @staticmethod
    def long_strands(data):
        """
        Calculates the dropchance based on the length of the sequence since longer strands have a higher chance to mutate.
        :param data: The sequence to check for the strand length.
        :return: The dropchance based on the strand length.
        """
        return shouldDropMax(
            data,
            [
                ("len(117)", 0.001),
                ("len(130)", 0.005),
                ("len(150)", 0.01),
                ("len(170)", 0.02),
                ("len(200)", 0.025),
                ("len(250)", 0.05),
                ("len(300)", 0.1),
                ("len(350)", 0.2),
            ],
        )

    @staticmethod
    def random_permutations(data):
        """
        Calculates the dropchance for simulated random mutations in the DNA-Data.
        :param data: The sequence to simulate random mutations for.
        :return: The dropchance based on random mutations.
        """
        return shouldDrop(data, [("*", 0.02)])

    @staticmethod
    def illegal_symbols(data):
        """
        Checks the DNA data for illegal symbols and returns a dropchance of 1.0 if the sequence contains them.
        :param data: The sequence to check for illegal symbols.
        :return: The dropchance based on the occurence of illegal symbols (0.0 or 1.0).
        """
        res = shouldDrop(data, [("strContainsIllegalChars(ACGT)", 1.0)])
        return res

    """ Roughly based on http://www.reinhardheckel.com/papers/DNA_channel_statistics.pdf """

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
        preamble = "charCountBiggerEqualThanX("
        return shouldDrop(
            data,
            [
                (preamble + ch + ",20)", 0.002),
                (preamble + ch + ",40)", 0.002),
                (preamble + ch + ",60)", 0.002),
                (preamble + ch + ",80)", 0.002),
                (preamble + ch + ",100)", 0.002),
                (preamble + ch + ",120)", 0.002),
                (preamble + ch + ",140)", 0.002),
                (preamble + ch + ",160)", 0.002),
                (preamble + ch + ",180)", 0.002),
                (preamble + ch + ",200)", 0.002),
            ],
        )

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
        preamble = "charCountBiggerEqualThanX("
        return shouldDrop(
            data,
            [
                (preamble + ch + ",20)", 0.001),
                (preamble + ch + ",40)", 0.001),
                (preamble + ch + ",60)", 0.001),
                (preamble + ch + ",80)", 0.001),
                (preamble + ch + ",100)", 0.001),
                (preamble + ch + ",120)", 0.001),
                (preamble + ch + ",140)", 0.001),
                (preamble + ch + ",160)", 0.001),
                (preamble + ch + ",180)", 0.001),
                (preamble + ch + ",200)", 0.001),
            ],
        )

    @staticmethod
    def a_permutation(data):
        """
        Calls @ch_at_permutation for 'A'
        :param data: Sequence
        :return: Dropchance
        """
        return DNARules.ch_at_permutation(data, "A")

    @staticmethod
    def t_permutation(data):
        """
        Calls @ch_at_permutation for 'T'
        :param data: Sequence
        :return: Dropchance
        """
        return DNARules.ch_at_permutation(data, "T")

    @staticmethod
    def g_permutation(data):
        """
        Calls @ch_cg_permutation for 'G'
        :param data: Sequence
        :return: Dropchance
        """
        return DNARules.ch_cg_permutation(data, "G")

    @staticmethod
    def c_permutation(data):
        """
        Calls @ch_cg_permutation for 'C'
        :param data: Sequence
        :return: Dropchance
        """
        return DNARules.ch_cg_permutation(data, "C")

    @staticmethod
    def overall_gc_content(data):
        """
        Calculates the dropchance of the sequence based on the GC-content since GC-rich DNA is more stable.
        :param data: The sequence to check for the GC-content.
        :return: The Dropchance based on the GC-content.
        """
        gc_con = gc_content(data)
        dropchance = (100 + (175 * gc_con) / 6 - (121 * gc_con ** 2) / 72 + (gc_con ** 3) / 36 - (
                gc_con ** 4) / 7200) / 100
        return max(0.0, min(1.0, dropchance))

    @staticmethod
    def windowed_gc_content(data, window_size=50):
        """

        :param data:
        :param window_size:
        :return:
        """
        chunks = [DNARules.overall_gc_content(data[i:i + window_size]) for i in range(0, len(data), window_size)]
        res = max(chunks)
        return min(1.0, res)

    @staticmethod
    def motif_search(data):
        """

        :param data:
        :return:
        """
        return shouldDrop(
            data,
            [
                # Promoter recognition motif (Euk).
                ("strContainsSub(TATAAA)", 0.01),
                # Promoter recognition motifs (Prok).
                ("strContainsSub(TTGACA)", 0.05),
                ("strContainsSub(TGTATAATG)", 0.05),
                # Polyadenylation signals (Euk).
                ("strContainsSub(AATAAA)", 0.01),
                ("strContainsSub(TTGTGTGTTG)", 0.01),
                # Lox sites.
                ("strContainsSub(ATAACTTCGTATAGCATACATTATACGAAGTTAT)", 1.01),
                ("strContainsSub(ATAACTTCGTATAGCATACATTATACGAACGGTA)", 1.01),
                ("strContainsSub(TACCGTTCGTATAGCATACATTATACGAAGTTAT)", 1.01),
                ("strContainsSub(TACCGTTCGTATAGCATACATTATACGAACGGTA)", 1.01),
                ("strContainsSub(TACCGTTCGTATATGGTATTATATACGAAGTTAT)", 1.01),
                ("strContainsSub(TACCGTTCGTATATTCTATCTTATACGAAGTTAT)", 1.01),
                ("strContainsSub(TACCGTTCGTATAGGATACTTTATACGAAGTTAT)", 1.01),
                ("strContainsSub(TACCGTTCGTATATACTATACTATACGAAGTTAT)", 1.01),
                ("strContainsSub(TACCGTTCGTATACTATAGCCTATACGAAGTTAT)", 1.01),
                ("strContainsSub(ATAACTTCGTATATGGTATTATATACGAACGGTA)", 1.01),
                ("strContainsSub(ATAACTTCGTATAGTATACCTTATACGAAGTTAT)", 1.01),
                # Lox site spacers not covered by the Lox sites.
                ("strContainsSub(AGGTATGC)", 1.01),
                ("strContainsSub(TTGTATGG)", 1.01),
                ("strContainsSub(GGATAGTA)", 1.01),
                ("strContainsSub(GTGTATTT)", 1.01),
                ("strContainsSub(GGTTACGG)", 1.01),
                ("strContainsSub(TTTTAGGT)", 1.01),
                ("strContainsSub(GTACACAT)", 1.01),
                # Restriction enzyme recognition motifs.
                # BpiI
                ("strContainsSub(GAAGAC)", 1.01),
                # inverse BpiI
                ("strContainsSub(CTTCTG)", 1.01),
                # BsaI
                ("strContainsSub(GGTCTC)", 1.01),
                # inverse BsaI
                ("strContainsSub(CCAGAG)", 1.01),

                ("strContainsSub(CGTCTC)", 0.01),
                ("strContainsSub(GCGATG)", 0.01),
                ("strContainsSub(GCTCTTC)", 0.01),
                # Oligo Adapters.
                ("strContainsSub(CTCGTAGACTGCGTACCA)", 0.01),
                ("strContainsSub(GACGATGAGTCCTGAGTA)", 0.01),
                # 5' extensions.
                ("strContainsSub(GGTTCCACGTAAGCTTCC)", 0.01),
                ("strContainsSub(GCGATTACCCTGTACACC)", 0.01),
                ("strContainsSub(GCCAGTACATCAATTGCC)", 0.01),
                # Twister Adapters:
                ("strContainsSub(GAAGTGCCATTCCGCCTGACCT)", 1.0),  # Twister 5' Adapter
                ("strContainsSub(AGGCTAGGTGGAGGCTCAGTG)", 1.0),  # Twister 3' Adapter

            ]
        )

    @staticmethod
    def motif_regex_search(data):
        """

        :param data:
        :return:
        """
        return shouldDrop(
            data,
            [
                # Promoter recognition motif (Euk).
                ("strContainsSubRegex(CANYYY)", 0.01),
                ("strContainsSubRegex(ANCCAATCA)", 0.01),
                ("strContainsSubRegex(KGGGCGGRRY)", 0.01),
                ("strContainsSubRegex(KRGGCGKRRY)", 0.01),
                # Promoter recognition motifs (Prok).
                ("strContainsSubRegex(AAAWWTWTTTTNNNAAA)", 0.05),
                # Ribosomal binding site (Euk).
                ("strContainsSubRegex(RCCACCATGG)", 0.05),
                # Ribosomal binding site (Prok).
                ("strContainsSubRegex(AGGAGGACAGCTAUG)", 0.05),
                # Lox sites.
                ("strContainsSubRegex(ATAACTTCGTATAGTAYACATTATACGAAGTTAT)", 0.01),
            ]
        )

    @staticmethod
    def add_complementary(input_lst):
        def revert(in_chr):
            return {"G": "C", "C": "G", "T": "A", "A": "T"}[in_chr]

        tmp = copy.deepcopy(input_lst)

        for x in input_lst:
            tmp.append("".join([revert(chr) for chr in x]))
        return tmp

    @staticmethod
    def simple_motif_search(data):
        return 1.0 if any([x in data for x in DNARules.add_complementary(["ATAACTTCGTATAGCATACATTATACGAAGTTAT",
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

    @staticmethod
    def apply_all_rules(packet):
        """

        :param packet:
        :return:
        """
        return DNARules.apply_all_rules_with_data(packet)[0]

    @staticmethod
    def apply_all_rules_with_data(packet):
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
            for x in [
                DNARules.a_permutation,
                DNARules.t_permutation,
                DNARules.c_permutation,
                DNARules.g_permutation,
                DNARules.dinucleotid_runs,
                DNARules.homopolymers,
                DNARules.overall_gc_content,
                #  DNARules.long_strands, DNARules.illegal_symbols,
                DNARules.trinucleotid_runs,
                DNARules.random_permutations,
                DNARules.motif_search,
                DNARules.motif_regex_search,
                DNARules.windowed_gc_content,
                DNARules.repeatRegion,
                DNARules.smallRepeatRegion
            ]
        ]
        return sum(res_arr), res_arr, packet


if __name__ == "__main__":
    print(DNARules.smallRepeatRegion("GCGCGCAATA"))
    """
    a = []
    for file in os.listdir("/home/michael/Code/norec4dna/adad_LT_Dorn"):
        with open("/home/michael/Code/norec4dna/adad_LT_Dorn/" + file, 'r') as f:
            content = f.read()
            # print(file)
            a.append(DNARules.simple_motif_search(content))
    print(sum(a))
    """
