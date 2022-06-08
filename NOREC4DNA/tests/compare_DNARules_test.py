import pytest

from norec4dna.rules.DNARules import DNARules
from norec4dna.rules.FastDNARules import FastDNARules


@pytest.mark.parametrize("params", [("AA", 0.001), ("AAAAA", 0.05), ("AAAAAAAA", 1.0)])
def test_homopolymers(params):
    assert fast_comp(DNARules.homopolymers(params[0]), 3) == fast_comp(FastDNARules.homopolymers(params[0])) == params[
        1]


@pytest.mark.parametrize("params", [("AAAAATTTTT", 1.0), ("GGGGGCCCCC", 1.0), ("ATGC", 0.00), ("AAAAAAGCGC", 0.0),
                                    ("GCGCGCCATA", 1.0)])
def test_gc_content(params):
    assert round(DNARules.overall_gc_content(params[0]), 2) == round(FastDNARules.overall_gc_content(params[0]), 2) == \
           params[1]


@pytest.mark.parametrize("params", [("ACCAAAAAAAAAAAAATTAGACCAAAAAAAAAAAAATTAG", 1.0),
                                    ("AAAGTAGATAGATAGAAACACACACAGTACACACA", 0.0), ("AAA", 0.0)])
def test_repeat_region(params):
    assert fast_comp(DNARules.repeatRegion(params[0])) == fast_comp(FastDNARules.repeatRegion(params[0])) == params[1]


@pytest.mark.parametrize("params", [("ACGACGTTT", 1.00), ("ACGTGCATACACGTGCATACACGTGCATAC", 1.00)])
def test_small_repeat_region(params):
    assert fast_comp(DNARules.smallRepeatRegion(params[0])) == fast_comp(FastDNARules.smallRepeatRegion(params[0])) == \
           params[1]


@pytest.mark.parametrize("params", [("AC" * 10, 0.001), ("AC" * 50, 0.045), ("AC" * 100, 0.19)])
def test_dinucleotid_runs(params):
    assert fast_comp(DNARules.dinucleotid_runs(params[0])) == fast_comp(FastDNARules.dinucleotid_runs(params[0])) == \
           params[1]


@pytest.mark.parametrize("params", [("ACG" * 10, 0.001), ("ACG" * 50, 0.045), ("ACG" * 100, 0.19)])
def test_trinucleotid_runs(params):
    assert fast_comp(DNARules.trinucleotid_runs(params[0])) == fast_comp(FastDNARules.trinucleotid_runs(params[0])) == \
           params[1]


@pytest.mark.parametrize("params", ["A" * 170, "A" * 200, "A" * 350, "A"])
def test_long_strands(params):
    assert abs(DNARules.long_strands(params) - FastDNARules.long_strands(params) < 0.01)


@pytest.mark.parametrize("params", [("ACGACGTACGACTG", 0.02), ("ACGACGTACGACGTAGCAGTGACGTACG", 0.02)])
def test_random_permutations(params):
    assert DNARules.random_permutations(params[0]) == FastDNARules.random_permutations(params[0]) == params[1]


@pytest.mark.parametrize("params", [("ACGACTGACTX", 1.0), ("AACGACTAGCG", 0.0)])
def test_illegal_symbols(params):
    assert DNARules.illegal_symbols(params[0]) == FastDNARules.illegal_symbols(params[0]) == params[1]


@pytest.mark.parametrize("params", [("A" * 20, 0.001), ("A" * 100, 0.005), ("A" * 200, 0.01), ("C" * 20, 0.000)])
def test_a_permutation(params):
    assert fast_comp(DNARules.a_permutation(params[0])) == fast_comp(FastDNARules.a_permutation(params[0])) == params[1]


@pytest.mark.parametrize("params", [("T" * 20, 0.001), ("T" * 100, 0.005), ("T" * 200, 0.01), ("G" * 20, 0.000)])
def test_t_permutation(params):
    assert fast_comp(DNARules.t_permutation(params[0])) == fast_comp(FastDNARules.t_permutation(params[0])) == params[1]


@pytest.mark.parametrize("params", [("G" * 20, 0.002), ("G" * 100, 0.01), ("G" * 200, 0.02), ("T" * 20, 0.000)])
def test_g_permutation(params):
    assert fast_comp(DNARules.g_permutation(params[0])) == fast_comp(FastDNARules.g_permutation(params[0])) == params[1]


@pytest.mark.parametrize("params", [("C" * 20, 0.002), ("C" * 100, 0.01), ("C" * 200, 0.02), ("A" * 20, 0.000)])
def test_c_permutation(params):
    assert fast_comp(DNARules.c_permutation(params[0])) == fast_comp(FastDNARules.c_permutation(params[0])) == params[1]


@pytest.mark.parametrize("params", [
    "ATTAGCGTATCCAATCAGCTGACACCAAAAAAAAAAAATTAGCGTATCCAATCAGCTGACACCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGCAAGCAAAAAAAAAAAAAAAAAAAAAAATTAGCGTATCCAATCAGCTGACACCAAAAAAAAAAAAAAAAAAAAAAAAAAGCAAGCAAAAAACCATGGTCGAAAAAAAAAAAAAAAAAAGCAAGCAAAAAAAAAAAAAAAAAAAAAAAGCAAGCAAAAAACCATGGTCG",
    "AGCCGGCGCTAGATGCCGACCCCACGCACGGCCGCCAGAACTTTCGCCATGAGCCGGCGCTAGATGCCGACGCATTCGCCCGACCGAGATACCGTGCTATAGCGAGTCATGTCGCAAGAACGTGCGCCAGCTCTGGCGCCCAAGCGAACGTGCGCACAGTATGGCGAGCGGGCTCCCGCGCGTGAGTAAGAACTTTCTATATGGATCGAGCACGGCCTTGCTCACGCCAACTTTCGCCATATTCGCCCTCTT"])
def test_windowed_gc_content(params):
    assert abs(DNARules.windowed_gc_content(params) - FastDNARules.windowed_gc_content(params)) == 0


@pytest.mark.parametrize("params", [("TTGACA", 0.05), ("A", 0.00), ("GCTCTTCAATAAA", 0.02)])
def test_motif_search(params):
    assert fast_comp(DNARules.motif_search(params[0])) == fast_comp(FastDNARules.motif_search(params[0])) == params[1]


@pytest.mark.parametrize("params", [("AGGAGGACAGCTAUG", 0.05), ("A", 0.00), ("GCCACCATGGCATCCC", 0.06)])
def test_motif_regex_search(params):
    assert fast_comp(DNARules.motif_regex_search(params[0])) == fast_comp(FastDNARules.motif_regex_search(params[0])) == \
           params[1]


@pytest.mark.parametrize("params", [("A", 0.00), ("GCATACAT", 1.0)])
def test_simple_motif_search(params):
    assert fast_comp(DNARules.simple_motif_search(params[0])) == fast_comp(
        FastDNARules.simple_motif_search(params[0])) == params[1]


def fast_comp(num, dec=3):
    p = float(10 ** dec)
    return int(num * p + 0.5) / p
