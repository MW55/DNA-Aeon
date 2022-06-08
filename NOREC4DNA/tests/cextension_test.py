import numpy as np
from norec4dna.GEPP import GEPP


def test_gepp():
    a_arr = [[False, True, False, True], [False, False, True, False], [False, True, True, False],
             [False, True, False, False], [True, True, True, False]]
    b_arr = [np.frombuffer(b"Hallo", dtype="uint8"),
             np.frombuffer(b"Welt ", dtype="uint8"),
             np.frombuffer(b"Test3", dtype="uint8"),
             np.frombuffer(b"Test1", dtype="uint8"),
             np.frombuffer(b"12345", dtype="uint8")]
    A = np.array([a_arr[0]], dtype=bool)
    b = np.array([b_arr[0]], dtype="uint8")

    gauss_elim_piv = GEPP(np.copy(A), np.copy(b))
    gauss_elim_piv_c = GEPP(np.copy(A), np.copy(b))
    for a_elem, b_elem in zip(a_arr[1:], b_arr[1:]):
        gauss_elim_piv.addRow(np.copy(a_elem), np.copy(b_elem))
        gauss_elim_piv_c.addRow(np.copy(a_elem), np.copy(b_elem))

    gauss_elim_piv.insert_tmp()
    gauss_elim_piv_c.insert_tmp()

    res_c = gauss_elim_piv_c._elimination()
    res = gauss_elim_piv._py_elimination()

    assert res == res_c
    assert all([(a == b).all() for a, b in zip([gauss_elim_piv.b[i][0] for i in gauss_elim_piv.result_mapping],
                                               [gauss_elim_piv_c.b[i][0] for i in gauss_elim_piv_c.result_mapping])])
    # for all we care, these tests could even fail without breaking the code:
    assert (gauss_elim_piv.result_mapping == gauss_elim_piv_c.result_mapping).all()
    assert (gauss_elim_piv.A == gauss_elim_piv_c.A).all()
    assert (gauss_elim_piv.b[:gauss_elim_piv.m] == gauss_elim_piv_c.b[:gauss_elim_piv.m]).all()


def test_bit_set():
    from norec4dna.helper.helper_cpu_single_core import bitSet
    from norec4dna.helper.fallback_code import bitSet as bitSet_fallback
    for _ in range(100):
        x = np.random.randint(0, 10000000)
        for b in range(16):
            assert bitSet(x, b) == bitSet_fallback(x, b)


def test_bits_set():
    from norec4dna.helper.helper_cpu_single_core import bitsSet
    from norec4dna.helper.fallback_code import bitsSet as bitsSet_fallback
    for _ in range(50):
        x = np.uint64(np.random.randint(0, 100000))
        assert bitsSet(x) == bitsSet_fallback(x)


def test_gray_code():
    from norec4dna.helper.helper_cpu_single_core import grayCode
    from norec4dna.helper.fallback_code import grayCode as grayCode_fallback
    for _ in range(40):
        x = np.uint64(np.random.randint(0, 1000000000))
        assert (grayCode(x) == grayCode_fallback(x))


def test_build_gray_sequence():
    x = np.random.randint(1000, 1500)
    b = np.random.randint(8, 16)
    from norec4dna.helper.helper_cpu_single_core import buildGraySequence
    from norec4dna.helper.fallback_code import buildGraySequence as buildGraySequence_fallback
    assert (buildGraySequence(x, b) == buildGraySequence_fallback(x, b)).all()


def test_r_region():
    from cdnarules import repeatRegion as rRegion
    from norec4dna.helper.fallback_code import r_region as fallback_rRegion
    repeat_length = 20
    for i in range(21):
        data = "A" * i + "ACCCGATCTGCCGCCCCCATCG" + "G" * (10 + i) + "TACACCCGTAGCTGCGCGATGCTAGCTCCGCCC"
        assert rRegion(data, repeat_length) == fallback_rRegion(data, repeat_length)


def test_small_r_region():
    from cdnarules import smallRepeatRegion as smallrRegion
    from norec4dna.helper.fallback_code import small_r_region as fallback_smallrRegion
    repeat_length = 9
    for i in range(21):
        data = "A" * i + "ACCCGATCTGCCGCCCCCATCG" + "G" * (10 + i) + "TACACCCGTAGCTGCGCGATGCTAGCTCCGCCC"
        assert np.round(smallrRegion(data, repeat_length), 5) == np.round(fallback_smallrRegion(data, repeat_length), 5)
    txt = "ACGACGTTT"
    assert np.round(smallrRegion(txt, 9), 5) == np.round(fallback_smallrRegion(txt, 9), 5)


def test_microsatellite():
    from cdnarules import microsatellite
    from norec4dna.helper.fallback_code import microsatellite_python
    length_to_look_for = np.random.randint(2, 5)
    for i in range(21):
        text = "A" * i + "ACCCGATCTGCCGCCCCCATCG" + "G" * (10 + i) + "TACACCCGTAGCTGCGCGATGCTAGCTCCGCCC"
        assert microsatellite(text, length_to_look_for) == microsatellite_python(text, length_to_look_for)
    text = "TACACCCGTAGCTGCGCGATGCTAGCTCCGCCCCCCCCCCCCCCCCCCCCC"
    assert microsatellite(text, length_to_look_for) == microsatellite_python(text, length_to_look_for)


def test_longest_sequence_of_char():
    from cdnarules import longestSequenceOfChar
    from norec4dna.helper.fallback_code import longestSequenceOfChar_python
    for i in range(21):
        for char in ["A", "C", "G", "T", "*"]:
            text = "A" * i + "ACCCGATCTGCCGCCCCCATCG" + "G" * (10 + i) + "TACACCCGTAGCTGCGCGATGCTAGCTCCGCCC"
            tmp_1, tmp_2 = longestSequenceOfChar(text, char)
            assert longestSequenceOfChar_python(text, char) == (tmp_1.decode(), tmp_2)


def test_str_contains_sub():
    from cdnarules import strContainsSub
    from norec4dna.helper.fallback_code import strContainsSub_python
    for i in range(21):
        for char in ["AAAAAAAAA", "ACGGTGC", "GAGAGAGAGAGAGAGAAT", "TTTTTTTT", "GGGGGGGGGGGGGGGGGGGGGG"]:
            text = "A" * i + "ACCCGATCTGCCGCCCCCATCG" + "G" * (10 + i) + "TACACCCGTAGCTGCGCGATGCTAGCTCCGCCC"
            assert strContainsSub_python(text, char) == strContainsSub(text, char)
