# Stores methods for fallback if c-extensions are not installed.
import warnings
import numpy
import typing


def bitSet(x: int, b: int) -> bool:
    return ((x >> b) & 1) == 1


def bitsSet(x: numpy.uint64) -> int:  # x is of type uint64 !
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        x -= (x >> numpy.uint64(1)) & numpy.uint64(0x5555555555555555)
        x = (x & numpy.uint64(0x3333333333333333)) + (
                (x >> numpy.uint64(2)) & numpy.uint64(0x3333333333333333)
        )
        x = (x + (x >> numpy.uint64(4))) & numpy.uint64(0x0F0F0F0F0F0F0F0F)
        res = numpy.int64((x * numpy.uint64(0x0101010101010101)) >> numpy.uint64(56))
        return res


def grayCode(x: int):
    return numpy.bitwise_xor((numpy.uint64(x) >> numpy.uint64(1)), numpy.uint64(x))


def buildGraySequence(length: int, b: int):
    s = numpy.empty(length, dtype=int)
    i = 0
    x = 0  # numpy.uint64(0)
    while True:
        g = grayCode(x)
        if bitsSet(g) == b:
            s[i] = g
            i += 1
            if i >= length:
                break
        x += 1
    return s


def xor_intern(n_p1, n_p2):
    return numpy.bitwise_xor(n_p1, n_p2)


def r_region(data: str, repeat_length: int = 20) -> float:
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


def small_r_region(data: str, repeat_length: int = 9) -> float:
    count = 1
    for i in range(len(data) - repeat_length):
        subseq = data[i:i + repeat_length]
        if data[i + 1:].find(subseq) >= 0:
            count += 1
    return 1.0 if 1.0 * count * repeat_length / len(data) > 0.44 else count * repeat_length / len(data) * 0.5


def microsatellite_python(text: typing.AnyStr, lengthToLookFor: int) -> typing.Tuple[int, str]:
    i = 0
    n = len(text)
    res = 1
    res_chars = text[:lengthToLookFor]
    max_lenght = 0
    while i <= n - 2 * lengthToLookFor:
        if text[i: i + lengthToLookFor] == text[i + lengthToLookFor: i + 2 * lengthToLookFor]:
            res += 1  # we found one
        else:
            if max_lenght < res:
                max_lenght = res
                res_chars = text[i: i + lengthToLookFor]
            res = 1
        i += lengthToLookFor
    if max_lenght < res:
        max_lenght = res
        res_chars = text[i: i + lengthToLookFor]
    return max_lenght, res_chars


def longestSequenceOfChar_python(text: typing.AnyStr, char_x="*") -> typing.Tuple[str, int]:
    n = len(text)
    c = 0
    res = char_x
    curr = 1
    i = 0
    while i < n:
        if i < n - 1 and text[i] == text[i + 1]:
            curr += 1
        else:
            if curr > c and (text[i] == char_x or char_x == "*"):
                c = curr
                res = text[i]
            curr = 1
        i += 1
    return res, c


def strContainsSub_python(text: typing.AnyStr, sequence: typing.AnyStr) -> bool:
    res = sequence in text
    return res
