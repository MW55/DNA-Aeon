#!/usr/bin/python
# -*- coding: latin-1 -*-
from random import random
import zlib, numpy
from functools import reduce, lru_cache
import typing

try:
    from cdnarules import xorArray as xor_intern
except ImportError:
    from norec4dna.helper.fallback_code import xor_intern


def xor_numpy(p1, p2):
    if (isinstance(p2, numpy.ndarray) and isinstance(p1, numpy.ndarray)) and (
            (p1.dtype == numpy.uint8 and p2.dtype == numpy.uint8)
            or (p1.dtype == numpy.int64 and p2.dtype == numpy.int64)
            or (p1.dtype == numpy.bool and p2.dtype == numpy.bool)
    ):
        n_p1 = p1
        n_p2 = p2
    else:
        n_p1 = numpy.frombuffer(p1, dtype="uint8")
        n_p2 = numpy.frombuffer(p2, dtype="uint8")
    return xor_intern(n_p1, n_p2)


def listXOR(plist):
    return reduce(xor_numpy, plist)


def logical_xor(plist):
    return reduce(numpy.logical_xor, plist)


def xor_pakets(packet1, packet2):
    assert len(packet1) == len(packet2)
    a = [a ^ b for (a, b) in zip(bytes(packet1, "utf-8"), bytes(packet2, "utf-8"))]
    return a


def should_drop_packet(rules, packet: 'Packet', upper_bound: float = 1.0, limit_only: bool = True) -> bool:
    rand = upper_bound * random()  # create number from [0, upper_bound)
    drop_chance = rules.apply_all_rules(packet)
    if type(drop_chance) == list:
        drop_chance = drop_chance[0]
    packet.set_error_prob(drop_chance)
    # print(str(rand) + " , " + str(drop_chance))
    # drop packet if rand bigger than the drop_chance for this Packet.
    return (drop_chance > upper_bound) if limit_only else (drop_chance > rand)


def calc_crc(data) -> int:
    return zlib.crc32(data) & 0xFFFFFFFF


@lru_cache(maxsize=1024)
def xor_mask(data: typing.Union[int, float, bytes, numpy.ndarray, typing.Iterable], len_format: str = "I",
             mask: int = 0b11111001110000110110111110011100):
    if len_format == "B":
        return data
    if len_format == "H":
        mask = 0b1111100111000011
    if len_format == "I":
        mask = 0b11111001110000110110111110011100
    if len_format == "Q":
        mask = 0b1111100111000011011011111001110011111001110000110110111110011100
    with numpy.errstate(over="ignore"):
        return numpy.bitwise_xor(data, mask)


try:
    from cdnarules import bitSet as bitSet_c


    def bitSet(x: int, b: int) -> bool:
        return bitSet_c(int(x), int(b))
except ImportError:
    print("BitSet - C Module failed to load, falling back to slow mode")
    from norec4dna.helper.fallback_code import bitSet

try:
    from cdnarules import bitsSet as bitsSet_c


    def bitsSet(x: numpy.uint64) -> int:
        return bitsSet_c(int(x))
except ImportError:
    print("BitsSet - C Module failed to load, falling back to slow mode")
    from norec4dna.helper.fallback_code import bitsSet

try:
    from cdnarules import grayCode as grayCode_c


    def grayCode(x: int):
        return grayCode_c(int(x))
except ImportError:
    print("Gray-Code - C Module failed to load, falling back to slow mode")

    from norec4dna.helper.fallback_code import grayCode

try:
    from cdnarules import buildGraySequence
except ImportError:
    print("Graysequence - C Module failed to load, falling back to slow mode")
    from norec4dna.helper.fallback_code import buildGraySequence
