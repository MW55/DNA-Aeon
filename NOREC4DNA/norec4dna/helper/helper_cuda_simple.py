#!/usr/bin/python
# -*- coding: latin-1 -*-
from random import random

from numba import jit, vectorize
import zlib, numpy
from functools import reduce
import os

os.environ[
    "NUMBAPRO_NVVM"
] = r"D:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v9.1\nvvm\bin\nvvm64_32_0.dll"
os.environ[
    "NUMBAPRO_LIBDEVICE"
] = r"D:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v9.1\nvvm\libdevice"


@jit(nogil=True)
def xor_numpy(p1, p2):
    if (isinstance(p2, numpy.ndarray) and isinstance(p1, numpy.ndarray)) and (
        (p1.dtype == numpy.uint8 and p2.dtype == numpy.uint8)
        or (p1.dtype == numpy.bool and p2.dtype == numpy.bool)
    ):
        n_p1 = p1
        n_p2 = p2
    else:
        n_p1 = numpy.fromstring(p1, dtype="uint8")
        n_p2 = numpy.fromstring(p2, dtype="uint8")
    return xor_intern(n_p1, n_p2)


@vectorize(["uint8(uint8,uint8)", "boolean(boolean,boolean)"], target="cuda")
def xor_intern(n_p1, n_p2):
    return n_p1 ^ n_p2


@jit(cache=True)
def listXOR(plist):
    return reduce(xor_numpy, plist)


@jit(cache=True)
def logical_xor(plist):
    return reduce(numpy.logical_xor, plist)


@jit(cache=True)
def xor_pakets(packet1, packet2):
    assert len(packet1) == len(packet2)
    a = [a ^ b for (a, b) in zip(bytes(packet1, "utf-8"), bytes(packet2, "utf-8"))]
    return a


@jit(cache=True)
def should_drop_packet(rules, packet):
    rand = random()  # create number from [0,1)
    drop_chance = rules.apply_all_rules(packet)
    packet.set_error_prob(drop_chance)
    # drop packet if rand bigger than the drop_chance for this Packet.
    return drop_chance > rand


@jit(cache=True)
def calc_crc(data):
    return zlib.crc32(data) & 0xFFFFFFFF


@jit(cache=True)
def xor_mask(data, mask=0b11110000101011000011101001101010):
    return numpy.bitwise_xor(data, mask)


# bitSet returns true if x has the b'th bit set
@vectorize(["boolean(uint32, uint32)", "boolean(uint64, uint64)"], target="cuda")
def bitSet(x, b):
    return ((x >> b) & 1) == 1


# bitsSet returns how many bits in x are set.
# This algorithm basically uses shifts and ANDs to sum up the bits in
# a tree fashion.
@vectorize(["int64(uint64)", "int64(uint32)"], target="cuda")
def bitsSet(x):  # x is of type uint64 !
    x -= (x >> numpy.uint64(1)) & numpy.uint64(0x5555555555555555)
    x = (x & numpy.uint64(0x3333333333333333)) + (
        (x >> numpy.uint64(2)) & numpy.uint64(0x3333333333333333)
    )
    x = (x + (x >> numpy.uint64(4))) & numpy.uint64(0x0F0F0F0F0F0F0F0F)
    res = numpy.int64((x * numpy.uint64(0x0101010101010101)) >> numpy.uint64(56))
    return res


# grayCode calculates the gray code representation of the input argument
# The Gray code is a binary representation in which successive values differ
# by exactly one bit. See http://en.wikipedia.org/wiki/Gray_code
@vectorize(["uint64(uint64)", "uint64(uint32)", "uint64(float64)"], target="cuda")
def grayCode(x):
    return (numpy.uint64(x) >> numpy.uint64(1)) ^ numpy.uint64(x)


# buildGraySequence returns a sequence (in ascending order) of "length" Gray numbers,
# all of which have exactly "b" bits set.
@jit("int32(int32,int32)")
def buildGraySequence(length, b):
    s = numpy.empty(length, dtype=int)
    i = 0
    x = numpy.uint64(0)
    while True:
        g = grayCode(x)
        if bitsSet(g) == b:
            s[i] = int(g)
            i += 1
            if i >= length:
                break
        x += 1
    return s
