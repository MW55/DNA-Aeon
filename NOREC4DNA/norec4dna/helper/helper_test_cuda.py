#!/usr/bin/python
# -*- coding: latin-1 -*-
import os
import numba
import zlib, numpy
from random import random
from functools import reduce
from numba import jit, vectorize, cuda

os.environ["NUMBAPRO_NVVM"] = r"D:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v9.1\nvvm\bin\nvvm64_32_0.dll"
os.environ["NUMBAPRO_LIBDEVICE"] = r"D:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v9.1\nvvm\libdevice"


@jit(nogil=True)
def xor_numpy(p1, p2=None):
    if p2 is not None:
        return xor_numpy(numpy.row_stack((p1, p2)))
    if isinstance(p1, list):
        r = []
        for x in p1:
            if isinstance(x, numpy.ndarray):
                r.append(x.copy(order="C"))
            else:
                if isinstance(x, list):
                    r.append(numpy.asarray(x, dtype="uint8"))
                else:
                    r.append(numpy.fromstring(x, dtype="uint8"))

        p1 = r
    if len(p1) == 1:
        return numpy.asarray(p1)
    if (isinstance(p1, numpy.ndarray)) and (
            (p1.dtype == numpy.uint8) or (p1.dtype == numpy.bool)
    ):
        n_p1 = numpy.asarray(p1)

    else:
        if isinstance(p1, str):
            n_p1 = numpy.fromstring(p1, dtype="uint8")
        else:
            n_p1 = numpy.asarray(p1)
    n_p1 = numpy.concatenate(n_p1).reshape((len(n_p1), -1))
    n_p1 = numpy.ascontiguousarray(n_p1)
    """
    threadsperblock = 32
    blockspergrid = (n_p1.size + (threadsperblock - 1)) // threadsperblock

    sm_size = n_p1.size * n_p1.dtype.itemsize
    """
    threadsperblock = (32, 32)
    blockspergrid_x = numpy.math.ceil(n_p1.shape[0] / threadsperblock[0])
    blockspergrid_y = numpy.math.ceil(n_p1.shape[1] / threadsperblock[1])
    blockspergrid = (blockspergrid_x, blockspergrid_y)
    stream = 0
    sm_size = n_p1.size * n_p1.dtype.itemsize
    arrXOR[blockspergrid, threadsperblock, stream, sm_size](n_p1)
    n_p1[n_p1.shape[0] - 1, :].tolist()

    for x in range(n_p1.shape[0]):
        print(str(x) + " - " + str(n_p1[x, :]))
    return n_p1[0, :]


@cuda.jit
def arrXOR(arr):
    sm = cuda.shared.array(shape=0, dtype=numba.uint8)

    arrshape = arr.shape[0]
    x = cuda.blockIdx.x * cuda.blockDim.x + cuda.threadIdx.x
    y = cuda.blockIdx.y * cuda.blockDim.y + cuda.threadIdx.y
    sm[x * cuda.blockDim.x + y] = arr[x, y]
    """
    # Thread id in a 1D block
    tx = cuda.threadIdx.x
    # Block id in a 1D grid
    ty = cuda.blockIdx.x
    # Block width, i.e. number of threads per block
    bw = cuda.blockDim.x
    # Compute flattened index inside the array
    j = tx + ty * bw
    """
    if y < arr.shape[1]:
        if (not arr.shape[0] % 2 == 0) and arr.shape[0] > 2:
            # remove last element..
            sm[(arr.shape[0] - 2) * cuda.blockDim.x + y] = (
                    sm[(arr.shape[0] - 2) * cuda.blockDim.x + y]
                    ^ sm[(arr.shape[0] - 1) * cuda.blockDim.x + y]
            )
            i = (arr.shape[0] - 1) // 2
            arrshape = arr.shape[0] - 1
        else:
            i = arr.shape[0] // 2
        i = i - x

        while i > 1:
            if i * 2 + x < arrshape:
                sm[(2 * i + x) * cuda.blockDim.x + y] = x

            i = i >> 1
            i = i - x
        arr[x, y] = sm[x * cuda.blockDim.x + y]


@jit(cache=True)
def listXOR(plist):
    if len(plist) == 1:
        return numpy.fromstring(plist[0], dtype=numpy.uint8)

    if isinstance(plist[0], numpy.ndarray):
        tmp = []
        for x in plist:
            if isinstance(x, numpy.ndarray):
                tmp.append(x.tolist())
            elif isinstance(x, str) or isinstance(x, bytes):
                tmp.append(numpy.fromstring(x, dtype="uint8").tolist())
            else:
                tmp.append([x])
        plist = tmp
    return xor_numpy(plist)


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
@vectorize(["uint64(uint64)", "uint64(uint32)", "uint64(float64)"], target="parallel")
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
