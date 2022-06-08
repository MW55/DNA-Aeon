#!/usr/bin/python
# -*- coding: latin-1 -*-
import os, typing

try:
    from cdnarules import byte2QUATS as b2quats
    from cdnarules import getQUAT
except ImportError:
    print("C Module failed to load, falling back to slow mode")

"""
Mapping:
A = 0
C = 1
G = 2
T = 3
"""


def getQUAT(bit1: bool, bit2: bool):
    if not (bit1 or bit2):
        return "A"
    elif (not bit1) and bit2:
        return "C"
    elif bit1 and (not bit2):
        return "G"
    elif bit1 and bit2:
        return "T"
    else:
        print("ERROR, this might never happen.")
        return "E"


def old_b2quats(byte):
    res = ""
    if not isinstance(byte, int):
        byte = byte[0]
    if not isinstance(byte, str):
        byt = iter(bin(byte)[2:].rjust(8, "0"))
    else:
        byt = iter(bin(ord(byte))[2:].rjust(8, "0"))
    for x, y in zip(byt, byt):
        res += getQUAT(str2bool(x), str2bool(y))
    return res


def byte2QUATS(byte: typing.Union[str, bytes]):
    try:
        return b2quats(byte)
    except:
        return old_b2quats(byte)


def bin2Quaternary(filename: str):
    with open(filename, "rb") as f:
        with open(filename + ".quat", "w") as ff:
            byte = f.read(1)
            while byte != b"":
                ff.write(byte2QUATS(byte))
                byte = f.read(1)


def string2QUATS(text: typing.Union[str, bytes, typing.List[typing.AnyStr]]):
    return [byte2QUATS(x) for x in text]


def str2bool(s: str):
    return s == "1"


def quads2dna(quads: typing.Union[typing.List[int], bytes]):
    translation = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
    return "".join(translation[x] for x in quads)


def main():
    folder = "RU10_b_lq.webm"
    for filename in os.listdir(folder):
        if filename.endswith(".RU10"):
            bin2Quaternary(folder + "/" + filename)


if __name__ == "__main__":
    # main()
    filename = (
        "vergleich bzgl ACGT-verteilung/mit dna rules/RU10_logo.jpg/0.RU10"
    )  # ""RU10_b_lq.webm/1.RU10"  # "raptor.pdf"
    x = [0, 0, 0, 1, 2, 3, 1, 1, 1, ]
    print(quads2dna(x))
    print(b2quats(x))


    def bitstring_to_bytes(s: typing.Union[str, bytes, bytearray]):
        v = int(s, 2)
        b = bytearray()
        while v:
            b.append(v & 0xff)
            v >>= 8
        return bytes(b[::-1])


    s = "001011000110100101011110001010110000110110101110011111000101000100000011111110110001011010010011101110001111101010110100"
    print("".join(string2QUATS(bitstring_to_bytes(s))))
