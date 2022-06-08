"""
Mapping:
A = 0
C = 1
G = 2
T = 3
"""
import typing
from io import BytesIO


def quaternary_to_bin(filename: str) -> None:
    with open(filename, "r") as f:
        with open("OUT_" + filename.split(".quat")[0], "wb") as ff:
            quats = f.read(4)
            while quats != "":
                ff.write(quats_to_bytes(quats))
                quats = f.read(4)


def quat_file_to_bin(file: str) -> BytesIO:
    # we could try to recover from various errors: e.g. appended newline: strip newline and retry
    # or we might want to replace all characters =! [A,C,G,T] with "" and retry reading
    # (optimally writing to a temporary file to not polute our working dir)
    out = BytesIO()
    with open(file, "r") as f:
        quats = f.read(4)
        while quats != "":
            out.write(quats_to_bytes(quats))
            quats = f.read(4)
    out.seek(0)
    return out


def quad_file_to_bytes(filename: str) -> BytesIO:
    with open(filename, "r") as f:
        seq = f.read()
        return BytesIO(dna2quads(seq))


def quats_to_bytes(quats: typing.AnyStr) -> bytearray:
    byte = bytearray(
        [
            (get_quarter_byte(quats[0]) << 6)
            + (get_quarter_byte(quats[1]) << 4)
            + (get_quarter_byte(quats[2]) << 2)
            + (get_quarter_byte(quats[3]))
        ]
    )
    return byte


def tranlate_quat_to_byte(in_txt):
    out = b''
    for i in range(0, len(in_txt), 4):
        out += quats_to_bytes(in_txt[i:i + 4])
    return out


def get_quarter_byte(quat: str) -> int:
    if quat == "A":
        return 0b00
    elif quat == "C":
        return 0b01
    elif quat == "G":
        return 0b10
    elif quat == "T":
        return 0b11
    else:
        print("ERROR, this should never happen. Does your inputfile contain characters other than A,C,G,T?")


def dna2quads(dna: typing.AnyStr) -> bytes:
    translation = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    return bytes([translation[x] for x in dna])


if __name__ == "__main__":
    filename = "raptor.pdf.quat"
    quaternary_to_bin(filename)
