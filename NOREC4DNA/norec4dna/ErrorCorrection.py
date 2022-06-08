import struct
import typing

from norec4dna.helper import calc_crc
from reedsolo import RSCodec
from bitarray import bitarray

rscodec: typing.Optional[RSCodec] = None
number_r_symbols: int = 2


def nocode(txt: typing.AnyStr) -> typing.AnyStr:
    return txt


def crc32(txt: typing.Union[bytes, str, bytearray], crc_len_format="L") -> bytes:
    crc = calc_crc(txt)
    packed = struct.pack("<" + str(len(txt)) + "s" + crc_len_format, txt, crc)
    return packed


def crc32_decode(txt: typing.Union[bytes, str, bytearray], crc_len_format="L") -> bytes:
    crc_len = -struct.calcsize("<" + crc_len_format)
    crc = struct.unpack("<" + crc_len_format, txt[crc_len:])[0]
    payload = txt[:crc_len]
    calced_crc = calc_crc(payload)
    assert crc == calced_crc, "CRC-Error - " + str(hex(crc)) + " != " + str(hex(calced_crc))
    return payload


def reed_solomon_encode(txt: typing.Union[bytes, str, bytearray], number_repair_symbols: int = 2) -> bytes:
    global rscodec, number_r_symbols
    if rscodec is None or number_r_symbols != number_repair_symbols:
        rscodec = RSCodec(number_repair_symbols)
        number_r_symbols = number_repair_symbols
    return bytes(rscodec.encode(txt))


def reed_solomon_decode(txt: typing.Union[bytes, str, bytearray], number_repair_symbols: int = 2) -> bytes:
    global rscodec, number_r_symbols
    if rscodec is None or number_r_symbols != number_repair_symbols:
        rscodec = RSCodec(number_repair_symbols)
        number_r_symbols = number_repair_symbols
    return bytes(rscodec.decode(txt))


def dna_reed_solomon_encode(txt: typing.Union[bytes, str, bytearray], number_repair_symbols: int = 2, c_exp: int = 2,
                            prim_poly: int = 7) -> bytes:
    """
    Warning: dna_reed_solomon_* requires a custom RSCodec version with support for custom c_exp and nsize.
    """
    tmp_rscodec = RSCodec(number_repair_symbols, c_exp=c_exp, prim=prim_poly, nsize=2 ** c_exp - 1)
    dec_list = bits_to_dec(txt)
    return bytes(tmp_rscodec.encode(dec_list))


def dna_reed_solomon_decode(txt: typing.Union[bytes, str, bytearray], number_repair_symbols: int = 2, c_exp: int = 2,
                            prim_poly: int = 7) -> bytes:
    """
    Warning: dna_reed_solomon_* requires a custom RSCodec version with support for custom c_exp and nsize.
    """
    tmp_rscodec = RSCodec(number_repair_symbols, c_exp=c_exp, prim=prim_poly, nsize=2 ** c_exp - 1)
    decoded: bytearray = tmp_rscodec.decode(txt)
    return bytes(bits_to_bytes(decoded))


# Helper functions.
def bits_to_bytes(decoded_bytes: typing.Union[bytearray, typing.List[int], bytes]) -> bytes:
    decoded_bits = dec_to_bits(decoded_bytes)
    return bitarray(decoded_bits).tobytes()


def str_to_bits(input_string: str) -> str:
    return ''.join(format(x, "08b") for x in input_string)


def bytes_to_bits(input_bytes: bytes) -> str:
    bits = bitarray(endian='big')
    bits.frombytes(input_bytes)
    return bits.to01()


def bits_to_dec(input_string: bytes) -> typing.List[int]:
    translation = {'00': 0, '01': 1, '10': 2, '11': 3}
    input_bits = bytes_to_bits(input_string)
    if len(input_bits) % 2 == 1:
        print(input_string)
        print(input_bits)
    return [translation[input_bits[i:i + 2]] for i in range(0, len(input_bits), 2)]


def dec_to_bits(decoded_bytes: typing.List[int]) -> str:
    translation = {0: '00', 1: '01', 2: '10', 3: '11'}
    return "".join([translation[bits] for bits in decoded_bytes])


def get_error_correction_decode(e_correction: str, repair_symbols: int):
    if e_correction == "nocode":
        error_correction = nocode
    elif e_correction == "crc":
        error_correction = crc32_decode
    elif e_correction == "reedsolomon":
        if repair_symbols != 2:
            error_correction = lambda x: reed_solomon_decode(x, repair_symbols)
        else:
            error_correction = reed_solomon_decode
    elif e_correction == "dna_reedsolomon":
        if repair_symbols != 2:
            error_correction = lambda x: dna_reed_solomon_decode(x, repair_symbols)
        else:
            error_correction = dna_reed_solomon_decode
    else:
        raise NotImplementedError(
            "Selected Error Correction not supported, choose: 'nocode', 'crc', 'reedsolomon' or 'dna_reedsolomon'")
    return error_correction


def get_error_correction_encode(e_correction: str, repair_symbols: int):
    if e_correction == "nocode":
        error_correction = nocode
    elif e_correction == "crc":
        error_correction = crc32
    elif e_correction == "reedsolomon":
        if repair_symbols != 2:
            error_correction = lambda x: reed_solomon_encode(x, repair_symbols)
        else:
            error_correction = reed_solomon_encode
    elif e_correction == "dna_reedsolomon":
        if repair_symbols != 2:
            error_correction = lambda x: dna_reed_solomon_encode(x, repair_symbols)
        else:
            error_correction = dna_reed_solomon_encode
    else:
        raise NotImplementedError(
            "Selected Error Correction not supported, choose: 'nocode', 'crc', 'reedsolomon' or 'dna_reedsolomon'")
    return error_correction
