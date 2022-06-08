#!/usr/bin/python
# -*- coding: latin-1 -*-
import struct
import typing
import numpy as np

from norec4dna.helper import xor_mask
from norec4dna.Packet import Packet
from bitstring import BitArray
from norec4dna.helper.RU10Helper import intermediate_symbols
from norec4dna.distributions.RaptorDistribution import RaptorDistribution
from norec4dna.ErrorCorrection import nocode


class RU10Packet(Packet):
    def __init__(self, data, used_packets: typing.Collection[int], total_number_of_chunks, id, dist=None,
                 read_only=False,
                 error_correction=nocode, packet_len_format="I", crc_len_format="L", number_of_chunks_len_format="L",
                 id_len_format="L", save_number_of_chunks_in_packet=True, method=None, window=None, prepend="",
                 append=""):
        self.id: int = id
        self.bool_arrayused_packets: typing.Optional[typing.List[bool]] = None
        self.total_number_of_chunks: int = total_number_of_chunks
        self.data: bytes = data
        self.used_packets: typing.Optional[typing.Iterable[int]] = None
        self.internal_hash: typing.Optional[int] = None
        self.set_used_packets(used_packets)
        self.degree: int = len(used_packets)
        self.dna_data: typing.Optional[str] = None
        self.packet_len_format: str = packet_len_format
        self.crc_len_format: str = crc_len_format
        self.number_of_chunks_len_format: str = number_of_chunks_len_format
        self.id_len_format: str = id_len_format
        self.error_correction: typing.Callable[[typing.Any], typing.Any] = error_correction
        self.save_number_of_chunks_in_packet: bool = save_number_of_chunks_in_packet
        self.error_prob: typing.Optional[float] = None
        if method:
            self.method: typing.Optional[str] = method
            self.window: typing.Optional[int] = window
            self.packedMethod: typing.Optional[bytes] = self.packMethod()
        else:
            self.packedMethod = None
        if not read_only and (len(self.data) > 0 or self.data != ""):
            self.packed_used_packets = self.prepare_and_pack()
            self.packed = self.calculate_packed_data()
        else:
            self.packed_used_packets = None
            self.packed = None
        if dist is None:
            self.dist = RaptorDistribution(total_number_of_chunks)
        else:
            self.dist = dist
        _, self.s, self.h = intermediate_symbols(total_number_of_chunks, self.dist)
        self.prepend = prepend
        self.append = append
        # super().__init__(data, used_packets, total_number_of_chunks, read_only, error_correction=error_correction)

    def set_used_packets(self, u_packets):
        self.used_packets = u_packets
        self.internal_hash = hash(frozenset(i for i in self.used_packets))
        tmp_lst = np.full((1, self.total_number_of_chunks), False)
        for x in self.used_packets:
            if x < self.total_number_of_chunks:
                tmp_lst[0, x] = True
        self.bool_arrayused_packets = tmp_lst[0]
        self.update_degree()
        # [
        #    x in self.used_packets for x in range(0, self.total_number_of_chunks)
        # ]

    def prepare_and_pack(self) -> bytes:
        # Format = Highest possible Packetnumber for this file,
        # number of used Packets for this File and the seed for the Indices of the used Packets
        if self.save_number_of_chunks_in_packet:
            return struct.pack("<" + self.number_of_chunks_len_format + self.id_len_format,
                               xor_mask(self.total_number_of_chunks, self.number_of_chunks_len_format),
                               xor_mask(self.id, self.id_len_format))
        else:
            return struct.pack("<" + self.id_len_format, xor_mask(self.id, self.id_len_format))

    def packMethod(self) -> bytes:
        if "window" not in self.method:
            if self.method == "even":
                data = BitArray(bin='00101010').tobytes()
            elif self.method == "odd":
                data = BitArray(bin='01101010').tobytes()
            else:
                raise RuntimeError("Unknown method: ", self.method)
        else:
            if self.method == "window_30":
                data_str = '10'
            elif self.method == "window_40":
                data_str = '01'
            else:
                raise RuntimeError("Unknown method: ", self.method)
            bin_win = bin(self.window)[2:]
            while len(data_str) + len(bin_win) < 8:
                data_str += '0'
            data_str += bin_win
            data = BitArray(bin=data_str).tobytes()
        return data

    def calculate_packed_data(self) -> bytes:
        # size of the packets + UsedPackets + Data + crc
        self.packed_data = struct.pack("<" + str(len(self.data)) + "s", bytes(self.data))
        if self.packedMethod:
            payload = struct.pack(
                "<" + str(len(self.packed_used_packets)) + "s" + str(len(self.packed_data)) + "s" + str(
                    len(self.packedMethod)) + "s",  # method data
                self.packed_used_packets, self.packed_data, self.packedMethod)
        else:
            payload = struct.pack("<" + str(len(self.packed_used_packets)) + "s" + str(len(self.packed_data)) + "s",
                                  self.packed_used_packets, self.packed_data)
        return self.error_correction(payload)  # proxy payload through dynamic error correction / detection

    def getId(self) -> int:
        return self.id

    def setId(self, id: int):
        self.id = id

    @classmethod
    def from_packet(cls, packet: 'RU10Packet', pseudo: bool = False) -> 'RU10Packet':
        if not pseudo:
            data = packet.get_data()
        else:
            data = ""
        used_packets = packet.get_used_packets()
        number_of_packets = packet.get_total_number_of_chunks()
        res = cls(data, used_packets, number_of_packets, packet.getId())
        res.error_correction = packet.get_error_correction()
        return res

    def get_number_of_half_blocks(self) -> int:
        return self.h

    def get_number_of_ldpc_blocks(self) -> int:
        return self.s

    def get_bool_array_used_packets(self) -> typing.Optional[typing.List[bool]]:
        return self.bool_arrayused_packets

    def get_bool_array_all_used_packets(self) -> typing.List[bool]:
        return [x in self.used_packets for x in
                range(self.total_number_of_chunks + self.get_number_of_ldpc_blocks() + self.get_number_of_half_blocks())]

    def get_bool_array_used_and_ldpc_packets(self) -> typing.List[bool]:
        # speedup candidate
        u_bound = self.total_number_of_chunks + self.get_number_of_ldpc_blocks()
        tmp_lst = np.full((1, u_bound), False)
        for x in self.used_packets:
            if x < u_bound:
                tmp_lst[0, x] = True
        res = tmp_lst[0]
        del tmp_lst
        return res

    def get_bool_array_ldpc_packets(self) -> typing.List[bool]:
        # speedup candidate
        return [x in self.used_packets for x in
                range(self.total_number_of_chunks, self.total_number_of_chunks + self.get_number_of_ldpc_blocks(), )]

    def get_bool_array_half_packets(self) -> typing.List[bool]:
        return [x in self.used_packets for x in
                range(self.total_number_of_chunks + self.get_number_of_ldpc_blocks(),
                      self.total_number_of_chunks + self.get_number_of_ldpc_blocks() + self.get_number_of_half_blocks(), )]

    def get_bool_array_repair_packets(self) -> typing.List[bool]:
        return [x in self.used_packets for x in
                range(self.total_number_of_chunks,
                      self.total_number_of_chunks + self.get_number_of_ldpc_blocks() + self.get_number_of_half_blocks(), )]

    def __str__(self) -> str:
        return "< used_packets: " + str(self.used_packets) + " , Data: " + str(self.data) + " >"


if __name__ == "__main__":
    print("This class must not be called by itself.")
