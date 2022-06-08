#!/usr/bin/python
# -*- coding: latin-1 -*-
from math import ceil
import struct
import numpy as np
import typing

from norec4dna.ErrorCorrection import nocode
from norec4dna.distributions.OnlineDistribution import OnlineDistribution
from norec4dna.Packet import Packet
from norec4dna.helper import xor_mask


class OnlinePacket(Packet):
    def __init__(self, data: bytes, total_number_of_chunks: int, quality: int, epsilon: float, check_block_number: int,
                 used_packets: typing.Optional[typing.Union[typing.List[int], typing.Set[int]]] = None,
                 dist: typing.Optional[OnlineDistribution] = None, read_only: bool = False,
                 error_correction: typing.Callable[[typing.Any], typing.Any] = nocode, crc_len_format: str = "L",
                 number_of_chunks_len_format: str = "I", quality_len_format: str = "I", epsilon_len_format: str = "f",
                 check_block_number_len_format: str = "I", save_number_of_chunks_in_packet: bool = True, prepend="",
                 append=""):
        self.data: bytes = data
        self.total_number_of_chunks: int = total_number_of_chunks
        self.quality: int = quality
        self.epsilon: float = epsilon
        self.error_correction: typing.Callable[[typing.Any], typing.Any] = error_correction
        self.check_block_number: int = check_block_number
        self.id: int = self.check_block_number
        self.dna_data: typing.Optional[str] = None
        self.used_packets: typing.Optional[typing.Set[int]] = None
        self.number_of_chunks_len_format: str = number_of_chunks_len_format
        self.quality_len_format: str = quality_len_format
        self.epsilon_len_format: str = epsilon_len_format
        self.check_block_number_len_format: str = check_block_number_len_format
        self.id_len_format: str = self.check_block_number_len_format
        self.crc_len_format: str = crc_len_format
        self.save_number_of_chunks_in_packet: bool = save_number_of_chunks_in_packet
        self.internal_hash: typing.Optional[int] = None
        self.error_prob: typing.Optional[int] = None
        self.packet_len_format: typing.Optional[str] = None
        self.prepend = prepend  # not tested for Online-Code
        self.append = append  # not tested for Online-Code
        if not read_only:
            self.packedInfo: bytes = self.prepare_and_pack()
            self.packed: bytes = self.calculate_packed_data()
        self.aux_number: int = -1
        if dist is not None:
            self.dist: OnlineDistribution = dist
        else:
            self.dist = OnlineDistribution(self.epsilon)
        if used_packets is None:
            self.calcUsedPackets()
        else:
            self.set_used_packets(used_packets)
            self.update_degree()

    def getQuality(self) -> int:
        return self.quality

    def getEpsilon(self) -> float:
        return self.epsilon

    def prepare_and_pack(self) -> bytes:
        # Format = Highest possible Packetnumber for this file, quality settings, epsilon, and a (hopefully) unique checkBlock-Number
        struct_format = "<" + (
            self.number_of_chunks_len_format if self.save_number_of_chunks_in_packet else "") + self.quality_len_format + self.epsilon_len_format + self.check_block_number_len_format
        if self.save_number_of_chunks_in_packet:
            return struct.pack(struct_format, xor_mask(self.total_number_of_chunks, self.number_of_chunks_len_format),
                               xor_mask(self.quality, self.quality_len_format), self.epsilon, self.check_block_number)
        else:
            return struct.pack(struct_format, xor_mask(self.quality, self.quality_len_format), self.epsilon,
                               self.check_block_number)

    def calculate_packed_data(self) -> bytes:
        # Länge des Packets + Infos + Data + crc
        payload = struct.pack("<" + str(len(self.packedInfo)) + "s" + str(len(self.data)) + "s", self.packedInfo,
                              bytes(self.data))
        return self.error_correction(payload)

    def get_data(self) -> bytes:
        return self.data

    def set_used_packets(self, used_packets: typing.Set[int]):
        self.used_packets = used_packets
        self.internal_hash = hash(frozenset(i for i in self.used_packets))
        self.update_degree()

    def calcUsedPackets(self):
        self.dist.set_seed(self.check_block_number)
        degree: int = self.dist.getNumber()
        self.set_used_packets(self.choose_packet_numbers(degree, self.check_block_number))

    def choose_packet_numbers(self, degree: int, seed: int) -> typing.Set[int]:
        rng = np.random
        rng.seed(seed)
        res: typing.Set[int] = set()
        for _ in range(0, degree):
            tmp = rng.choice(range(0, self.total_number_of_chunks + self.getNumberOfAuxBlocks()))  # +1 for HeaderChunk
            while tmp in res:
                tmp = rng.choice(
                    range(0, self.total_number_of_chunks + self.getNumberOfAuxBlocks()))  # +1 for HeaderChunk
            res.add(tmp)
        return res

    def getNumberOfAuxBlocks(self) -> int:
        return ceil(0.55 * self.quality * self.epsilon * self.total_number_of_chunks)

    def get_bool_array_used_packets(self) -> typing.List[bool]:
        return [x in self.used_packets for x in range(self.total_number_of_chunks)]

    def getBoolArrayAuxPackets(self) -> typing.List[bool]:
        return [x in self.used_packets
                for x in range(self.total_number_of_chunks, self.total_number_of_chunks + self.getNumberOfAuxBlocks(), )
                ]

    def __str__(self) -> str:
        return ("< used_packets: " + str(self.used_packets) + " , Data: " + str(
            self.data) + " , Error Correction: " + str(self.error_correction) + " >")

    def __hash__(self) -> int:
        return hash(frozenset(hash(i) for i in self.used_packets))
