#!/usr/bin/python
# -*- coding: latin-1 -*-
import typing

from norec4dna.Packet import Packet

"""
    Aux Packets only have Data (used_packets are inferred by pseudo-random Number)
    -> they have same length as the real Datachunks
"""


class OnlineAuxPacket(Packet):
    def __init__(self, data: bytes, used_packets: typing.Optional[typing.Set[int]] = None,
                 aux_number: typing.Optional[int] = None, total_number_of_chunks: int = 0):
        super().__init__(data, used_packets, total_number_of_chunks)
        self.data: bytes = data
        self.total_number_of_chunks: int = total_number_of_chunks
        self.used_packets: typing.Optional[typing.Set[int]] = used_packets
        self.update_degree()
        self.error_correction: typing.Callable[[typing.Any], typing.Any] = lambda x: x  # "AUX-Packet - NO CRC"
        self.aux_number: int = aux_number
        self.dna_data: typing.Optional[str] = None
        self.error_prob: typing.Optional[int] = None

    def get_struct(self, split_to_multiple_files: bool = False) -> bytes:  # split_to_multiple_files not used here.
        return self.get_data()

    def get_data(self) -> bytes:
        return self.data

    def set_data(self, data: bytes):
        self.data = data

    def __str__(self):
        return ("< aux_number : " + str(self.aux_number) + ", used_packets: " + str(self.used_packets)
                + " , Data: " + str(self.data) + " , Error Correction: "
                + str(self.error_correction) + " >")

    def __eq__(self, other) -> bool:
        return hash(self) == hash(other)

    def __hash__(self) -> int:
        return hash(frozenset(hash(i) for i in self.used_packets))
