#!/usr/bin/python
# -*- coding: latin-1 -*-
import typing

from norec4dna.ErrorCorrection import nocode
from norec4dna.Packet import Packet
from norec4dna.helper import xor_numpy


class DecodePacket(Packet):
    def __init__(self, data, used_packets, error_correction=nocode, number_of_chunks: int = -1):
        self.set_used_packets(used_packets)
        self.data = data
        self.error_correction = error_correction
        self.did_change = False
        self.internal_hash = None
        self.update_degree()
        self.total_number_of_chunks = number_of_chunks
        self.error_prob = None

    @classmethod
    def from_packet(cls, packet: Packet, pseudo: bool = False):
        if not pseudo:
            data = packet.get_data()
        else:
            data = ""
        used_packets = packet.get_used_packets()
        error_correction = packet.get_error_correction()
        res = cls(data, used_packets, number_of_chunks=packet.total_number_of_chunks,
                  error_correction=error_correction)
        res.total_number_of_chunks = packet.get_total_number_of_chunks()
        return res

    def get_used_packets(self):
        return self.used_packets

    def get_data(self):
        return self.data

    def remove_packets(self, packet_set: typing.Collection):
        assert self.used_packets is not None
        self.set_used_packets(self.used_packets.difference(packet_set))
        self.update_degree()
        self.did_change = True  # CRC is no longer valid

    def xor_and_remove_packet(self, packet):
        self.remove_packets(packet.get_used_packets())
        self.data = xor_numpy(self.data, packet.get_data())
        self.did_change = True  # CRC is no longer valid
