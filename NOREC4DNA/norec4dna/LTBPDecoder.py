#!/usr/bin/python
# -*- coding: latin-1 -*-
import os
import struct
import time
import typing
import numpy as np
from collections import deque

from norec4dna.BPDecoder import BPDecoder
from norec4dna.DecodePacket import DecodePacket
from norec4dna.ErrorCorrection import nocode, crc32
from norec4dna.HeaderChunk import HeaderChunk
from norec4dna.Packet import Packet
from norec4dna.distributions.Distribution import Distribution
from norec4dna.helper import calc_crc, xor_mask
from norec4dna.helper.quaternary2Bin import quat_file_to_bin

"""
    |  len(packed)  | total_number_of_chunks  | len(used_packets) |         used_packets        |        Data       |           CRC32           |
    |  I (4 bytes)  |      I (4 byte)      |     I (4 byte)   | len(used_packets)*H (2byte) | len(packed_data)*s | L (unsiged long) (4 byte) |
                    |__________________________________________________________________________________________|  -----------^
            ^------ |______________________________________________________________________________________________________________________|
"""


class LTBPDecoder(BPDecoder):
    def __init__(self, file: typing.Optional[str] = None, error_correction: typing.Callable = nocode,
                 use_headerchunk: bool = True, static_number_of_chunks: typing.Optional[int] = None,
                 implicit_mode: bool = True, dist: typing.Optional[Distribution] = None):
        super().__init__(file, error_correction, use_headerchunk, static_number_of_chunks)
        self.implicit_mode: bool = implicit_mode
        self.use_headerchunk: bool = use_headerchunk
        self.file: str = file
        self.decodedPackets: typing.Dict = {}
        self.degreeToPacket: typing.Dict[int, typing.Set] = {}
        if file is not None:
            self.isFolder: bool = os.path.isdir(file)
            if not self.isFolder:
                self.f = open(self.file, "rb")
        self.correct: int = 0
        self.corrupt: int = 0
        self.number_of_chunks: int = 1000000
        self.headerChunk: typing.Optional[HeaderChunk] = None
        self.queue: deque = deque()
        self.error_correction: typing.Callable = error_correction
        self.static_number_of_chunks: int = static_number_of_chunks
        self.dist: typing.Optional[Distribution] = dist  # if implicit_mode is True, dist MUST be != None

    def decodeFolder(self, packet_len_format: str = "I", crc_len_format: str = "L",
                     number_of_chunks_len_format: str = "I", degree_len_format: str = "I", seed_len_format: str = "I",
                     last_chunk_len_format: str = "I") -> typing.Optional[int]:
        decoded: bool = False
        self.EOF: bool = False
        if self.static_number_of_chunks is not None:
            self.number_of_chunks = self.static_number_of_chunks
            number_of_chunks_len_format = ""  # if we got static number_of_chunks we do not need it in struct string
        if self.implicit_mode:
            degree_len_format = ""
        for file in os.listdir(self.file):
            if file.endswith(".LT") or file.endswith("DNA"):
                self.EOF = False
                if file.endswith("DNA"):
                    self.f = quat_file_to_bin(self.file + "/" + file)
                else:
                    self.f = open(self.file + "/" + file, "rb")
                new_pack = self.getNextValidPacket(True, packet_len_format=packet_len_format,
                                                   crc_len_format=crc_len_format,
                                                   number_of_chunks_len_format=number_of_chunks_len_format,
                                                   degree_len_format=degree_len_format, seed_len_format=seed_len_format,
                                                   last_chunk_len_format=last_chunk_len_format)
                if new_pack is not None:
                    ## koennte durch input_new_packet ersetzt werden:
                    self.addPacket(new_pack)
                    decoded = self.updatePackets(new_pack)
                if decoded:
                    break
                ##
        print("Decoded Packets: " + str(self.correct))
        print("Corrupt Packets : " + str(self.corrupt))
        if hasattr(self, "f"):
            self.f.close()
        if not decoded and self.EOF:
            print("Unable to retrieve File from Chunks. Too much errors?")
            return -1

    def decodeFile(self, packet_len_format: str = "I", crc_len_format: str = "L",
                   number_of_chunks_len_format: str = "I", degree_len_format: str = "I", seed_len_format: str = "I",
                   last_chunk_len_format: str = "I") -> typing.Optional[int]:
        decoded: bool = False
        self.EOF: bool = False
        if self.static_number_of_chunks is not None:
            self.number_of_chunks = self.static_number_of_chunks
            number_of_chunks_len_format = ""  # if we got static number_of_chunks we do not need it in struct string
        if self.implicit_mode:
            degree_len_format = ""
        while not (decoded or self.EOF):
            new_pack: Packet = self.getNextValidPacket(False, packet_len_format=packet_len_format,
                                                       crc_len_format=crc_len_format,
                                                       number_of_chunks_len_format=number_of_chunks_len_format,
                                                       degree_len_format=degree_len_format,
                                                       seed_len_format=seed_len_format,
                                                       last_chunk_len_format=last_chunk_len_format)
            if new_pack is None:
                break
            ## koennte durch input_new_packet ersetzt werden:
            self.addPacket(new_pack)
            decoded = self.updatePackets(new_pack)
            ##
        print("Decoded Packets: " + str(self.correct))
        print("Corrupt Packets : " + str(self.corrupt))
        self.f.close()
        if not decoded and self.EOF:
            print("Unable to retrieve File from Chunks. Too much errors?")
            return -1

    def input_new_packet(self, packet: Packet) -> bool:
        """ Used for easy PseudoDecode """
        if not isinstance(packet, DecodePacket):
            packet: Packet = DecodePacket.from_packet(packet)
        self.addPacket(packet)
        return self.updatePackets(packet)

    def addPacket(self, packet: Packet) -> None:
        if (not packet.get_degree() in self.degreeToPacket) or (
                not isinstance(self.degreeToPacket[packet.get_degree()], set)):
            self.degreeToPacket[packet.get_degree()] = set()
        self.number_of_chunks = packet.get_total_number_of_chunks()
        self.degreeToPacket[packet.get_degree()].add(packet)

    def updatePackets(self, packet: Packet) -> bool:
        self.queue.append(packet)
        finished: bool = False
        while len(self.queue) > 0 and not finished:
            finished = self.reduceAll(self.queue.popleft())
        return finished

    def compareAndReduce(self, packet: Packet, other: Packet) -> typing.Union[bool, int]:
        if self.file is None:  # In case of PseudoDecode: DO NOT REALLY COMPUTE XOR
            packet.remove_packets(other.get_used_packets())
        else:
            packet.xor_and_remove_packet(other)
        degree = packet.get_degree()
        if (degree not in self.degreeToPacket) or (not isinstance(self.degreeToPacket[degree], set)):
            self.degreeToPacket[degree] = set()
        self.degreeToPacket[degree].add(packet)
        if self.is_decoded():
            return True
        self.queue.append(packet)
        return degree

    """def reduceAll(self, packet: Packet) -> bool:
        # looup all packets for this to solve with ( when this packet has a subset of used Packets)
        fin: bool = False

        lookup: typing.List[int] = [i for i in self.degreeToPacket.keys() if packet.get_degree() < i]
        for i in lookup:
            if not isinstance(self.degreeToPacket[i], set):
                self.degreeToPacket[i] = set()
            for p in self.degreeToPacket[i].copy():
                p_used = p.get_used_packets()
                pack_used = packet.get_used_packets()
                if len(pack_used) < len(p_used) and pack_used.issubset(p_used):
                    self.degreeToPacket[i].remove(p)
                    degree = self.compareAndReduce(p, packet)
                    if isinstance(degree, bool) and degree:
                        return degree
        degree = packet.get_degree()
        lookup = [i for i in self.degreeToPacket.keys() if packet.get_degree() > i]
        for i in lookup:
            if not isinstance(self.degreeToPacket[i], set):
                self.degreeToPacket[i] = set()
            for p in self.degreeToPacket[i].copy():
                p_used = p.get_used_packets()
                pack_used = packet.get_used_packets()
                if len(pack_used) > len(p_used) and p_used.issubset(pack_used):
                    try:
                        self.degreeToPacket[degree].remove(packet)
                        degree = self.compareAndReduce(packet, p)
                        if isinstance(degree, bool) and degree:
                            return degree
                    except Exception:
                        continue
        return fin or self.is_decoded()"""

    def is_decoded(self) -> bool:
        return (1 in self.degreeToPacket and len(
            self.degreeToPacket[1]) == self.number_of_chunks)  # self.number_of_chunks

    def getSolvedCount(self) -> int:
        return len(self.degreeToPacket[1])

    def getNextValidPacket(self, from_multiple_files: bool = False, packet_len_format: str = "I",
                           crc_len_format: str = "L", number_of_chunks_len_format: str = "I",
                           degree_len_format: str = "I", seed_len_format: str = "I",
                           last_chunk_len_format: str = "I") -> typing.Optional[Packet]:
        if not from_multiple_files:
            packet_len: typing.Union[bytes, int] = self.f.read(struct.calcsize("<" + packet_len_format))
            packet_len = struct.unpack("<" + packet_len_format, packet_len)[0]
            packet: bytes = self.f.read(int(packet_len))
        else:
            packet = self.f.read()
            packet_len = len(packet)
        if not packet or not packet_len:  # EOF
            self.EOF: bool = True
            self.f.close()
            return None
        crc_len: typing.Optional[int] = -struct.calcsize("<" + crc_len_format)
        if self.error_correction.__code__.co_name == crc32.__code__.co_name:
            payload = packet[:crc_len]
            # instead of typing.Any we would have _SupportsIndex:
            crc: typing.Union[int, typing.Any] = struct.unpack("<" + crc_len_format, packet[crc_len:])[0]
            calced_crc = calc_crc(payload)
            if crc != calced_crc:  # If the Packet is corrupt, try next one
                print("[-] CRC-Error - " + str(hex(crc)) + " != " + str(hex(calced_crc)))
                self.corrupt += 1
                return self.getNextValidPacket(from_multiple_files, packet_len_format=packet_len_format,
                                               crc_len_format=crc_len_format,
                                               number_of_chunks_len_format=number_of_chunks_len_format,
                                               degree_len_format=degree_len_format, seed_len_format=seed_len_format,
                                               last_chunk_len_format=last_chunk_len_format)
        else:
            crc_len = None
            try:
                packet = self.error_correction(packet)
            except:
                self.corrupt += 1
                return self.getNextValidPacket(from_multiple_files)

        struct_str = "<" + number_of_chunks_len_format + degree_len_format + seed_len_format
        struct_len = struct.calcsize(struct_str)
        len_data = struct.unpack(struct_str, packet[0:struct_len])
        degree = None
        if self.static_number_of_chunks is None:
            if self.implicit_mode:
                number_of_chunks, seed = len_data
            else:
                number_of_chunks, degree, seed = len_data
            self.number_of_chunks = xor_mask(number_of_chunks, number_of_chunks_len_format)
        else:
            if self.implicit_mode:
                seed, = len_data
            else:
                degree, seed = len_data
        seed = xor_mask(seed, seed_len_format)
        if degree is None:
            self.dist.set_seed(seed)
            degree = self.dist.getNumber()
        else:
            degree = xor_mask(degree, degree_len_format)
        used_packets = self.choose_packet_numbers(degree, seed=seed)
        data = packet[struct_len:crc_len]

        self.correct += 1
        res = DecodePacket(data, used_packets, error_correction=self.error_correction,
                           number_of_chunks=self.number_of_chunks)
        if used_packets.issubset({0}) and self.headerChunk is None and self.use_headerchunk:
            self.headerChunk = HeaderChunk(res)
        return res

    def choose_packet_numbers(self, degree: int, seed: int = 0) -> typing.Set:
        assert degree <= self.number_of_chunks
        res = set()
        rng = np.random
        rng.seed(seed)
        for _ in range(0, degree):
            tmp = rng.choice(range(0, self.number_of_chunks))
            while tmp in res:
                tmp = rng.choice(range(0, self.number_of_chunks))
            res.add(tmp)
        return res

    def saveDecodedFile(self, last_chunk_len_format: str = "I", null_is_terminator: bool = False,
                        print_to_output: bool = True) -> None:
        assert self.is_decoded(), "Can not save File: Unable to reconstruct."
        file_name: str = "DEC_" + self.file.split("\x00")[0]  # split is needed for weird  MAC / Windows bugs...
        sort_list: typing.List = sorted(self.degreeToPacket[1])
        if 0 in sort_list[0].get_used_packets() and self.use_headerchunk:
            self.headerChunk = HeaderChunk(sort_list[0])
        output_concat: bytes = b""
        if self.headerChunk is not None:
            file_name = self.headerChunk.get_file_name().decode("utf-8")
        with open(file_name, "wb") as f:
            for decoded in sort_list:
                if 0 in decoded.get_used_packets() and self.use_headerchunk:
                    self.headerChunk = HeaderChunk(decoded, last_chunk_len_format=last_chunk_len_format)
                else:
                    if self.number_of_chunks - 1 in decoded.get_used_packets() and self.use_headerchunk:
                        output = decoded.get_data()[0: self.headerChunk.get_last_chunk_length()]
                        if type(output) == bytes:
                            output_concat += output
                        else:
                            output_concat += output.tobytes()
                        f.write(output)
                    else:
                        if null_is_terminator:
                            data = decoded.get_data()
                            if type(data) == bytes:
                                splitter = data.decode().split("\x00")
                            else:
                                splitter = data.tostring().decode().split("\x00")
                            output = splitter[0].encode()
                            output_concat += output
                            f.write(output)
                            if len(splitter) > 1:
                                break  # since we are in null-terminator mode, we exit once we see the first 0-byte
                        else:
                            output = decoded.get_data()
                            if type(output) == np.ndarray or type(output) != bytes:
                                output_concat += output.tobytes()
                            else:
                                output_concat += output
                            f.write(output)

        print("Saved file as '" + str(file_name) + "'")
        if print_to_output:
            print("Result:")
            print(output_concat.decode("utf-8"))


if __name__ == "__main__":
    filename = "LT_logo.jpg"
    p_start = time.time()
    x = LTBPDecoder(filename)
    x.decode()
    p_end = time.time() - p_start
    print(
        "LT_Approx_Decode,"
        + str(x.number_of_chunks)
        + ","
        + str(len(x.decodedPackets))
        + ","
        + str(p_end)
        + ","
        + filename
    )
    x.saveDecodedFile()
