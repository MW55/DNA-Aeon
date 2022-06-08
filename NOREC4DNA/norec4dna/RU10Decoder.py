#!/usr/bin/python
# -*- coding: latin-1 -*-
import argparse
import io
import os
import struct
import typing
import numpy as np
from math import floor, ceil
from PIL import Image
from io import BytesIO

from norec4dna.helper.RU10Helper import intermediate_symbols, from_true_false_list, choose_packet_numbers
from norec4dna.distributions.Distribution import Distribution
from norec4dna.ErrorCorrection import get_error_correction_decode, nocode
from norec4dna.Packet import Packet
from norec4dna.HeaderChunk import HeaderChunk
from norec4dna.RU10IntermediatePacket import RU10IntermediatePacket
from norec4dna.RU10Packet import RU10Packet
from norec4dna.helper import logical_xor, xor_mask, buildGraySequence, bitSet
from norec4dna.distributions.RaptorDistribution import RaptorDistribution
from norec4dna.GEPP import GEPP
from norec4dna.Decoder import Decoder
from norec4dna.helper.quaternary2Bin import quat_file_to_bin, quad_file_to_bytes, tranlate_quat_to_byte
from zipfile import ZipFile

DEBUG = False


class RU10Decoder(Decoder):
    def __init__(self, file: typing.Optional[str] = None, error_correction=nocode, use_headerchunk: bool = True,
                 static_number_of_chunks: typing.Optional[int] = None, use_method: bool = False):
        self.debug = False
        super().__init__()
        self.isPseudo: bool = False
        self.file: typing.Optional[str] = file
        self.degreeToPacket: dict = {}
        self.use_method: bool = use_method
        if file is not None:
            self.isFolder = os.path.isdir(file)
            self.isZip = file.endswith(".zip")
            if not self.isFolder:
                self.f = open(self.file, "rb")
        self.correct: int = 0
        self.corrupt: int = 0
        self.number_of_chunks: int = 1000000
        self.headerChunk: typing.Optional[HeaderChunk] = None
        self.GEPP: typing.Optional[GEPP] = None
        self.pseudoCount: int = 0
        self.ldpcANDhalf: typing.Dict[int, RU10IntermediatePacket] = dict()
        self.repairBlockNumbers: dict = dict()
        self.s: int = -1
        self.h: int = -1
        self.distribution: typing.Optional[Distribution] = None
        self.EOF: bool = False
        self.counter: dict = dict()
        self.count: bool = True
        self.error_correction: typing.Callable = error_correction
        self.use_headerchunk: bool = use_headerchunk
        self.static_number_of_chunks: typing.Optional[int] = static_number_of_chunks

    def decodeZip(self, packet_len_format: str = "I", crc_len_format: str = "I",
                  number_of_chunks_len_format: str = "I", id_len_format: str = "I"):
        if hasattr(self, "f"):
            self.f.close()
        decoded = False
        self.EOF = False
        if self.static_number_of_chunks is not None:
            self.number_of_chunks = self.static_number_of_chunks
            number_of_chunks_len_format = ""  # if we got static number_of_chunks we do not need it in struct string
        archive = ZipFile(self.file, 'r')
        namelist = archive.namelist()
        try:
            nam = [x.split("_") for x in namelist]
            sorted_by_second = sorted(nam, key=lambda tup: float(tup[1]), reverse=False)
            namelist = [x[0] + "_" + x[1] for x in sorted_by_second]
        except Exception:
            pass
        for name in namelist:
            self.f = io.BytesIO(archive.read(name))
            new_pack = self.getNextValidPacket(True, packet_len_format=packet_len_format,
                                               crc_len_format=crc_len_format,
                                               number_of_chunks_len_format=number_of_chunks_len_format,
                                               id_len_format=id_len_format)
            if hasattr(self, "f"):
                self.f.close()
            if new_pack is None:
                break
            # koennte durch input_new_packet ersetzt werden:
            # self.addPacket(new_pack)
            if new_pack != "CORRUPT":
                decoded = self.input_new_packet(new_pack)
            else:
                if DEBUG: print(f"Packet with name={name} corrupt.")
            if decoded:
                break
            if self.progress_bar is not None:
                self.progress_bar.update(self.correct, Corrupt=self.corrupt)
            ##
        print("Decoded Packets: " + str(self.correct))
        print("Corrupt Packets: " + str(self.corrupt))

        if self.GEPP is None:
            print("No Packet was correctly decoded. Check your configuration.")
            return -1
        if self.GEPP.isPotentionallySolvable():
            decoded = self.GEPP.solve()
        if not decoded and self.EOF:
            print("Unable to retrieve File from Chunks. Too many errors?")
            return -1

    def decodeFolder(self, packet_len_format: str = "I", crc_len_format: str = "I",
                     number_of_chunks_len_format: str = "I", id_len_format: str = "I"):
        """
        Decodes the information from a folder if self.file represents a folder and the packets were saved
        in multiple files and prints the number of decoded and corrupted packets.
        :param packet_len_format: Format of the packet length
        :param crc_len_format:  Format of the crc length
        :param number_of_chunks_len_format: Format of the number of chunks length
        :param id_len_format: Format of the ID length
        :return: -1 if the decoding wasn't successful
        """
        decoded = False
        self.EOF = False
        if self.static_number_of_chunks is not None:
            self.number_of_chunks = self.static_number_of_chunks
            number_of_chunks_len_format = ""  # if we got static number_of_chunks we do not need it in struct string
        for file_in_folder in os.listdir(self.file):
            if file_in_folder.endswith(".RU10") or file_in_folder.endswith("DNA"):
                self.EOF = False
                if file_in_folder.endswith("DNA"):
                    if self.error_correction.__name__ == 'dna_reed_solomon_decode':
                        try:
                            self.f = quad_file_to_bytes(self.file + "/" + file_in_folder)
                        except TypeError:
                            print("skipping CORRUPT file - contains illegal character(s)")
                            self.corrupt += 1
                            continue
                    else:
                        try:
                            self.f = quat_file_to_bin(self.file + "/" + file_in_folder)
                        except TypeError:
                            print("skipping CORRUPT file - contains illegal character(s)")
                            self.corrupt += 1
                            continue
                else:
                    self.f = open(self.file + "/" + file_in_folder, "rb")
                new_pack = self.getNextValidPacket(True, packet_len_format=packet_len_format,
                                                   crc_len_format=crc_len_format,
                                                   number_of_chunks_len_format=number_of_chunks_len_format,
                                                   id_len_format=id_len_format)
                if new_pack is not None:
                    # koennte durch input_new_packet ersetzt werden:
                    # self.addPacket(new_pack)
                    decoded = self.input_new_packet(new_pack)
                if decoded:
                    break
            if self.progress_bar is not None:
                self.progress_bar.update(self.correct, Corrupt=self.corrupt)
            ##
        print("Decoded Packets: " + str(self.correct))
        print("Corrupt Packets: " + str(self.corrupt))
        if hasattr(self, "f"):
            self.f.close()
        if self.GEPP is None:
            print("No Packet was correctly decoded. Check your configuration.")
            return -1
        if self.GEPP.isPotentionallySolvable():
            decoded = self.GEPP.solve()
        if not decoded and self.EOF:
            print("Unable to retrieve File from Chunks. Too many errors?")
            return -1

    def decodeFile(self, packet_len_format: str = "I", crc_len_format: str = "L",
                   number_of_chunks_len_format: str = "I", id_len_format: str = "I"):
        """
        Decodes the information from a file if self.file represents a file and the packets were saved in a single file.
        :param packet_len_format: Format of the packet length
        :param crc_len_format:  Format of the crc length
        :param number_of_chunks_len_format: Format of the number of chunks length
        :param id_len_format: Format of the ID length
        :return: -1 if the decoding wasn't successful
        """
        decoded = False
        self.EOF = False
        if self.file.lower().endswith("dna"):
            try:
                self.f.close()
                self.f = quat_file_to_bin(self.file)
            except TypeError:
                print("skipping CORRUPT file - contains illegal character(s)")
                self.corrupt += 1
        if self.static_number_of_chunks is not None:
            self.number_of_chunks = self.static_number_of_chunks
            number_of_chunks_len_format = ""  # if we got static number_of_chunks we do not need it in struct string
        if self.file.lower().endswith("fasta"):
            self.f.close()
            self.f = open(self.file, "r")
            raw_packet_list = []
            while not (decoded or self.EOF):
                line = self.f.readline()
                if not line:
                    self.EOF = True
                    break
                try:
                    error_prob, seed = line[1:].replace("\n", "").split("_")
                except:
                    error_prob, seed = "0", "0"
                line = self.f.readline()
                if not line:
                    self.EOF = True
                    break
                dna_str = line.replace("\n", "")
                raw_packet_list.append((error_prob, seed, dna_str))
                try:
                    new_pack = self.parse_raw_packet(BytesIO(tranlate_quat_to_byte(dna_str)).read(),
                                                     crc_len_format=crc_len_format,
                                                     number_of_chunks_len_format=number_of_chunks_len_format,
                                                     packet_len_format=packet_len_format,
                                                     id_len_format=id_len_format)
                except Exception:
                    new_pack = "CORRUPT"
                if new_pack != "CORRUPT":
                    decoded = self.input_new_packet(new_pack)
                    if self.progress_bar is not None:
                        self.progress_bar.update(self.correct, Corrupt=self.corrupt)
        else:
            while not (decoded or self.EOF):
                new_pack = self.getNextValidPacket(False, packet_len_format=packet_len_format,
                                                   crc_len_format=crc_len_format,
                                                   number_of_chunks_len_format=number_of_chunks_len_format,
                                                   id_len_format=id_len_format)
                if new_pack is None:
                    break
                # koennte durch input_new_packet ersetzt werden:
                # self.addPacket(new_pack)
                if new_pack != "CORRUPT":
                    decoded = self.input_new_packet(new_pack)
                #
        print("Decoded Packets: " + str(self.correct))
        print("Corrupt Packets : " + str(self.corrupt))
        if self.GEPP.isPotentionallySolvable():
            return self.GEPP.solve()
        if not decoded and self.EOF:
            print("Unable to retrieve File from Chunks. Too much errors?")
            return -1
        # self.f.close()

    def getNumberOfLDPCBlocks(self):
        return self.s

    def getNumberOfHalfBlocks(self):
        return self.h

    def getNumberOfRepairBlocks(self):
        return self.getNumberOfHalfBlocks() + self.getNumberOfLDPCBlocks()

    def input_new_packet(self, packet: RU10Packet):
        """
        Removes auxpackets (LDPC and Half) and adds the remaining data to the GEPP matrix.
        :param packet: A Packet to add to the GEPP matrix
        :return: True: If solved. False: Else.
        """
        if self.ldpcANDhalf == dict() and self.distribution is None:  # self.isPseudo and
            self.distribution = RaptorDistribution(self.number_of_chunks)
            self.number_of_chunks = packet.get_total_number_of_chunks()
            _, self.s, self.h = intermediate_symbols(self.number_of_chunks, self.distribution)
            self.createAuxBlocks()
            self.progress_bar = self.create_progress_bar(self.number_of_chunks + 0.02 * self.number_of_chunks)
        # we need to do it twice sine half symbols may contain ldpc symbols (which by definition are repair codes.)
        if self.debug:
            print("----")
            print("Id = " + str(packet.id))
            print(packet.used_packets)
        removed = self.removeAndXorAuxPackets(packet)
        if self.debug:
            print(from_true_false_list(removed))
            print(packet.get_error_correction())
            print("----")
        if self.count:
            for i in range(len(removed)):
                if i in self.counter.keys():
                    if removed[i]:
                        self.counter[i] += 1
                else:
                    self.counter[i] = 1
        if self.GEPP is None:
            self.GEPP = GEPP(np.array([removed], dtype=bool), np.array([[packet.get_data()]], dtype=bytes), )
        else:
            self.GEPP.addRow(np.array(removed, dtype=bool), np.frombuffer(packet.get_data(), dtype="uint8"), )
        if (self.isPseudo or not self.read_all_before_decode) and self.GEPP.isPotentionallySolvable():
            # and self.GEPP.n % 5 == 0:  # Nur alle 5 Packete versuch starten
            if self.debug:
                print("current size: " + str(self.GEPP.n))
            return self.GEPP.solve(partial=False)
        return False

    # Correct
    def removeAndXorAuxPackets(self, packet: RU10Packet):
        """
        Removes auxpackets (LDCP and Half) from a given packet to get the packets data.
        :param packet: Packet to remove auxpackets from
        :return: The data without the auxpackets
        """
        aux_mapping = self.getHalfPacketListFromPacket(packet)  # Enthaelt Data + LDPC Nummern
        aux_mapping.append(packet.get_bool_array_used_and_ldpc_packets())
        xored_list = logical_xor(aux_mapping)
        del aux_mapping
        tmp = from_true_false_list(xored_list)  # Nur noch Data + LDPC sind vorhanden
        if self.debug:
            print(tmp)
        tmp = RU10Packet("", tmp, self.number_of_chunks, packet.id, packet.dist, read_only=True)
        aux_mapping = self.getAuxPacketListFromPacket(tmp)
        aux_mapping.append(tmp.get_bool_array_used_packets())  # [-len(self.auxBlocks):])
        res = logical_xor(aux_mapping)
        del tmp, aux_mapping
        return res

    def createAuxBlocks(self):
        """
        Reconstructs the auxblocks to be able to remove them afterwards.
        :return:
        """
        assert (self.number_of_chunks is not None), "createAuxBlocks can only be called AFTER first Packet"
        if self.debug: print("We should have " + str(self.getNumberOfLDPCBlocks()) + " LDPC-Blocks, " + str(
            self.getNumberOfHalfBlocks()) + " Half-Blocks and " + str(
            self.number_of_chunks) + " normal Chunks (including 1 HeaderChunk)")
        for i in range(0, self.getNumberOfRepairBlocks()):
            self.repairBlockNumbers[i] = set()
        i = 0
        for group in self.generateIntermediateBlocksFormat(self.number_of_chunks):
            for elem in group:
                self.repairBlockNumbers[i] = elem
                i += 1
        # XOR all Chunks into the corresponding AUX-Block
        for aux_number in self.repairBlockNumbers.keys():
            self.ldpcANDhalf[aux_number] = RU10IntermediatePacket("", self.repairBlockNumbers[aux_number],
                                                                  total_number_of_chunks=self.number_of_chunks,
                                                                  id=aux_number, dist=self.distribution)
            if self.debug:
                print(str(aux_number) + " : " + str(self.ldpcANDhalf[aux_number].used_packets))

    # Correct
    def getAuxPacketListFromPacket(self, packet: RU10Packet):
        """
        Creates a list for a packet with information about whether auxpackets have been used for that packet.
        :param packet: The packet to check.
        :return: Information about used auxpackets.
        """
        res = []
        aux_used_packets = packet.get_bool_array_repair_packets()
        for i in range(len(aux_used_packets)):
            if aux_used_packets[i]:
                res.append((self.ldpcANDhalf[i].get_bool_array_used_packets()))
        return res

    def getHalfPacketListFromPacket(self, packet: RU10Packet) -> typing.List[typing.List[bool]]:
        """
        Generates a list of halfpackets from a packet.
        :param packet: The packet to get the list from
        :return: List of halfpackets
        """
        res: typing.List[typing.List[bool]] = []
        aux_used_packets = packet.get_bool_array_half_packets()
        for i in range(len(aux_used_packets)):
            if aux_used_packets[i]:
                res.append(
                    (self.ldpcANDhalf[packet.get_number_of_ldpc_blocks() + i].get_bool_array_used_and_ldpc_packets()))
        return res

    def solve(self, partial=False) -> bool:
        """
        Calls GEPP.solve()
        :return: True: If GEPP was able to solve the matrix. False: Else.
        """
        return self.GEPP.solve(partial=partial)

    def getSolvedCount(self) -> int:
        return self.GEPP.getSolvedCount()

    def is_decoded(self) -> bool:
        """
        Checks if the data is decoded.
        :return: True: If the decoding was successfull. False: Else.
        """
        return self.GEPP is not None and self.GEPP.isPotentionallySolvable() and self.GEPP.isSolved()

    def getNextValidPacket(self, from_multiple_files: bool = False, packet_len_format: str = "I",
                           crc_len_format: str = "L", number_of_chunks_len_format: str = "I",
                           id_len_format: str = "I") -> typing.Optional[RU10Packet]:
        """
        Takes a raw packet from a file and calls @parse_raw_packet to get a RU10 packet. If the packet is corrupt the
        next one will be taken.
        :param from_multiple_files: True: The packets were saved in multiple files. False: Packets were saved in one file.
        :param packet_len_format: Format of the packet length
        :param crc_len_format:  Format of the crc length
        :param number_of_chunks_len_format: Format of the number of chunks length
        :param id_len_format: Format of the ID length
        :return: RU10Packet
        """
        if not from_multiple_files:
            packet_len = self.f.read(struct.calcsize("<" + packet_len_format))
            try:
                packet_len = struct.unpack("<" + packet_len_format, packet_len)[0]
                packet = self.f.read(int(packet_len))
            except:
                return None
        else:
            packet = self.f.read()
            packet_len = len(packet)
        if not packet or not packet_len:  # EOF
            self.EOF = True
            try:
                self.f.close()
            except:
                return None
            return None
        res = self.parse_raw_packet(packet, crc_len_format=crc_len_format,
                                    packet_len_format=packet_len_format,
                                    number_of_chunks_len_format=number_of_chunks_len_format,
                                    id_len_format=id_len_format)
        if res == "CORRUPT" and not from_multiple_files:
            res = self.getNextValidPacket(from_multiple_files, packet_len_format=packet_len_format,
                                          crc_len_format=crc_len_format,
                                          number_of_chunks_len_format=number_of_chunks_len_format,
                                          id_len_format=id_len_format)
        return res

    def parse_raw_packet(self, packet, crc_len_format: str = "L", number_of_chunks_len_format: str = "L",
                         packet_len_format: str = "I", id_len_format: str = "L") -> typing.Union[RU10Packet, str]:
        """
        Creates a RU10 packet from a raw given packet. Also checks if the packet is corrupted. If any method was used to
        create packets from specific chunks, set self.use_method = True. This will treat the last byte of the raw packet
        data as the byte that contains the information about the used method ("even", "odd", "window_30 + window" or
        "window_40 + window". See RU10Encoder.create_new_packet_from_chunks for further information.
        :param packet: A raw packet
        :param packet_len_format: Format of the packet length
        :param crc_len_format:  Format of the crc length
        :param number_of_chunks_len_format: Format of the number of chunks length
        :param id_len_format: Format of the ID length
        :return: RU10Packet or an error message
        """
        struct_str = "<" + number_of_chunks_len_format + id_len_format
        struct_len = struct.calcsize(struct_str)
        try:
            packet = self.error_correction(packet)
        except:
            self.corrupt += 1
            return "CORRUPT"

        data = packet[struct_len:]
        chunk_lst = []
        if self.use_method:
            method_data = bin(data[-1])[2:]
            while len(method_data) < 8:
                method_data = '0' + method_data
            data = data[:-1]
            if method_data.startswith('00'):
                chunk_lst = [ch for ch in range(0, self.number_of_chunks + 1) if ch % 2 == 0]
            elif method_data.startswith('01'):
                chunk_lst = [ch for ch in range(0, self.number_of_chunks + 1) if ch % 2 != 0]
            elif method_data.startswith('10'):
                window = int(method_data[2:], 2)
                window_size = 30
                start = window * (window_size - 10)
                chunk_lst = [ch for ch in range(start, start + window_size) if ch <= self.number_of_chunks]
            elif method_data.startswith('11'):
                window = int(method_data[2:], 2)
                window_size = 40
                start = window * (window_size - 10)
                chunk_lst = [ch for ch in range(start, start + window_size) if ch <= self.number_of_chunks]
            else:
                raise RuntimeError("Not a valid start:", method_data)
        len_data = struct.unpack(struct_str, packet[0:struct_len])
        if self.static_number_of_chunks is None:
            self.number_of_chunks = xor_mask(len_data[0], number_of_chunks_len_format)
            unxored_id = xor_mask(len_data[1], id_len_format)
        else:
            unxored_id = xor_mask(len_data[0], id_len_format)
        if self.distribution is None:
            self.distribution = RaptorDistribution(self.number_of_chunks)
            _, self.s, self.h = intermediate_symbols(self.number_of_chunks, self.distribution)
            self.progress_bar = self.create_progress_bar(self.number_of_chunks + 0.02 * self.number_of_chunks)

        if self.correct == 0:
            self.createAuxBlocks()
        self.correct += 1
        if self.use_method:
            numbers = choose_packet_numbers(len(chunk_lst), unxored_id, self.distribution, systematic=False,
                                            max_l=len(chunk_lst))
            used_packets = [chunk_lst[i] for i in numbers]
        else:
            used_packets = choose_packet_numbers(self.number_of_chunks, unxored_id, self.distribution, systematic=False)
        res = RU10Packet(data, used_packets, self.number_of_chunks, unxored_id, read_only=True,
                         packet_len_format=packet_len_format, crc_len_format=crc_len_format,
                         number_of_chunks_len_format=number_of_chunks_len_format, id_len_format=id_len_format,
                         save_number_of_chunks_in_packet=self.static_number_of_chunks is None)
        return res

    def generateIntermediateBlocksFormat(self, number_of_chunks: int) -> typing.List[typing.List[typing.List[int]]]:
        """
        Generates the format of the intermediate blocks from the number of used chunks.
        :param number_of_chunks: The number of used chunks.
        :return:
        """
        compositions: typing.List[typing.List[int]] = [[] for _ in range(self.s)]
        for i in range(0, number_of_chunks):
            a = 1 + (int(floor(np.float64(i) / np.float64(self.s))) % (self.s - 1))
            b = int(i % self.s)
            compositions[b].append(i)
            b = (b + a) % self.s
            compositions[b].append(i)
            b = (b + a) % self.s
            compositions[b].append(i)

        hprime: int = int(ceil(np.float64(self.h) / 2))
        m = buildGraySequence(number_of_chunks + self.s, hprime)
        hcompositions: typing.List[typing.List[int]] = [[] for _ in range(self.h)]
        for i in range(0, self.h):
            hcomposition = []
            for j in range(0, number_of_chunks + self.s):
                if bitSet(np.uint32(m[j]), np.uint32(i)):
                    hcomposition.append(j)
            hcompositions[i] = hcomposition
        res = [compositions, hcompositions]
        return res

    def saveDecodedFile(self, last_chunk_len_format: str = "I", null_is_terminator: bool = False,
                        print_to_output: bool = True, return_file_name=False, partial_decoding: bool = True) -> \
            typing.Union[bytes, str]:
        """
        Saves the file - if decoded. The filename is either taken from the headerchunk or generated based on the input
        filename.
        :param partial_decoding: perform partial decoding if full decoding failed, missing parts will be filled with "\x00"
        :param return_file_name: if set to true, this function will return the filename under which the file as been saved
        :param last_chunk_len_format: Format of the last chunk length
        :param null_is_terminator: True: The file is handled as null-terminated C-String.
        :param print_to_output: True: Result we be printed to the command line.
        :return:
        """
        assert self.is_decoded() or partial_decoding, "Can not save File: Unable to reconstruct. You may try saveDecodedFile(partial_decoding=True)"
        if partial_decoding:
            self.solve(partial=True)
        dirty = False
        if self.use_headerchunk:
            header_row = self.GEPP.result_mapping[0]
            if header_row >= 0:
                self.headerChunk = HeaderChunk(
                    Packet(self.GEPP.b[header_row], {0}, self.number_of_chunks, read_only=True),
                    last_chunk_len_format=last_chunk_len_format)
        file_name = "DEC_" + os.path.basename(self.file) if self.file is not None else "RU10.BIN"
        output_concat = b""
        if self.headerChunk is not None:
            try:
                file_name = self.headerChunk.get_file_name().decode("utf-8")
            except Exception as ex:
                print("Warning:", ex)
        file_name = file_name.split("\x00")[0]
        with open(file_name, "wb") as f:
            for x in self.GEPP.result_mapping:
                if x < 0:
                    f.write(b"\x00" * len(self.GEPP.b[x][0]))
                    dirty = True
                    continue
                if 0 != x or not self.use_headerchunk:
                    if self.number_of_chunks - 1 == x and self.use_headerchunk:
                        output = self.GEPP.b[x][0][0: self.headerChunk.get_last_chunk_length()]
                        output_concat += output.tobytes()
                        f.write(output)
                    else:
                        if null_is_terminator:
                            splitter = self.GEPP.b[x].tostring().decode().split("\x00")
                            output = splitter[0].encode()
                            if type(output) == bytes:
                                output_concat += output
                            else:
                                output_concat += output.tobytes()
                            f.write(output)
                            if len(splitter) > 1:
                                break  # since we are in null-terminator mode, we exit once we see the first 0-byte
                        else:
                            output = self.GEPP.b[x]
                            if type(output) == bytes:
                                output_concat += output
                            else:
                                output_concat += output.tobytes()
                            f.write(output)
        print("Saved file as '" + str(file_name) + "'")
        if dirty:
            print("Some parts could not be restored, file WILL contain sections with \\x00 !")
        if print_to_output:
            print("Result:")
            print(output_concat.decode("utf-8"))
        if self.progress_bar is not None:
            self.progress_bar.update(self.number_of_chunks, Corrupt=self.corrupt)
        if return_file_name:
            return file_name
        return output_concat

    def mode_1_bmp_decode(self, last_chunk_len_format: str = "I"):
        dec_out = self.saveDecodedFile(last_chunk_len_format=last_chunk_len_format, null_is_terminator=False,
                                       print_to_output=False)
        return self.bytes_to_bitmap(dec_out)

    def bytes_to_bitmap(self, img_byt: bytes):
        width, height = struct.unpack('>H', img_byt[:2])[0], struct.unpack('>H', img_byt[2:4])[0]
        unpack = np.unpackbits(
            np.frombuffer(img_byt, dtype=np.uint8, count=int((width * height) / 8), offset=4)).reshape(height,
                                                                                                       width).transpose()
        flip_bits = np.logical_not(unpack).astype(int)
        new_img = self.draw_img(flip_bits, width, height)
        tmp_file_name = os.path.basename(self.file) + ".bmp"
        file_name = "DEC_" + tmp_file_name if self.file is not None else "RU10.BIN.bmp"
        new_img.save(file_name)
        return file_name

    @staticmethod
    def draw_img(unpacked_flipped_bits, width: int, height: int) -> Image:
        new_img = Image.new('1', (width, height))
        pixels = new_img.load()

        for i in range(new_img.size[0]):
            for j in range(new_img.size[1]):
                pixels[i, j] = int(unpacked_flipped_bits[i, j])
        return new_img


def main(file: str, number_of_chunks: int, error_correction: typing.Callable = nocode, insert_header: bool = False,
         mode_1_bmp: bool = False):
    print("Pure Gauss-Mode")
    x = RU10Decoder(file, use_headerchunk=insert_header, error_correction=error_correction,
                    static_number_of_chunks=number_of_chunks)
    x.decode(id_len_format="I", number_of_chunks_len_format="I")
    x.saveDecodedFile(null_is_terminator=False, print_to_output=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "filename", metavar="file", type=str, help="the file / folder to Decode"
    )
    parser.add_argument("--error_correction", metavar="error_correction", type=str, required=False,
                        default="nocode",
                        help="Error Correction Method to use; possible values: \
                            nocode, crc, reedsolomon, dna_reedsolomon (default=nocode)")
    parser.add_argument("--insert_header", required=False, action="store_true", default=False)
    parser.add_argument("--number_of_chunks", metavar="number_of_chunks", required=True, type=int)
    parser.add_argument("--repair_symbols", metavar="repair_symbols", type=int, required=False, default=2,
                        help="number of repair symbols for ReedSolomon (default=2)")
    parser.add_argument("--as_mode_1_bmp", required=False, action="store_true",
                        help="convert to a header-less B/W BMP format. (use only for image/bmp input)")
    args = parser.parse_args()
    _file = args.filename
    _repair_symbols = args.repair_symbols
    _insert_header = args.insert_header
    _mode_1_bmp = args.as_mode_1_bmp
    _number_of_chunks = args.number_of_chunks
    _error_correction = get_error_correction_decode(args.error_correction, _repair_symbols)
    print("File / Folder to decode: " + str(_file))
    main(_file, _number_of_chunks, _error_correction, _insert_header, _mode_1_bmp)
    print("Decoding finished.")
