#!/usr/bin/python
# -*- coding: latin-1 -*-
import argparse
import struct
import os
import typing
import numpy as np
from io import BytesIO

from norec4dna.Decoder import Decoder
from norec4dna.ErrorCorrection import crc32, nocode, reed_solomon_decode
from norec4dna.GEPP import GEPP
from norec4dna.HeaderChunk import HeaderChunk
from norec4dna.Packet import Packet
from norec4dna.distributions.Distribution import Distribution
from norec4dna.distributions.ErlichZielinskiRobustSolitonDisribution import ErlichZielinskiRobustSolitonDistribution
from norec4dna.helper import calc_crc, xor_mask
from norec4dna.helper.quaternary2Bin import quat_file_to_bin, tranlate_quat_to_byte


class LTDecoder(Decoder):
    def __init__(self, file: typing.Optional[str] = None, error_correction: typing.Callable = nocode,
                 use_headerchunk: bool = True, static_number_of_chunks: typing.Optional[int] = None,
                 implicit_mode: bool = True, dist: typing.Optional[Distribution] = None):
        super().__init__(file)
        self.use_headerchunk: bool = use_headerchunk
        self.isPseudo: bool = False
        self.file: typing.Optional[str] = file
        self.degreeToPacket: typing.Dict[int, Packet] = {}
        if file is not None:
            self.isFolder: bool = os.path.isdir(file)
            if not self.isFolder:
                self.f = open(self.file, "rb")
        self.correct: int = 0
        self.corrupt: int = 0
        self.number_of_chunks: int = 1000000
        self.headerChunk: typing.Optional[HeaderChunk] = None
        self.GEPP: typing.Optional[GEPP] = None
        self.pseudoCount: int = 0
        self.read_all_before_decode: bool = True
        self.count: bool = True
        self.counter: typing.Dict[int, int] = dict()
        self.error_correction: typing.Callable = error_correction
        self.static_number_of_chunks: int = static_number_of_chunks
        self.implicit_mode: bool = implicit_mode
        self.dist: typing.Optional[Distribution] = dist
        self.EOF: bool = False

    def decodeFolder(self, packet_len_format: str = "I", crc_len_format: str = "L",
                     number_of_chunks_len_format: str = "I", degree_len_format: str = "I", seed_len_format: str = "I",
                     last_chunk_len_format: str = "I") -> typing.Optional[int]:
        decoded: bool = False
        self.EOF: bool = False
        if self.static_number_of_chunks is not None:
            self.number_of_chunks = self.static_number_of_chunks
            number_of_chunks_len_format = ""  # if we got static number_of_chunks we do not need it in struct string
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
                    decoded = self.input_new_packet(new_pack)
                if decoded:
                    break
        self.EOF = True
        print("Decoded Packets: " + str(self.correct))
        print("Corrupt Packets : " + str(self.corrupt))
        if self.GEPP.isPotentionallySolvable():
            decoded = self.GEPP.solve()
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
                new_pack = self.parse_raw_packet(BytesIO(tranlate_quat_to_byte(dna_str)).read(),
                                                 crc_len_format=crc_len_format,
                                                 number_of_chunks_len_format=number_of_chunks_len_format,
                                                 degree_len_format=degree_len_format,
                                                 seed_len_format=seed_len_format)
                decoded = self.input_new_packet(new_pack)
                if self.progress_bar is not None:
                    self.progress_bar.update(self.correct, Corrupt=self.corrupt)
            else:
                while not (decoded or self.EOF):
                    new_pack = self.getNextValidPacket(False, packet_len_format=packet_len_format,
                                                       crc_len_format=crc_len_format,
                                                       number_of_chunks_len_format=number_of_chunks_len_format,
                                                       degree_len_format=degree_len_format,
                                                       seed_len_format=seed_len_format,
                                                       last_chunk_len_format=last_chunk_len_format)
                    if new_pack is None:
                        break
                    # koennte durch input_new_packet ersetzt werden:
                    # self.addPacket(new_pack)
                    decoded = self.input_new_packet(new_pack)
                    ##
        print("Decoded Packets: " + str(self.correct))
        print("Corrupt Packets : " + str(self.corrupt))
        self.f.close()
        if self.GEPP.isPotentionallySolvable():
            return self.GEPP.solve()
        if not decoded and self.EOF:
            print("Unable to retrieve File from Chunks. Too much errors?")
            return -1

    def input_new_packet(self, packet: Packet) -> bool:
        self.pseudoCount += 1
        packets: typing.List[bool] = packet.get_bool_array_used_packets()
        if self.count:
            for i in range(len(packets)):
                if i in self.counter.keys():
                    if packets[i]:
                        self.counter[i] += 1
                else:
                    self.counter[i] = 1
        if self.GEPP is None:
            self.GEPP = GEPP(np.array([packet.get_bool_array_used_packets()], dtype=bool),
                             np.array([[packet.get_data()]], dtype=bytes), )
        else:
            self.GEPP.addRow(packet.get_bool_array_used_packets(), np.frombuffer(packet.get_data(), dtype="uint8"), )
        if self.isPseudo and not self.read_all_before_decode and self.GEPP.isPotentionallySolvable():
            return self.GEPP.solve(partial=False)
        return False

    def solve(self) -> bool:
        return self.GEPP.solve()

    def getSolvedCount(self) -> int:
        return self.GEPP.getSolvedCount()

    def choose_packet_numbers(self, degree: int, seed: int = 0) -> typing.Set[int]:
        assert degree <= self.number_of_chunks
        res: typing.Set[int] = set()
        rng = np.random
        rng.seed(seed)
        for _ in range(0, degree):
            tmp = rng.choice(range(0, self.number_of_chunks))
            while tmp in res:
                tmp = rng.choice(range(0, self.number_of_chunks))
            res.add(tmp)
        return res

    def is_decoded(self) -> bool:
        return self.GEPP is not None and self.GEPP.isPotentionallySolvable() and self.GEPP.isSolved()

    def getNextValidPacket(self, from_multiple_files: bool = False, packet_len_format: str = "I",
                           crc_len_format: str = "L", number_of_chunks_len_format: str = "I",
                           degree_len_format: str = "I", seed_len_format: str = "I",
                           last_chunk_len_format: str = "I") -> typing.Optional[Packet]:
        if not from_multiple_files:
            packet_len: typing.Union[int, bytes] = self.f.read(struct.calcsize("<" + packet_len_format))
            packet_len = struct.unpack("<" + packet_len_format, packet_len)[0]
            packet: bytes = self.f.read(int(packet_len))
        else:
            packet = self.f.read()
            packet_len = len(packet)
        if not packet or not packet_len:  # EOF
            self.EOF: bool = True
            self.f.close()
            return None
        res = self.parse_raw_packet(packet, number_of_chunks_len_format=number_of_chunks_len_format,
                                    degree_len_format=degree_len_format,
                                    seed_len_format=seed_len_format)
        if res == "CORRUPT":
            res = self.getNextValidPacket(from_multiple_files=from_multiple_files,
                                          number_of_chunks_len_format=number_of_chunks_len_format,
                                          degree_len_format=degree_len_format,
                                          seed_len_format=seed_len_format)
        return res

    def saveDecodedFile(self, last_chunk_len_format: str = "I", null_is_terminator: bool = False,
                        print_to_output: bool = True) -> None:
        assert self.is_decoded(), "Can not save File: Unable to reconstruct."
        if self.use_headerchunk:
            self.headerChunk = HeaderChunk(Packet(self.GEPP.b[0], {0}, self.number_of_chunks, read_only=True),
                                           last_chunk_len_format=last_chunk_len_format)
        file_name = "DEC_" + os.path.basename(self.file) if self.file is not None else "LT.BIN"
        if self.headerChunk is not None:
            file_name = self.headerChunk.get_file_name().decode("utf-8")
        output_concat: bytes = b""
        file_name: str = file_name.split("\x00")[0]
        try:
            with open(file_name, "wb") as f:
                for x in self.GEPP.result_mapping:
                    if 0 != x or not self.use_headerchunk:
                        if self.number_of_chunks - 1 == x and self.use_headerchunk:
                            output: typing.Union[bytes, np.array] = self.GEPP.b[x][0][
                                                                    0: self.headerChunk.get_last_chunk_length()]
                            output_concat += output.tobytes()
                            f.write(output)
                        else:
                            if null_is_terminator:
                                splitter: str = self.GEPP.b[x].tostring().decode().split("\x00")
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
                                output_concat += output.tobytes()
                                f.write(output)
            print("Saved file as '" + str(file_name) + "'")
        except Exception as ex:
            raise ex
        if print_to_output:
            print("Result:")
            print(output_concat.decode("utf-8"))

    def parse_raw_packet(self, packet: bytes, crc_len_format: str = "L", number_of_chunks_len_format: str = "I",
                         degree_len_format: str = "I", seed_len_format: str = "I") -> typing.Union[str, Packet]:
        crc_len = -struct.calcsize("<" + crc_len_format)
        if self.error_correction.__code__.co_name == crc32.__code__.co_name:
            payload: bytes = packet[:crc_len]
            crc: int = struct.unpack("<" + crc_len_format, packet[crc_len:])[0]
            calced_crc: int = calc_crc(payload)
            if crc != calced_crc:  # If the Packet is corrupt, try next one
                print("[-] CRC-Error - " + str(hex(crc)) + " != " + str(hex(calced_crc)))
                self.corrupt += 1
                return "CORRUPT"
        else:
            crc_len = None
            try:
                packet = self.error_correction(packet)
            except:
                self.corrupt += 1
                return "CORRUPT"
        if self.implicit_mode:
            degree_len_format = ""
        struct_str: str = "<" + number_of_chunks_len_format + degree_len_format + seed_len_format
        struct_len: int = struct.calcsize(struct_str)
        len_data: typing.Union[int, typing.Tuple[int, int], typing.Tuple[int, int, int]] = struct.unpack(struct_str,
                                                                                                         packet[
                                                                                                         0:struct_len])
        degree: typing.Optional[int] = None
        if self.static_number_of_chunks is None:
            if self.implicit_mode:
                number_of_chunks, seed = len_data
            else:
                number_of_chunks, degree, seed = len_data
            self.number_of_chunks = xor_mask(number_of_chunks, number_of_chunks_len_format)
        else:
            if self.implicit_mode:
                seed = len_data
            else:
                degree, seed = len_data
        seed: int = xor_mask(seed, seed_len_format)
        if degree is None:
            self.dist.set_seed(seed)
            degree: int = self.dist.getNumber()
        else:
            degree: int = xor_mask(degree, degree_len_format)
        used_packets = self.choose_packet_numbers(degree, seed)
        data = packet[struct_len:crc_len]
        self.correct += 1

        return Packet(data, used_packets, self.number_of_chunks, read_only=True, error_correction=self.error_correction,
                      save_number_of_chunks_in_packet=self.static_number_of_chunks is None)


def main(file: str, number_of_chunks: int, error_correction: typing.Callable, insertheader: bool):
    dist = ErlichZielinskiRobustSolitonDistribution(number_of_chunks, seed=2)

    decoder = LTDecoder(file, error_correction=error_correction, use_headerchunk=insertheader,
                        static_number_of_chunks=number_of_chunks, implicit_mode=False, dist=dist)
    decoder.decode(number_of_chunks_len_format="H", seed_len_format="I", degree_len_format="H")
    decoder.saveDecodedFile(null_is_terminator=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", metavar="file", type=str, help="the file / folder to Decode")
    parser.add_argument("--error_correction", metavar="error_correction", type=str, required=False,
                        default="nocode",
                        help="Error Correction Method to use; possible values: \
                                    nocode, crc, reedsolomon (default=nocode)")
    parser.add_argument("--repair_symbols", metavar="repair_symbols", type=int, required=False, default=2,
                        help="number of repair symbols for ReedSolomon (default=2)")
    parser.add_argument("--insert_header", metavar="insert_header", required=False, type=bool, default=False)
    parser.add_argument("--number_of_chunks", metavar="number_of_chunks", required=True, type=int)
    args = parser.parse_args()
    filename = args.filename
    e_correction_str = args.error_correction
    _repair_symbols = args.repair_symbols
    _insert_header = args.insert_header
    _number_of_chunks = args.number_of_chunks
    if e_correction_str == "nocode":
        e_correction = nocode
    elif e_correction_str == "crc":
        e_correction = crc32
    elif e_correction_str == "reedsolomon":
        if _repair_symbols != 2:
            e_correction = lambda x: reed_solomon_decode(x, _repair_symbols)
        else:
            e_correction = reed_solomon_decode
    else:
        print("Selected Error Correction not supported, choose: 'nocode', 'crc' or 'reedsolomon'")
        e_correction = None
        exit()
    print("File / Folder to decode: " + str(filename))
    main(filename, _number_of_chunks, e_correction, _insert_header)
    print("Decoding finished.")
