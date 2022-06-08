#!/usr/bin/python
# -*- coding: latin-1 -*-
import argparse
import configparser
import datetime
import glob
import struct
import time
import os
import numpy as np
import typing
from math import ceil

from norec4dna.Decoder import Decoder
from norec4dna.Encoder import Encoder
from norec4dna.ErrorCorrection import nocode, crc32, reed_solomon_encode
from norec4dna.Packet import Packet
from norec4dna.distributions.Distribution import Distribution
from norec4dna.distributions.ErlichZielinskiRobustSolitonDisribution import ErlichZielinskiRobustSolitonDistribution
from norec4dna.helper import should_drop_packet, listXOR
from norec4dna.rules.DNARules import DNARules
from norec4dna.rules.DNARules2 import DNARules2
from norec4dna.rules.DNARules_ErlichZielinski import DNARules_ErlichZielinski
from norec4dna.rules.FastDNARules import FastDNARules


class LTEncoder(Encoder):
    def __init__(self, file: str, number_of_chunks: int, distribution: Distribution, insert_header: bool = True,
                 pseudo_decoder: typing.Optional[Decoder] = None, prioritized_packets=None,
                 chunk_size: int = 0, error_correction: typing.Callable = nocode, rules: typing.Optional[
                typing.Union[DNARules, DNARules2, FastDNARules, DNARules_ErlichZielinski]] = None,
                 implicit_mode: bool = True, packet_len_format: str = "I", crc_len_format: str = "L",
                 number_of_chunks_len_format: str = "I", used_packets_len_format: str = "I", id_len_format: str = "I",
                 last_chunk_len_format: str = "I", save_number_of_chunks_in_packet: bool = True, drop_upper_bound=1.0,
                 sequential_seed=True):
        super().__init__(file, number_of_chunks, distribution, insert_header, pseudo_decoder, chunk_size)
        if prioritized_packets is None:
            prioritized_packets = []
        assert number_of_chunks == distribution.get_size()
        self.out_file: typing.Optional[str] = None
        self.dist: Distribution = distribution
        self.rng: np.random = np.random
        self.insert_header: bool = insert_header
        self.chunk_size: int = chunk_size
        self.file: str = file
        self.rules: typing.Optional = rules
        self.upper_bound: float = drop_upper_bound
        if self.chunk_size == 0:
            self.number_of_chunks: int = number_of_chunks
        else:
            self.set_no_chunks_from_chunk_size()
        self.chunks: typing.List[bytes] = []
        self.overhead_limit: float = 0.20
        self.encodedPackets: typing.Set[Packet] = set()
        self.setOfEncodedPackets: typing.Set[int] = set()
        self.pseudo_decoder: typing.Optional[Decoder] = pseudo_decoder
        self.prioritized_packets: typing.List = prioritized_packets
        self.error_correction: typing.Callable = error_correction
        self.implicit_mode: bool = implicit_mode
        # Struct-Strings:
        self.packet_len_format: str = packet_len_format
        self.crc_len_format: str = crc_len_format
        self.number_of_chunks_len_format: str = number_of_chunks_len_format
        self.used_packets_len_format: str = used_packets_len_format
        self.id_len_format: str = id_len_format
        self.last_chunk_len_format: str = last_chunk_len_format
        self.save_number_of_chunks_in_packet: bool = save_number_of_chunks_in_packet
        self.ruleDrop: int = 0
        self.next_checkblock_id = -1
        self.sequential_seed = sequential_seed

    def encode_file(self, split_to_multiple_files: bool = False):
        self.encode_to_packets()
        self.save_packets(split_to_multiple_files)

    def prepareEncoder(self):
        self.ruleDrop = 0
        file_size: int = self.get_file_size(self.file)
        if self.insert_header:
            self.chunk_size = ceil(1.0 * file_size / (self.number_of_chunks - 1))
            self.chunks = self.create_chunks(self.chunk_size)
            self.number_of_chunks += 1
            # First Chunk is a Header
            self.chunks.insert(0, self.encode_header_info())
        else:
            self.chunk_size = ceil(1.0 * file_size / self.number_of_chunks)
            self.chunks = self.create_chunks(self.chunk_size)
        self.dist.update_number_of_chunks(self.number_of_chunks)
        self.fill_last_chunk()

    def encode_to_packets(self) -> float:
        self.prepareEncoder()
        start: float = time.time()
        if self.pseudo_decoder is not None:
            while not self.pseudo_decoder.is_decoded():
                # This process continues until the receiver signals that the
                # message has been received and successfully decoded.
                new_pack: Packet = self.create_new_packet()
                if self.rules is not None:
                    while should_drop_packet(self.rules, new_pack, self.upper_bound):
                        new_pack = self.create_new_packet()
                        self.ruleDrop += 1

                self.pseudo_decoder.input_new_packet(new_pack)
                self.encodedPackets.add(new_pack)
        else:
            while (len(self.encodedPackets) < (self.number_of_chunks + (self.number_of_chunks * self.overhead_limit))
                   or self.number_of_packets_encoded_already() < self.number_of_chunks):  # 20% Aufschlag
                # This process continues until the receiver signals that the
                # message has been received and successfully decoded.
                pack: Packet = self.create_new_packet()
                if self.rules is not None:
                    while should_drop_packet(self.rules, pack, self.upper_bound):
                        pack = self.create_new_packet()
                        self.ruleDrop += 1
                self.encodedPackets.add(pack)
        return start

    def set_overhead_limit(self, n: float):
        self.overhead_limit = n

    def encodePriotizedPackets(self, error_correction: typing.Callable = nocode):
        for num in self.prioritized_packets:
            new_pack = Packet(self.chunks[num], {num}, self.number_of_chunks, read_only=False,
                              implicit_mode=self.implicit_mode, error_correction=error_correction,
                              packet_len_format=self.packet_len_format, crc_len_format=self.crc_len_format,
                              number_of_chunks_len_format=self.number_of_chunks_len_format,
                              used_packets_len_format=self.used_packets_len_format, id_len_format=self.id_len_format)
            if self.pseudo_decoder is not None:
                self.pseudo_decoder.input_new_packet(new_pack)
            self.encodedPackets.add(new_pack)

    def encode_header_info(self) -> bytes:
        # Size of last Chunk
        # Filename
        # PAD-Bytes
        last_chunk: bytes = self.chunks[-1]
        file_name_length: int = len(self.file)
        assert file_name_length + 4 < self.chunk_size, "Chunks too small for HeaderInfo"
        # -4 for bytes to store length of last_chunk (I)
        last_chunk_len_struct_size = struct.calcsize("<" + self.last_chunk_len_format)
        struct_str = "<" + self.last_chunk_len_format + str(file_name_length) + "s" + \
                     str(self.chunk_size - file_name_length - last_chunk_len_struct_size) + "x"
        return struct.pack(struct_str, len(last_chunk), bytes(self.file, encoding="utf-8"))

    def fill_last_chunk(self):
        last: bytes = self.chunks[-1]
        assert (len(last) <= self.chunk_size), "Error, last Chunk ist bigger than ChunkSize"
        if len(last) < self.chunk_size:
            struct_str = ("<" + str(len(last)) + "s" + str(self.chunk_size - len(last)) + "x")
            self.chunks[-1] = struct.pack(struct_str, bytes(last))

    def number_of_packets_encoded_already(self) -> int:
        return len(self.setOfEncodedPackets)

    def save_packets(self, split_to_multiple_files: bool, out_file: typing.Optional[int] = None,
                     save_as_dna: bool = False) -> None:
        fileending: str = ".LT" + ("_DNA" if save_as_dna else "")
        if not split_to_multiple_files:
            if out_file is None:
                out_file = self.file + fileending
            with open(out_file, "wb" if not save_as_dna else "w") as f:
                for packet in self.encodedPackets:
                    f.write(packet.get_dna_struct(split_to_multiple_files) if save_as_dna else packet.get_struct(
                        split_to_multiple_files))
        else:
            # Folder:
            if out_file is None:
                full_dir, file_name = os.path.split(os.path.realpath(self.file))
                file_name = "LT_" + file_name
                out_file = os.path.join(full_dir, file_name)
                if not out_file.endswith("/"):
                    files = glob.glob(out_file + "/*")
                else:
                    files = glob.glob(out_file + "*")
                for f in files:
                    os.remove(f)
            i = 0
            if not os.path.exists(out_file):
                os.makedirs(out_file)
            for packet in self.encodedPackets:
                with open(
                        out_file + "/" + str(packet.error_prob) + "_" + str(i) + fileending,
                        "wb" if not save_as_dna else "w",
                ) as f:
                    f.write(
                        packet.get_dna_struct(split_to_multiple_files)
                        if save_as_dna
                        else packet.get_struct(split_to_multiple_files)
                    )
                i += 1
            self.out_file = os.path.relpath(out_file)
            print("Config: " + self.getConfigStr(out_file))

    def get_file_size(self, file: str) -> int:
        return os.stat(file).st_size

    def generate_new_checkblock_id(self, sequential=True) -> int:
        max_num: float = Encoder.calc_max_size(struct.calcsize("<" + self.id_len_format))
        if sequential:
            self.next_checkblock_id += 1
            if self.next_checkblock_id > max_num:
                raise RuntimeError("sequential checkblock_id > max allowed number!")
            return self.next_checkblock_id
        random_generator: np.random.RandomState = np.random.RandomState()
        return random_generator.randint(0, max_num, dtype=np.uint32)

    def create_new_packet(self, seed=None) -> Packet:
        if seed is None:
            seed: int = self.generate_new_checkblock_id(self.sequential_seed)
        if self.implicit_mode:  # in implicit mode we want to be able derive the used chunks from having only the seed
            self.dist.set_seed(seed)
        degree: int = self.dist.getNumber()

        packet_numbers: typing.Set[int] = self.choose_packet_numbers(degree, seed=seed)
        chunks: typing.List[bytes] = [self.chunks[i] for i in packet_numbers]
        self.setOfEncodedPackets |= set(packet_numbers)
        return Packet(
            listXOR(chunks),
            packet_numbers,
            self.number_of_chunks,
            read_only=False,
            seed=seed,
            error_correction=self.error_correction,
            implicit_mode=self.implicit_mode,
            packet_len_format=self.packet_len_format,
            crc_len_format=self.crc_len_format,
            number_of_chunks_len_format=self.number_of_chunks_len_format,
            used_packets_len_format=self.used_packets_len_format, id_len_format=self.id_len_format,
            save_number_of_chunks_in_packet=self.save_number_of_chunks_in_packet
        )

    def create_and_add_new_packet(self, error_correction=nocode):
        packet = self.create_new_packet()
        self.encodedPackets.add(packet)
        return packet

    def choose_packet_numbers(self, degree: int, seed: int = 0) -> typing.Set[int]:
        assert degree <= len(self.chunks)
        res: typing.Set[int] = set()
        self.rng.seed(seed)
        for _ in range(0, degree):
            tmp = self.rng.choice(range(0, len(self.chunks)))
            while tmp in res:
                tmp = self.rng.choice(range(0, len(self.chunks)))
            res.add(tmp)
        return res

    def save_config_file(self, default_map=None, section_name=None):
        if default_map is None:
            default_map = {}
        if section_name is None:
            section_name = self.out_file
        config = configparser.ConfigParser()
        config[section_name] = {'algorithm': 'Online', 'error_correction': self.error_correction.__code__,
                                'insert_header': self.insert_header,
                                'savenumberofchunks': self.save_number_of_chunks_in_packet,
                                'mode_1_bmp': self.mode_1_bmp, 'upper_bound': self.upper_bound,
                                'number_of_chunks': self.number_of_chunks, 'config_str': self.getConfigStr(),
                                'id_len_format': self.id_len_format,
                                'number_of_chunks_len_format': self.number_of_chunks_len_format,
                                'packet_len_format': self.packet_len_format, 'crc_len_format': self.crc_len_format,
                                'master_seed': 0, 'distribution': self.dist.get_config_string(),
                                'rules': [rule for rule in self.rules.active_rules],
                                'chunk_size': self.chunk_size, 'dropped_packets': self.ruleDrop,
                                'created_packets': len(self.encodedPackets)}
        for key, val in default_map.items():
            config[section_name][str(key)] = str(val)
        config_file_name = "{}_{}.ini".format(self.file, datetime.datetime.now().ctime().replace(" ", "_").replace(":",
                                                                                                                   "_"))
        with open(config_file_name, "w") as config_file:
            config.write(config_file)
        return config_file_name

    def getConfigStr(self, out_file=""):
        res = "USE_HEADER_CHUNK: " + str(self.insert_header) + ", NUMBER_OF_CHUNKS: " + str(self.number_of_chunks) + \
              " NUMBER_OF_CHUNKS_LEN_FORMAT: " + self.number_of_chunks_len_format + \
              " ID_LEN_FORMAT: " + self.id_len_format + " ERROR_CORRECTION: " + self.error_correction.__code__.co_name + \
              " CRC_LEN_FORMAT(Optional): " + self.crc_len_format + " FILE: " + self.file + " OUT_FILE: " + out_file + \
              " Distribution: " + self.dist.get_config_string()
        return res


def main(file, number_of_chunks: int = 0, chunk_size: int = 0, error_correction: typing.Callable = nocode,
         as_dna: bool = False, insert_header: bool = False):
    if chunk_size != 0:
        number_of_chunks = Encoder.get_number_of_chunks_for_file_with_chunk_size(file, chunk_size)
    if as_dna:
        rules = FastDNARules()
    else:
        rules = None
    dist = ErlichZielinskiRobustSolitonDistribution(number_of_chunks, seed=2)
    encoder = LTEncoder(file, number_of_chunks, dist, insert_header=insert_header, rules=rules,
                        error_correction=error_correction, number_of_chunks_len_format="H", id_len_format="I",
                        used_packets_len_format="H", save_number_of_chunks_in_packet=False,
                        implicit_mode=False)
    encoder.encode_to_packets()
    print("Number of Chunks=%s" % encoder.number_of_chunks)
    encoder.save_packets(split_to_multiple_files=True, save_as_dna=as_dna)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", metavar="file", type=str, help="the file to Encode")
    parser.add_argument("--error_correction", metavar="error_correction", required=False, type=str, default="nocode",
                        help="Error Correction Method to use; possible values: \
                        nocode, crc, reedsolomon (default=nocode)")
    parser.add_argument("--repair_symbols", metavar="repair_symbols", required=False, type=int, default=2,
                        help="number of repair symbols for ReedSolomon (default=2)")
    parser.add_argument("--chunk_size", metavar="chunk_size", required=False, type=int, default=0,
                        help="size of chunks to split the file into")
    parser.add_argument("--number_of_chunks", metavar="number_of_chunks", required=False, type=int, default=0,
                        help="number of chunks to split the file into,"
                             "only used if no chunk_size is given")
    parser.add_argument("--insert_header", metavar="insert_header", required=False, type=bool, default=False)
    parser.add_argument("--as_dna", help="convert packets to dna and use dna rules", action="store_true",
                        required=False)
    args = parser.parse_args()
    filename = args.filename
    e_correction_str = args.error_correction
    _repair_symbols = args.repair_symbols
    _chunk_size = args.chunk_size
    _number_of_chunks = args.number_of_chunks
    _insert_header = args.insert_header
    _as_dna = args.as_dna
    if _chunk_size == _number_of_chunks == 0:
        print("Please set either a chunk_size or a number_of_chunks")
        exit()
    if e_correction_str == "nocode":
        e_correction = nocode
    elif e_correction_str == "crc":
        e_correction = crc32
    elif e_correction_str == "reedsolomon":
        if _repair_symbols != 2:
            e_correction = lambda x: reed_solomon_encode(x, _repair_symbols)
        else:
            e_correction = reed_solomon_encode
    else:
        print("Selected Error Correction not supported, choose: 'nocode', 'crc' or 'reedsolomon'")
        e_correction = None
        exit()
    filename = args.filename
    print("File to encode: " + str(filename))
    main(filename, _number_of_chunks, _chunk_size, e_correction, _as_dna, _insert_header)
    print("File encoded.")
