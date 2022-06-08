#!/usr/bin/python
# -*- coding: latin-1 -*-
import argparse
import configparser
import datetime
import math
import struct
import time
import numpy
import os
import glob
import typing
from math import ceil

from norec4dna.distributions.Distribution import Distribution
from norec4dna.helper import should_drop_packet, listXOR
from norec4dna.rules.FastDNARules import FastDNARules
from norec4dna.Decoder import Decoder
from norec4dna.ErrorCorrection import get_error_correction_encode, nocode
from norec4dna.distributions.OnlineDistribution import OnlineDistribution
from norec4dna.OnlinePacket import OnlinePacket
from norec4dna.OnlineAuxPacket import OnlineAuxPacket
from norec4dna.Encoder import Encoder


class OnlineEncoder(Encoder):
    def __init__(self, file: str, number_of_chunks: int, distribution: Distribution, epsilon: float, quality: int,
                 insert_header: bool = True, pseudo_decoder: typing.Optional[Decoder] = None, chunk_size: int = 0,
                 rules=None, from_overhead: bool = True, packet_len_format: str = "I", crc_len_format: str = "L",
                 number_of_chunks_len_format: str = "I", quality_len_format: str = "I", epsilon_len_format: str = "f",
                 check_block_number_len_format: str = "I",
                 error_correction: typing.Callable[[typing.Any], typing.Any] = nocode,
                 save_number_of_chunks_in_packet=True, drop_upper_bound=0.0):
        super().__init__(file, number_of_chunks, distribution, insert_header, pseudo_decoder, chunk_size)
        assert (number_of_chunks >= distribution.get_size()), "Epsilon too small for desired number_of_chunks"
        self.out_file: typing.Optional[str] = None
        self.distribution: typing.Optional[Distribution] = distribution
        self.insert_header: bool = insert_header
        self.chunk_size: int = chunk_size
        self.rules = rules
        self.file: str = file
        if self.chunk_size == 0:
            self.number_of_chunks = number_of_chunks
        else:
            self.set_no_chunks_from_chunk_size()
        self.rng: numpy.random = numpy.random
        self.rng.seed(self.number_of_chunks)
        self.chunks: typing.List[bytes] = []
        self.auxBlockNumbers: typing.Dict[int, typing.Set[int]] = dict()
        self.auxBlocks: typing.Dict[int, OnlineAuxPacket] = dict()
        self.encodedPackets: typing.Set[OnlinePacket] = set()
        self.setOfEncodedPackets: typing.Set[int] = set()
        self.pseudo_decoder = pseudo_decoder
        self.epsilon: float = epsilon
        self.quality: int = quality
        self.overhead_limit: float = 0.20
        self.fromOverhead: bool = from_overhead
        self.debug: bool = False
        self.error_correction: typing.Callable[[typing.Any], typing.Any] = error_correction
        # Struct-Strings:
        self.packet_len_format: str = packet_len_format
        self.crc_len_format: str = crc_len_format
        self.number_of_chunks_len_format: str = number_of_chunks_len_format
        self.quality_len_format: str = quality_len_format
        self.epsilon_len_format: str = epsilon_len_format
        self.check_block_number_len_format: str = check_block_number_len_format
        self.save_number_of_chunks_in_packet: bool = save_number_of_chunks_in_packet
        self.upper_bound: float = drop_upper_bound
        self.ruleDrop: int = 0

    def set_overhead_limit(self, n: float):
        self.overhead_limit = n

    def encode_file(self, split_to_multiple_files: bool = False, save_as_dna: bool = True):
        self.encode_to_packets()
        self.save_packets(split_to_multiple_files, save_as_dna=save_as_dna)

    def encode_to_packets(self) -> float:
        def _encode(pseudo):
            pack = self.create_new_packet()
            if self.rules is not None:
                while should_drop_packet(self.rules, pack):
                    pack = self.create_new_packet()
                    self.ruleDrop += 1
            self.encodedPackets.add(pack)
            if pseudo:
                self.pseudo_decoder.input_new_packet(pack)
            self.update_progress_bar()

        self.prepare()
        start = time.time()
        print("Number of Chunks: " + str(self.number_of_chunks))
        limit = (self.number_of_chunks if self.fromOverhead else self.getEstimatedDecodeBlocksNeeded()) * (
                1.0 + self.overhead_limit)
        if self.pseudo_decoder is not None:
            while (not self.pseudo_decoder.is_decoded()) or len(self.encodedPackets) < limit:
                # This process continues until the receiver signals that the
                # message has been received and successfully decoded.
                _encode(True)
        else:
            while len(self.encodedPackets) < limit:
                # This process would continue until the receiver signals that
                # the message has been received and successfully decoded.
                _encode(False)
        return start

    def prepare(self):
        self.ruleDrop: int = 0
        file_size: int = self.get_file_size(self.file)
        if self.insert_header:
            self.chunk_size = ceil(1.0 * file_size / (self.number_of_chunks - 1))
            self.chunks = self.create_chunks(self.chunk_size)
            # First Chunk is a Header
            self.chunks.insert(0, self.encode_header_info())
            self.number_of_chunks += 1  # since we updated number_of_chunks during self.create_chunks...
        else:
            self.chunk_size = ceil(1.0 * file_size / self.number_of_chunks)
            self.chunks = self.create_chunks(self.chunk_size)
        self.fill_last_chunk()
        self.createAuxBlocks()

    def createAuxBlocks(self):
        """ Fuer jeden Chunk eine Anzahl an Aux-Bloecken auswaehlen, in welchen der jeweilige Chunk eingefuegt wird """
        if self.debug:
            print("Using " + str(self.getNumberOfAuxBlocks()) + " Aux-Blocks")
        self.rng.seed(self.number_of_chunks)
        for i in range(0, self.getNumberOfAuxBlocks()):
            self.auxBlockNumbers[i] = set()
        for chunk_num in range(0, self.number_of_chunks):
            # Insert this Chunk into quality different Aux-Packets
            for i in range(0, self.quality):
                # uniform choose a number of aux blocks
                aux_no = self.rng.randint(0, self.getNumberOfAuxBlocks())
                self.auxBlockNumbers[aux_no].add(chunk_num)

        # XOR all Chunks into the corresponding AUX-Block
        for aux_number in self.auxBlockNumbers.keys():
            self.auxBlocks[aux_number] = OnlineAuxPacket(
                listXOR([self.chunks[i] for i in self.auxBlockNumbers[aux_number]]),
                self.auxBlockNumbers[aux_number], aux_number=aux_number, )

    def number_of_packets_encoded_already(self) -> int:
        return len(self.setOfEncodedPackets)

    def save_packets(self, split_to_multiple_files: bool, out_file: typing.Optional[str] = None,
                     save_as_dna: bool = False):
        fileending = ".ONLINE" + ("_DNA" if save_as_dna else "")
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
                fulldir, filename = os.path.split(os.path.realpath(self.file))
                filename = "ONLINE_" + filename
                out_file = os.path.join(fulldir, filename)
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
                with open(out_file + "/" + str(i) + fileending, "wb" if not save_as_dna else "w", ) as f:
                    f.write(packet.get_dna_struct(split_to_multiple_files) if save_as_dna else packet.get_struct(
                        split_to_multiple_files))
                i += 1

    def create_new_packet(self, seed=None) -> OnlinePacket:
        """ Creates a new CheckBlock """
        if seed is None:
            check_block_id = self.generate_new_checkblock_id()
        else:
            check_block_id = seed
        self.distribution.set_seed(check_block_id)
        degree = self.distribution.getNumber()
        packet_numbers: typing.Set[int] = self.choose_packet_numbers(degree, check_block_id)
        packets = []
        for i in packet_numbers:
            if i < len(self.chunks):
                # Add Chunk to list
                packets.append(self.chunks[i])
            else:
                # Add AUX-Block to list
                packets.append(self.auxBlocks[i - len(self.chunks)].get_data())
        self.setOfEncodedPackets |= set(packet_numbers)
        return OnlinePacket(listXOR(packets), self.number_of_chunks, self.quality, self.epsilon, check_block_id,
                            packet_numbers, dist=self.distribution, read_only=False,
                            error_correction=self.error_correction,
                            crc_len_format=self.crc_len_format,
                            number_of_chunks_len_format=self.number_of_chunks_len_format,
                            quality_len_format=self.quality_len_format, epsilon_len_format=self.epsilon_len_format,
                            check_block_number_len_format=self.check_block_number_len_format,
                            save_number_of_chunks_in_packet=self.save_number_of_chunks_in_packet)

    def create_and_add_new_packet(self) -> OnlinePacket:  # , error_correction=nocode
        packet = self.create_new_packet()
        self.encodedPackets.add(packet)
        self.update_progress_bar()
        return packet

    def generate_new_checkblock_id(self) -> int:
        random_generator = numpy.random.RandomState()
        max_num = Encoder.calc_max_size(struct.calcsize("<" + self.check_block_number_len_format))
        return random_generator.randint(0, max_num, dtype=numpy.uint32)

    def choose_packet_numbers(self, degree, seed) -> typing.Set[int]:
        assert degree <= len(self.chunks) + len(self.auxBlocks), (
                str(degree) + ">" + str(len(self.chunks) + len(self.auxBlocks)))
        self.rng.seed(seed)
        res = set()
        for _ in range(0, degree):
            tmp = self.rng.choice(range(0, len(self.chunks) + len(self.auxBlocks)))
            while tmp in res:
                tmp = self.rng.choice(range(0, len(self.chunks) + len(self.auxBlocks)))
            res.add(tmp)
        return res

    def getNumberOfAuxBlocks(self) -> int:
        return int(ceil(0.55 * self.quality * self.epsilon * self.number_of_chunks))

    def getEstimatedDecodeBlocksNeeded(self) -> int:
        """ Optimal Lower-Bound needed for Decoding """
        return ceil((1 + self.epsilon) * (self.number_of_chunks + self.getNumberOfAuxBlocks()))

    def save_config_file(self, default_map=None, section_name=None):
        if default_map is None:
            default_map = {}
        if section_name is None:
            section_name = self.out_file
        config = configparser.ConfigParser()
        config[section_name] = {'algorithm': 'Online', 'error_correction': self.error_correction.__code__,
                                'insert_header': self.insert_header,
                                'savenumberofchunks': self.save_number_of_chunks_in_packet,
                                'upper_bound': self.upper_bound, 'number_of_chunks': self.number_of_chunks,
                                'config_str': self.getConfigStr(), 'id_len_format': self.check_block_number_len_format,
                                'number_of_chunks_len_format': self.number_of_chunks_len_format,
                                'packet_len_format': self.packet_len_format, 'crc_len_format': self.crc_len_format,
                                'quality_len_format': self.quality_len_format,
                                'epsilon_len_format': self.epsilon_len_format,
                                'master_seed': 0, 'distribution': self.distribution.get_config_string(),
                                'rules': [rule for rule in self.rules.active_rules],
                                'chunk_size': self.chunk_size, 'dropped_packets': self.ruleDrop,
                                'created_packets': len(self.encodedPackets)}
        for key, val in default_map.items():
            config[section_name][str(key)] = str(val)
        config_file_name = "{}_{}.ini".format(self.file,
                                              datetime.datetime.now().ctime().replace(" ", "_").replace(":", "_"))
        with open(config_file_name, "w") as config_file:
            config.write(config_file)
        return config_file_name

    def getConfigStr(self, out_file=""):
        res = "USE_HEADER_CHUNK: " + str(self.insert_header) + ", NUMBER_OF_CHUNKS: " + str(self.number_of_chunks) + \
              " NUMBER_OF_CHUNKS_LEN_FORMAT: " + self.number_of_chunks_len_format + \
              " ID_LEN_FORMAT: " + self.check_block_number_len_format + " PACKET_LEN_FORMAT: " + self.packet_len_format + \
              " QUALITY_LEN_FORMAT: " + self.quality_len_format + " EPSILON_LEN_FORMAT: " + self.quality_len_format + \
              " ERROR_CORRECTION: " + self.error_correction.__code__.co_name + \
              " CRC_LEN_FORMAT(Optional): " + self.crc_len_format + " FILE: " + self.file + " OUT_FILE: " + out_file + \
              " Distribution: " + self.distribution.get_config_string()
        return res


def roundup(x) -> int:
    return int(math.ceil(x / 20.0)) * 20


def main(file: str, error_correction: typing.Callable[[typing.Any], typing.Any], asdna: bool = True,
         epsilon: float = 0.06, insert_header: bool = False):
    dist = OnlineDistribution(epsilon)
    number_of_chunks = dist.get_size()
    quality = 7
    if asdna:
        rules = FastDNARules()
    else:
        rules = None
    encoder = OnlineEncoder(
        file, number_of_chunks, dist, epsilon, quality, error_correction=error_correction, quality_len_format="B",
        insert_header=insert_header, check_block_number_len_format="H", number_of_chunks_len_format="H", rules=rules,
        save_number_of_chunks_in_packet=False)
    encoder.set_overhead_limit(1.70)
    encoder.encode_file(split_to_multiple_files=True, save_as_dna=asdna)
    encoder.save_packets(True, save_as_dna=asdna)


if __name__ == "__main__":
    # try:
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", metavar="file", type=str, help="the file to Encode")
    parser.add_argument("--error_correction", metavar="error_correction", required=False, type=str, default="nocode",
                        help="Error Correction Method to use; possible values: \
                        nocode, crc, reedsolomon (default=nocode)")
    parser.add_argument("--repair_symbols", metavar="repair_symbols", required=False, type=int, default=2,
                        help="number of repair symbols for ReedSolomon (default=2)")
    parser.add_argument("--insert_header", metavar="insert_header", required=False, type=bool, default=False)
    parser.add_argument("--as_dna", help="convert packets to dna and use dna rules", action="store_true",
                        required=False)
    parser.add_argument("--epsilon", metavar="epsilon", required=False, type=float, default=0.06,
                        help="epsilon to use for the distribution")
    args = parser.parse_args()
    _file = args.filename
    _repair_symbols = args.repair_symbols
    _insert_header = args.insert_header
    _as_dna = args.as_dna
    _error_correction = get_error_correction_encode(args.error_correction, _repair_symbols)
    print("File to encode: " + str(_file))
    main(_file, _error_correction, _as_dna, insert_header=_insert_header)
    print("File encoded.")
