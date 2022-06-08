#!/usr/bin/python
# -*- coding: latin-1 -*-
import argparse
import configparser
import datetime
import struct
import os
import glob
import typing
import numpy as np
from math import ceil, floor

from norec4dna.helper.RU10Helper import int31, choose_packet_numbers, intermediate_symbols
from norec4dna.distributions.Distribution import Distribution
from norec4dna.distributions.RaptorDistribution import RaptorDistribution
from norec4dna.helper import should_drop_packet, listXOR, calc_crc, buildGraySequence, bitSet
from norec4dna.rules.FastDNARules import FastDNARules
from norec4dna.ErrorCorrection import get_error_correction_encode, nocode
from norec4dna.Encoder import Encoder
from norec4dna.RU10Packet import RU10Packet


class RU10Encoder(Encoder):
    def __init__(self, file, number_of_chunks, distribution: Distribution, insert_header=True, pseudo_decoder=None,
                 chunk_size=0, rules=None, error_correction=nocode, packet_len_format="I", crc_len_format="L",
                 number_of_chunks_len_format="L", id_len_format="L", save_number_of_chunks_in_packet=True,
                 mode_1_bmp=False, prepend="", append="", drop_upper_bound=1.0, keep_all_packets=False):
        super().__init__(file, number_of_chunks, distribution, insert_header, pseudo_decoder,
                         chunk_size, mode_1_bmp)
        self.success_packets = 0
        self.out_file = None
        self.esi: int = 0
        self.debug: bool = False
        self.file: str = file
        self.dist: Distribution = distribution
        self.chunk_size: int = chunk_size
        self.insert_header: bool = insert_header
        self.rules = rules
        self.mode_1_bmp: bool = mode_1_bmp
        self.upper_bound = drop_upper_bound
        if self.chunk_size == 0:
            self.number_of_chunks: int = number_of_chunks
        else:
            self.set_no_chunks_from_chunk_size()
            # we have to update the Dist noOfChunks..
            self.dist.S = self.number_of_chunks
            self.dist.rng.seed(self.number_of_chunks)
        self.chunks: typing.List[bytes] = []
        self.pseudo_decoder = pseudo_decoder
        self.encodedPackets: typing.Set[RU10Packet] = set()
        self.overhead_limit: float = 2.50
        self.error_correction: typing.Callable[[typing.Any], typing.Any] = error_correction
        self.packet_len_format: str = packet_len_format
        self.crc_len_format: str = crc_len_format
        self.number_of_chunks_len_format: str = number_of_chunks_len_format
        self.id_len_format: str = id_len_format
        self.save_number_of_chunks_in_packet: bool = save_number_of_chunks_in_packet
        # if we use None, we will mutate the RNG each call,
        # otherwise we will use the same RNG to generate different values
        self.random_state: np.random.RandomState = np.random.RandomState()
        if self.random_state is not None:
            self.__masterseed = self.random_state.get_state()[1][0]
            # print("Master-Seed used: " + str(self.__masterseed))
        self.prepend = prepend
        self.append = append
        self.ruleDrop: int = 0

    """ 
    Creates Chunks from the given File, fills last Chunk with padding, 
    adds a header packet and saves them to self.chunks
    """

    def prepare(self, systematic: bool = False):
        """
        Creates chunks for the given file, inserts a header as first chunk and fills the last chunk with padding.
        The chunks are saved to self.chunks.
        :param systematic:
        :return:
        """
        self.ruleDrop: int = 0
        file_size: int = self.get_file_size(self.file)
        if self.insert_header:
            self.chunk_size = ceil(1.0 * file_size / (self.number_of_chunks - 1))
            self.chunks = self.create_chunks(self.chunk_size)
            # First Chunk is a Header
            self.chunks.insert(0, self.encode_header_info())
            self.number_of_chunks += 1  # since the update for number_of_chunks happend inside create_chunks.
        else:
            self.chunk_size = ceil(1.0 * file_size / self.number_of_chunks)
            self.chunks = self.create_chunks(self.chunk_size)
        self.fill_last_chunk()
        # create Intermediate Blocks and save them as self.chunks
        self.generate_intermediate_blocks()

    def encode_to_packets(self):
        """
        Prepares the chunks with @prepare and calls @doEncode
        :return:
        """
        self.prepare()
        self.do_encode()

    def do_encode(self):
        """
        If a pseudodecoder is used packets will be generated until the pseudodecoder would be able to decode the file.
        Else packets are generated until there are n*overhead successfully generated packets and every chunk is filled.
        :return:
        """

        def _encode(pseudo):
            new_pack = self.create_new_packet()
            if self.rules is not None:
                while should_drop_packet(self.rules, new_pack, self.upper_bound):
                    del new_pack
                    new_pack = self.create_new_packet()
                    self.ruleDrop += 1
            if pseudo:
                self.pseudo_decoder.input_new_packet(new_pack)
            self.encodedPackets.add(new_pack)
            self.update_progress_bar()

        if self.pseudo_decoder is not None:
            _encode(pseudo=True)
            while not self.pseudo_decoder.is_decoded():
                # This process continues until the receiver signals that the
                # message has been received and successfully decoded.
                _encode(pseudo=True)
        else:
            while (self.number_of_packets_encoded_already() < (self.number_of_chunks + (
                    self.number_of_chunks * self.overhead_limit)) or self.number_of_packets_encoded_already() < self.number_of_chunks):  # % Aufschlag
                # This process continues until the receiver signals that the
                # message has been received and successfully decoded.
                _encode(pseudo=False)

    def generate_new_id(self, systematic: bool = False) -> int:
        """
        Generates a random ID for a new packet.
        :param systematic:
        :return: Random ID
        """
        if systematic and self.esi < self.number_of_chunks:
            self.esi += 1
            return self.esi - 1
        if self.random_state is None:
            random_generator = np.random.RandomState()
        else:
            random_generator = self.random_state
        # RU10 is only defined for max int31!
        max_num = min(Encoder.calc_max_size(struct.calcsize("<" + self.id_len_format)), int31)
        return random_generator.randint(0, max_num, dtype=np.uint32)

    def create_new_packet(self, systematic: bool = False, seed: typing.Optional[int] = None) -> RU10Packet:
        """
        Creates a new RU10Packet with an ID, the number of participated packets and the XORed payload.
        :param seed:
        :param systematic:
        :return: A new RU10Packet
        """
        if seed is None:
            seed = self.generate_new_id()
        packet_numbers = choose_packet_numbers(self.number_of_chunks, seed, self.dist, systematic=systematic)
        packets = [self.chunks[i] for i in packet_numbers]
        if self.debug:
            print("----")
            print("Id = " + str(seed))
            print(packet_numbers)
            print(calc_crc(listXOR(packets)))
            print("----")
        return RU10Packet(listXOR(packets), packet_numbers, self.number_of_chunks, seed, dist=self.dist,
                          read_only=False, error_correction=self.error_correction,
                          packet_len_format=self.packet_len_format, crc_len_format=self.crc_len_format,
                          number_of_chunks_len_format=self.number_of_chunks_len_format,
                          id_len_format=self.id_len_format,
                          save_number_of_chunks_in_packet=self.save_number_of_chunks_in_packet, prepend=self.prepend,
                          append=self.append)

    def create_new_packet_from_chunks(self, method: str, systematic: bool = False, seed: typing.Optional[int] = None,
                                      window: int = 0) -> typing.Optional[RU10Packet]:
        """
        Creates a new packet only using chunks based on the choosen method. You can use "even" or "odd" which
        basically uses chunks with even or odd numbers. "window_30" and "window_40" are methods, that split the whole
        chunklist into windows with 30 or 40 chunks overlapping by 10 chunks for each window. "window" says what window
        the packet should be generated for, starting with 0 for the window that contains the first 30 or 40 chunks.
        The packets will be added to the encoders packets and can be saved with the save_packets method like normal
        packets. To decode the packets set use_method = True while creating a decoder.
        :param method:
        :param systematic:
        :param seed:
        :param window:
        :return: A new RU10Packet
        """
        if method == "window_30":
            window_size = 30
            start = window * (window_size - 10)
            if start > self.number_of_chunks:
                print("Please select a valid window.")
                return None
            chunk_lst: typing.List[int] = [ch for ch in range(start, start + window_size) if
                                           ch <= self.number_of_chunks]
        elif method == "window_40":
            window_size = 40
            start = window * (window_size - 10)
            if start > self.number_of_chunks:
                print("Please select a valid window.")
                return None
            chunk_lst = [ch for ch in range(start, start + window_size) if ch <= self.number_of_chunks]
        elif method == "even":
            chunk_lst = [ch for ch in range(0, self.number_of_chunks + 1) if ch % 2 == 0]
        elif method == "odd":
            chunk_lst = [ch for ch in range(0, self.number_of_chunks + 1) if ch % 2 != 0]
        else:
            raise RuntimeError("Please select a valid method (even, odd, window (including window_size and start_ind)")
        if seed is None:
            seed = self.generate_new_id()
        packet_numbers = choose_packet_numbers(len(chunk_lst), seed, self.dist, systematic=systematic,
                                               max_l=len(chunk_lst))
        packets = [self.chunks[chunk_lst[i]] for i in packet_numbers]
        packet = RU10Packet(listXOR(packets), [chunk_lst[i] for i in packet_numbers], self.number_of_chunks, seed,
                            read_only=False, error_correction=self.error_correction,
                            packet_len_format=self.packet_len_format, crc_len_format=self.crc_len_format,
                            number_of_chunks_len_format=self.number_of_chunks_len_format,
                            id_len_format=self.id_len_format,
                            save_number_of_chunks_in_packet=self.save_number_of_chunks_in_packet,
                            method=method, window=window, prepend=self.prepend, append=self.append)
        self.encodedPackets.add(packet)
        return packet

    def generate_intermediate_blocks(self) -> typing.List[bytes]:
        """
        Generates intermediate blocks used to generate the complete auxblocks afterwards.
        :return: Self.chunks containing the intermediate blocks
        """
        _, s, h = intermediate_symbols(self.number_of_chunks, self.dist)

        k = self.number_of_chunks
        compositions: typing.List[typing.List[int]] = [[] for _ in range(s)]
        for i in range(0, k):
            a = 1 + (int(floor(np.float64(i) / np.float64(s))) % (s - 1))
            b = int(i % s)
            compositions[b].append(i)
            b = (b + a) % s
            compositions[b].append(i)
            b = (b + a) % s
            compositions[b].append(i)
        for i in range(0, s):
            b = listXOR([self.chunks[x] for x in compositions[i]])
            self.chunks.append(b)
            if self.debug:
                print(str(len(self.chunks) - 1) + " : " + str(compositions[i]))
        hprime = int(ceil(np.float64(h) / 2))
        m = buildGraySequence(k + s, hprime)
        for i in range(0, h):
            hcomposition: typing.List[int] = []
            for j in range(0, k + s):
                if bitSet(np.uint32(m[j]), np.uint32(i)):
                    hcomposition.append(j)
            b = listXOR([self.chunks[x] for x in hcomposition])
            self.chunks.append(b)
            if self.debug:
                print(str(len(self.chunks) - 1) + " : " + str(hcomposition))
        return self.chunks

    def encode_file(self, split_to_multiple_files: bool = False):
        """
        Encodes a filed to packets (@encode_to_packets) and saves the generated packets (@save_packets)
        :param split_to_multiple_files:
        :return:
        """
        self.encode_to_packets()
        self.save_packets(split_to_multiple_files)

    def save_packets(self, split_to_multiple_files: bool, out_file: typing.Optional[str] = None,
                     save_as_dna: bool = False, clear_output: bool = True, seed_is_filename: bool = False):
        """
        Saves the generated packets either to multiple files or to a single one. It's possible to save the packets
        either as DNA or binary.
        :param split_to_multiple_files: True: Saves the packets in multiple files. False: Saves all packets in one file.
        :param out_file: The location of the output file
        :param save_as_dna: True: Saves the information in bases. False: Saves the information binary.
        :param clear_output: Clears the location of the output file.
        :param seed_is_filename: True: Sets the seed as filename.
        :return:
        """
        file_ending = ".RU10" + ("_DNA" if save_as_dna else "")
        if not split_to_multiple_files:
            if out_file is None:
                out_file = self.file + file_ending
            with open(out_file, "wb" if not save_as_dna else "w") as f:
                for packet in self.encodedPackets:
                    f.write(
                        packet.get_dna_struct(split_to_multiple_files)
                        if save_as_dna
                        else packet.get_struct(split_to_multiple_files)
                    )
        else:
            # Folder:
            if out_file is None:
                fulldir, filename = os.path.split(os.path.realpath(self.file))
                filename = "RU10_" + filename
                out_file = os.path.join(fulldir, filename)
                if not out_file.endswith("/"):
                    files = glob.glob(out_file + "/*")
                else:
                    files = glob.glob(out_file + "*")
                if clear_output:
                    for f in files:
                        os.remove(f)
            i = 0
            e_prob = ""
            if not os.path.exists(out_file):
                os.makedirs(out_file)
            for packet in sorted(self.encodedPackets, key=lambda elem: (elem.error_prob, elem.__hash__())):
                if seed_is_filename:
                    i = packet.id
                    e_prob = (str(ceil(packet.error_prob * 100)) + "_") if packet.error_prob is not None else ""
                with open(out_file + "/" + e_prob + str(i) + file_ending, "wb" if not save_as_dna else "w") as f:
                    f.write(
                        packet.get_dna_struct(split_to_multiple_files)
                        if save_as_dna
                        else packet.get_struct(split_to_multiple_files)
                    )
                i += 1
            self.out_file = os.path.relpath(out_file)
            print("Config: " + self.getConfigStr(out_file))

    def getConfigStr(self, out_file=""):
        res = "USE_HEADER_CHUNK: " + str(self.insert_header) + ", NUMBER_OF_CHUNKS: " + str(self.number_of_chunks) + \
              " NUMBER_OF_CHUNKS_LEN_FORMAT: " + self.number_of_chunks_len_format + \
              " ID_LEN_FORMAT: " + self.id_len_format + " ERROR_CORRECTION: " + self.error_correction.__code__.co_name + \
              " CRC_LEN_FORMAT(Optional): " + self.crc_len_format + " FILE: " + self.file + " OUT_FILE: " + out_file + \
              " Distribution: " + self.dist.get_config_string()
        return res

    def save_config_file(self, default_map=None, section_name=None):
        if default_map is None:
            default_map = {}
        if section_name is None:
            section_name = self.out_file
        config = configparser.ConfigParser()
        config[section_name] = {'algorithm': "RU10",
                                'error_correction': self.error_correction.__code__,
                                'insert_header': self.insert_header,
                                'savenumberofchunks': self.save_number_of_chunks_in_packet,
                                'mode_1_bmp': self.mode_1_bmp, 'upper_bound': self.upper_bound,
                                'number_of_chunks': self.number_of_chunks, 'config_str': self.getConfigStr(),
                                'id_len_format': self.id_len_format,
                                'number_of_chunks_len_format': self.number_of_chunks_len_format,
                                'packet_len_format': self.packet_len_format, 'crc_len_format': self.crc_len_format,
                                'master_seed': self.__masterseed, 'distribution': self.dist.get_config_string(),
                                'rules': [rule for rule in self.rules.active_rules] if self.rules is not None else [],
                                'chunk_size': self.chunk_size, 'dropped_packets': self.ruleDrop,
                                'created_packets': len(self.encodedPackets)}
        for key, val in default_map.items():
            config[section_name][str(key)] = str(val)
        config_file_name = "{}_{}.ini".format(self.file, datetime.datetime.now().ctime().replace(" ", "_")).replace(":",
                                                                                                                    "_")
        with open(config_file_name, "w") as config_file:
            config.write(config_file)
        return config_file_name


def main(in_file: str, num_chunks=0, chunk_size=0, as_dna=True,
         err_correction: typing.Callable[[typing.Any], typing.Any] = nocode, insert_header=False,
         save_number_of_chunks_in_packet=False, mode_1_bmp=False):
    if chunk_size != 0:
        num_chunks = Encoder.get_number_of_chunks_for_file_with_chunk_size(in_file, chunk_size)
        dist = RaptorDistribution(num_chunks)
    elif num_chunks != 0:
        dist = RaptorDistribution.RaptorDistribution(num_chunks)
    else:
        print("Aborting. Please set either chunk_size or number_of_chunks!")
        return
    if as_dna:
        rules = FastDNARules()
    else:
        rules = None
    x = RU10Encoder(in_file, num_chunks, dist, chunk_size=chunk_size, insert_header=insert_header, rules=rules,
                    error_correction=err_correction, id_len_format="H", number_of_chunks_len_format="B",
                    save_number_of_chunks_in_packet=save_number_of_chunks_in_packet, mode_1_bmp=mode_1_bmp)
    x.encode_to_packets()
    x.save_packets(True, save_as_dna=as_dna, seed_is_filename=False)
    conf = {'error_correction': e_correction, 'repair_symbols': norepair_symbols, 'asdna': as_dna,
            'number_of_splits': 0, 'find_minimum_mode': True, 'seq_seed': False}
    x.save_config_file(conf, section_name="RU10_" + in_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--as_dna", help="convert packets to dna and use dna rules", action="store_true",
                        required=False)
    parser.add_argument("filename", metavar="file", type=str, help="the file to Encode")
    parser.add_argument("--chunk_size", metavar="chunk_size", required=False, type=int, default=0,
                        help="size of chunks to split the file into")
    parser.add_argument("--number_of_chunks", metavar="number_of_chunks", required=False, type=int, default=0,
                        help="number of chunks to split the file into,"
                             "only used if no chunk_size is given")
    parser.add_argument("--error_correction", metavar="error_correction", required=False, type=str, default="nocode",
                        help="Error Correction Method to use; possible values: \
                        nocode, crc, reedsolomon, dna_reedsolomon (default=nocode)")
    parser.add_argument("--repair_symbols", metavar="repair_symbols", required=False, type=int, default=2,
                        help="number of repair symbols for ReedSolomon (default=2)")
    parser.add_argument("--insert_header", metavar="insert_header", required=False, type=bool, default=False)
    parser.add_argument("--save_number_of_chunks", metavar="save_number_of_chunks", required=False, type=bool,
                        default=False)
    parser.add_argument("--as_mode_1_bmp", required=False, action="store_true")
    overhead = 6.0
    args = parser.parse_args()
    arg_file = args.filename
    arg_chunk_size = args.chunk_size
    arg_number_of_chunks = args.number_of_chunks
    arg_as_dna = args.as_dna
    e_correction = args.error_correction
    norepair_symbols = args.repair_symbols
    arg_insert_header = args.insert_header
    save_number_of_chunks = args.save_number_of_chunks
    arg_mode_1_bmp = args.as_mode_1_bmp
    if arg_chunk_size == arg_number_of_chunks == 0:
        print("Please set either a chunk_size or a number_of_chunks")
        exit()
    arg_error_correction = get_error_correction_encode(e_correction, norepair_symbols)
    print("File to encode: " + str(arg_file))
    main(arg_file, arg_number_of_chunks, arg_chunk_size, arg_as_dna, arg_error_correction, arg_insert_header,
         save_number_of_chunks, arg_mode_1_bmp)
    print("File encoded.")
