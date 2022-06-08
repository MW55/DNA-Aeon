import math
import random
import secrets
import struct
from bitarray import bitarray
from Cryptodome.Protocol.SecretSharing import Shamir

from norec4dna.RU10Encoder import RU10Encoder
from norec4dna.Encoder import Encoder
from norec4dna.RU10Decoder import RU10Decoder
from norec4dna.ErrorCorrection import nocode, reed_solomon_encode, crc32
from norec4dna.distributions.RaptorDistribution import RaptorDistribution
from norec4dna.helper.RU10Helper import intermediate_symbols, from_true_false_list


class Multiplex:
    def __init__(self, file, chunk_size, no_channel, no_secure_channel, error_correction=nocode, header=False,
                 id_len_format='H'):
        self.file = file
        self.chunk_size = chunk_size
        self.header = header
        self.no_chunks = Encoder.get_number_of_chunks_for_file_with_chunk_size(self.file, self.chunk_size,
                                                                               insert_header=self.header)
        self.error_correction = error_correction
        self.id_len_format = id_len_format
        if no_channel < 3:
            print("You need to use at least 3 channel")
            return
        self.no_channel = no_channel
        if no_channel < no_secure_channel:
            print("Number of secure channels needs to be <= number of channels")
            self.no_secure_channel = self.no_channel
            print("New number of secure channels is " + str(self.no_channel))
        else:
            self.no_secure_channel = no_secure_channel
        self.no_public_channel = self.no_channel - self.no_secure_channel
        print("Number of public channels is " + str(self.no_public_channel))
        self.decoder = self.init_decoder()
        self.channel, self.used_methods = self.init_channel()
        print("Created decoder and " + str(no_channel) + " channel.")
        self.packets = list()

    def init_decoder(self):
        """
        Initializes the decoder to check if the receiver is potentially able to decode the file if he got all packets
        created for the channels.
        :return:
        """
        decoder = RU10Decoder.pseudo_decoder(self.no_chunks, False)
        if decoder.distribution is None:
            decoder.distribution = RaptorDistribution(self.no_chunks)
            decoder.number_of_chunks = self.no_chunks
            _, decoder.s, decoder.h = intermediate_symbols(self.no_chunks, decoder.distribution)
            decoder.createAuxBlocks()
        return decoder

    def get_packets_for_channel(self, no_packets, channel_number):
        """
        Uses the channel method to create a list of packets, adds them to the list of all created packets and returns
        them.
        :param no_packets: Number of packets to create
        :param channel_number: Channel to create the packets for
        :return:
        """
        try:
            packets = self.channel[channel_number].create_packets(no_packets)
            self.packets.extend(packets)
            for pack in packets:
                self.decoder.input_new_packet(pack)
            return packets
        except:
            print("Channel not available.")

    def file_potentially_decodable(self):
        """
        Checks if the receiver is able to decode the file if he got all packets created.
        :return:
        """
        if self.decoder.GEPP is not None:
            return self.decoder.is_decoded()
        else:
            return False

    def get_channel(self):
        """
        Returns the channel objects stored in the Multiplexer. Not needed anymore.
        :return:
        """
        return self.channel

    def init_channel(self):
        """
        Initializes the channels used to multiplex the file. no_channel channels are created with no_secure_channel of
        them being "secure" channels which means that they are allowed to generate packets with grades greater than 1
        that can't ensure that no single chunk is decodable with these channels packets. The other channels are "public"
        channels and are only allowed to send packets with even grades which means, that no single chunk can be decoded
        with the packets of these channels.
        :return:
        """
        channel = list()
        methods = self.get_channel_methods()
        for x in range(0, self.no_secure_channel):
            channel.append(MultiplexChannel(self.file, methods[x][0], self.chunk_size, methods[x][1], secure=True,
                                            id_len_format=self.id_len_format, error_correction=self.error_correction))
        for y in range(0, self.no_public_channel):
            channel.append(MultiplexChannel(self.file, methods[y + self.no_secure_channel][0], self.chunk_size,
                                            methods[y + self.no_secure_channel][1], secure=False,
                                            id_len_format=self.id_len_format, error_correction=self.error_correction))
        return channel, methods

    def get_channel_methods(self):
        """
        Gets the methods used to create packets from specific chunks for every channel. The first channel always only
        uses chunks with even numbers, the second chunks with odd numbers and to be able to decode the file with
        packets from all channels, the third channel always uses chunks from the first 30 chunk wide window. This makes
        sure that there is an intersection in the chunks used for the channels.
        :return:
        """
        methods = list()
        methods.append(("even", 0))
        methods.append(("odd", 0))
        methods.append(("window_30", 0))
        max_window = math.floor(self.no_chunks / 20)
        for x in range(1, self.no_channel - 2):
            if x > max_window:
                x = x % max_window
            methods.append(("window_30", x))
        return methods

    def add_channel(self, secure=False):
        """
        Adds a new channel to the Multiplexer. The channel can be secure or public. The function checks if any method is
        not used for a channel. If all are used at least once it chooses a random window for the channel.
        :param secure:
        :return:
        """
        window = 0
        if ("even", 0) not in self.used_methods:
            method = "even"
        elif ("odd", 0) not in self.used_methods:
            method = "odd"
        else:
            method = "window_30"
            max_window = math.floor(self.no_chunks / 20)
            for x in range(0, max_window):
                if (method, x) not in self.used_methods:
                    window = x
                    break
                else:
                    window = random.randint(0, max_window)
        channel = MultiplexChannel(self.file, method, self.chunk_size, window=window, secure=secure,
                                   id_len_format=self.id_len_format, error_correction=self.error_correction,
                                   header=self.header)
        self.channel.append(channel)
        self.used_methods.append((method, window))
        return True

    def remove_channel(self, channel_no):
        """
        Deletes a channel and checks if the receiver could be able to decode the file with packets from the remaining
        channels (one secure channel is needed). Maybe we should check if all chunks are used on the remaining channels.
        :param channel_no:
        :return:
        """
        secure_channel = 0
        for chan in self.channel:
            if chan.secure is True:
                secure_channel += 1
        if self.channel[channel_no].secure is True and secure_channel == 1:
            print("Channel is not deletable. The receiver would not be able to decode the file anymore if you delete "
                  "the last secure channel.")
            return False
        else:
            try:
                del self.used_methods[channel_no]
                del self.channel[channel_no]
                print("Channel " + str(channel_no) + " deleted.")
                return True
            except:
                print("Channel " + str(channel_no) + " not available.")
                return False

    def all_chunks_used(self):
        """
        Checks if all chunks are used by channels. If not, the decoder will never be able to decode the file with the
        packets from the channels.
        :return:
        """
        chunk_list = [x for x in range(0, self.no_chunks)]
        for method in self.used_methods:
            if method[0] == "even":
                chunk_list = [x for x in chunk_list if x % 2 != 0]
            elif method[0] == "odd":
                chunk_list = [x for x in chunk_list if x % 2 == 0]
            elif method[0] == "window":
                window_size = 30
                window = method[1]
                start = window * (window_size - 10)
                chunks = [y for y in range(start, start + window_size)]
                chunk_list = [x for x in chunk_list if x not in chunks]
            if len(chunk_list) == 0:
                return True
        return False

    def get_shamir_shares(self, min_shares):
        """
        16 Byte secret:
        number of chunks 4 byte + header and error_correction encoding 3 bit + 5 0 bits = 1 byte + 9 random bytes
        :param min_shares: Minimum shares needed to reconstruct the secret
        :return:
        """
        if self.no_channel < min_shares:
            print("Please decrease the minimum since there are less channel.")
            return None
        noc = struct.pack('I', self.no_chunks)
        arr = bitarray()
        arr.append(self.header)
        if self.error_correction is nocode:
            arr.extend([False, False])
        elif self.error_correction is crc32:
            arr.extend([False, True])
        elif self.error_correction is reed_solomon_encode:
            arr.extend([True, False])
        while len(arr) < 8:
            arr.append(random.getrandbits(1))
        fill = secrets.token_bytes(11)
        secret = noc + bytes(arr) + fill
        return Shamir.split(min_shares, self.no_channel, secret)

    @staticmethod
    def combine_shamir_shares(shares):
        """
        Combines the shamir shares to get the no_chunks, header and error_correction. Needs at least the minimum shares
        set in @get_shamir_shares.
        :param shares:
        :return:
        """
        byte_str = Shamir.combine(shares)
        no_chunks = struct.unpack('I', byte_str[:4])[0]
        bool_byte = byte_str[4:5]
        header = bool((bool_byte[0] >> 7) & 1)
        err_cor_1 = bool((bool_byte[0] >> 6) & 1)
        err_cor_2 = bool((bool_byte[0] >> 5) & 1)
        if not err_cor_1 and not err_cor_2:
            error_correction = nocode
        elif err_cor_1 is False and err_cor_2 is True:
            error_correction = crc32
        elif err_cor_1 is True and err_cor_2 is False:
            error_correction = reed_solomon_encode
        else:
            raise RuntimeError("illegal state")
        return no_chunks, header, error_correction


class MultiplexChannel:
    def __init__(self, file, method, chunk_size, window=0, secure=False, id_len_format='H', error_correction=nocode,
                 header=False):
        self.file = file
        self.method = method
        self.window = window
        self.chunk_size = chunk_size
        self.secure = secure
        self.id_len_format = id_len_format
        self.error_correction = error_correction
        self.header = header
        self.encoder = self.init_encoder()
        self.decoder = self.init_decoder()
        self.packets = list()

    def init_encoder(self):
        """
        Initializes the encoder to generate the packets for the channel.
        :return: Encoder instance
        """
        no_chunks = Encoder.get_number_of_chunks_for_file_with_chunk_size(self.file, self.chunk_size,
                                                                          insert_header=self.header)
        dist = RaptorDistribution(no_chunks)
        rules = None
        enc = RU10Encoder(self.file, no_chunks, dist, chunk_size=self.chunk_size, insert_header=self.header,
                          rules=rules,
                          error_correction=self.error_correction, id_len_format=self.id_len_format,
                          number_of_chunks_len_format="B",
                          save_number_of_chunks_in_packet=False, mode_1_bmp=False)
        enc.prepare()
        return enc

    def init_decoder(self):
        dec = RU10Decoder.pseudo_decoder(self.encoder.number_of_chunks, False)
        if dec.distribution is None:
            dec.distribution = RaptorDistribution(self.encoder.number_of_chunks)
            dec.number_of_chunks = self.encoder.number_of_chunks
            _, dec.s, dec.h = intermediate_symbols(self.encoder.number_of_chunks, dec.distribution)
            dec.createAuxBlocks()
        return dec

    def create_packets(self, no_packets=1):
        """
        Creates no_packets packets based on the method given for the channel and checks for the packets grade. Secure
        channels are allowed to use packets with grades greater than 1 and not-secure channels are only allowed to use
        packets with even grades.
        :param no_packets:
        :return:
        """
        packets = []
        discarded_packets = 0
        while len(packets) < no_packets:
            packet = self.encoder.create_new_packet_from_chunks(method=self.method, window=self.window)
            packet_chunks = from_true_false_list(self.decoder.removeAndXorAuxPackets(packet))
            if len(packet_chunks) > 1 and self.secure or (not self.secure and len(packet_chunks) % 2 == 0):
                packets.append(packet)
                self.packets.append(packet)
                self.decoder.input_new_packet(packet)
            else:
                discarded_packets += 1
        print("Generated " + str(no_packets + discarded_packets) + " packets and discarded " + str(discarded_packets))
        return packets

    def is_decodable(self):
        """
        Checks if the packets created for the channel allow decoding the file.
        :return:
        """
        if self.decoder.GEPP is not None:
            return self.decoder.is_decoded()


if __name__ == '__main__':
    file = '.INFILES/Dorn'
    cmp_file = 'tests/cmp_dorn'


    def do_test_decodable(error_correction=nocode, header=True, no_channel=5):
        mlt = Multiplex(file, 50, no_channel, 1, error_correction=error_correction, header=header)
        for ind, ch in enumerate(mlt.channel):
            mlt.get_packets_for_channel(200, ind)
        assert mlt.file_potentially_decodable() is True
        for ch in mlt.channel:
            assert ch.is_decodable() is False


    do_test_decodable()
