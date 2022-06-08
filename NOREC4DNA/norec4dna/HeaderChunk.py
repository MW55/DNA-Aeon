import struct
import numpy
import typing

from norec4dna.Packet import Packet


class HeaderChunk:
    def __init__(self, packet: Packet, last_chunk_len_format: str = "I"):
        assert packet.get_used_packets().issubset({0}), "only first packet can be HeaderPacket"
        if isinstance(packet.data, numpy.ndarray):
            self.data: bytes = packet.data.tobytes()
        else:
            self.data: bytes = packet.data
        self.last_chunk_len_format: str = last_chunk_len_format
        self.last_chunk_length, self.file_name = self.decode_header_info()

    def get_last_chunk_length(self) -> int:
        return self.last_chunk_length

    def get_file_name(self) -> typing.Union[str, bytes]:
        return self.file_name

    def decode_header_info(self) -> typing.Tuple[int, typing.Union[bytes, str]]:
        # Size of last Chunk
        # Filename
        # PAD-Bytes
        last_chunk_struct_len: int = struct.calcsize("<" + self.last_chunk_len_format)
        last_chunk_length: int = struct.unpack("<" + self.last_chunk_len_format,
                                               bytes(self.data[0:last_chunk_struct_len]))[0]
        data: bytes = self.data
        end_of_file_name: int = data.find(0x00, last_chunk_struct_len + 1)
        file_name: typing.Union[str, bytes] = \
            struct.unpack("<" + str(len(data[last_chunk_struct_len:end_of_file_name])) + "s",
                          data[last_chunk_struct_len:end_of_file_name])[0]
        return last_chunk_length, file_name

    def __str__(self) -> str:
        return "< last_chunk_length: " + str(self.last_chunk_length) + " , file_name: " + str(self.file_name) + " >"

    def __repr__(self):
        return self.__str__()
