import struct
import typing

from norec4dna.ErrorCorrection import nocode
from norec4dna.helper import xor_mask, xor_numpy
from norec4dna.helper.bin2Quaternary import string2QUATS, quads2dna


class Packet:
    def __init__(self, data, used_packets: typing.Collection[int], total_number_of_chunks: int, read_only: bool = False,
                 seed: int = 0, implicit_mode: bool = True, error_correction: typing.Callable = nocode,
                 packet_len_format: str = "I", crc_len_format: str = "L", number_of_chunks_len_format: str = "I",
                 used_packets_len_format: str = "I", id_len_format: str = "I",
                 save_number_of_chunks_in_packet: bool = True, prepend="", append=""):
        self.data: bytes = data
        self.error_correction: typing.Callable = error_correction
        self.total_number_of_chunks: int = total_number_of_chunks
        self.used_packets: typing.Optional[typing.Collection[int]] = None
        self.set_used_packets(used_packets)
        self.internal_hash = None
        self.degree: int = 0  # just a stub
        self.update_degree()
        self.dna_data: typing.Optional[str] = None
        self.id: int = seed
        self.packet_len_format: str = packet_len_format
        self.crc_len_format: str = crc_len_format
        self.number_of_chunks_len_format: str = number_of_chunks_len_format
        self.used_packets_len_format: str = used_packets_len_format
        self.id_len_format: str = id_len_format
        self.save_number_of_chunks_in_packet: bool = save_number_of_chunks_in_packet
        self.error_prob: typing.Optional[int] = None
        self.implicit_mode: bool = implicit_mode
        self.did_change: bool = False
        self.packed_data: typing.Optional[bytes] = None
        if not read_only:
            self.packed_used_packets: bytes = self.prepare_and_pack()
            self.packed: bytes = self.calculate_packed_data()
        self.prepend = prepend
        self.append = append

    @classmethod
    def from_packet(cls, packet, pseudo=False):
        if not pseudo:
            data = packet.get_data()
        else:
            data = ""
        used_packets = packet.get_used_packets()
        number_of_packets = packet.get_total_number_of_chunks()
        res = cls(data, used_packets, number_of_packets)
        res.error_correction = packet.get_error_correction()
        res.dna_data = None
        return res

    def get_org_class(self):
        return self.__module__.split(".")[1]

    def set_used_packets(self, used_packets: typing.Iterable[int]):
        self.used_packets = used_packets
        self.internal_hash = None

    def prepare_and_pack(self):
        # Format = Highest possible Packetnumber for this file, number of used Packets for this
        # File and the Indizies of the used Packets
        struct_str = "<" + (self.number_of_chunks_len_format if self.save_number_of_chunks_in_packet else "") + (
            self.used_packets_len_format if not self.implicit_mode else "") + self.id_len_format
        if self.save_number_of_chunks_in_packet:
            if self.implicit_mode:
                return struct.pack(struct_str, xor_mask(self.total_number_of_chunks, self.number_of_chunks_len_format),
                                   xor_mask(self.id, self.id_len_format))
            else:
                return struct.pack(struct_str, xor_mask(self.total_number_of_chunks, self.number_of_chunks_len_format),
                                   xor_mask(len(self.used_packets), self.used_packets_len_format),
                                   xor_mask(self.id, self.id_len_format))
        else:
            if self.implicit_mode:
                return struct.pack(struct_str, xor_mask(self.id, self.id_len_format))
            else:
                return struct.pack(struct_str, xor_mask(len(self.used_packets), self.used_packets_len_format),
                                   xor_mask(self.id, self.id_len_format))

    def calculate_packed_data(self):
        # Laenge des Packets + UsedPackets + Data + crc
        self.packed_data = struct.pack("<" + str(len(self.data)) + "s", bytes(self.data))
        payload = struct.pack(
            "<"
            + str(len(self.packed_used_packets))
            + "s"
            + str(len(self.packed_data))
            + "s",
            self.packed_used_packets,
            self.packed_data,
        )
        return self.error_correction(payload)

    def get_struct(self, split_to_multiple_files: bool) -> bytes:
        packed = self.packed
        if not split_to_multiple_files:
            return struct.pack("<" + self.packet_len_format + str(len(packed)) + "s", len(packed), packed)
        else:
            return packed

    def get_dna_struct(self, split_to_multiple_files: bool) -> str:
        if self.dna_data is None and self.error_correction.__name__ == 'dna_reed_solomon_encode':
            self.dna_data = self.prepend + quads2dna(self.get_struct(split_to_multiple_files)) + self.append
        elif self.dna_data is None:
            self.dna_data = self.prepend + "".join(string2QUATS(self.get_struct(split_to_multiple_files))) + self.append
        self.internal_hash = None  # enforce recalculation of hash
        return self.dna_data

    def get_data(self) -> bytes:
        return self.data

    def set_data(self, data: bytes):
        self.data = data

    def get_used_packets(self):
        return self.used_packets

    def get_bool_array_used_packets(self) -> typing.List[bool]:
        return [x in self.used_packets for x in range(self.total_number_of_chunks)]

    def set_bool_array_used_packet(self, b_array: typing.List[bool]):
        assert len(b_array) == self.total_number_of_chunks, "Problem"
        self.set_used_packets(set([k for k in range(self.total_number_of_chunks) if b_array[k]]))

    def get_total_number_of_chunks(self) -> int:
        return self.total_number_of_chunks

    def get_error_correction(self) -> typing.Callable:
        return self.error_correction

    def update_degree(self):
        self.degree = len(self.used_packets)

    def get_degree(self) -> int:
        return self.degree

    def remove_packets(self, packet_set: typing.Set):
        self.set_used_packets(self.used_packets.difference(packet_set))
        self.update_degree()
        self.did_change = True  # CRC is no longer valid

    def xor_and_remove_packet(self, packet):
        self.remove_packets(packet.get_used_packets())
        self.data = xor_numpy(self.data, packet.get_data())
        self.did_change = True  # CRC is no longer valid

    def set_error_prob(self, error_prob: typing.Optional[int] = None):
        self.error_prob = error_prob

    def __str__(self) -> str:
        return (
                "< used_packets: "
                + str(self.used_packets)
                + " , Data: "
                + str(self.data)
                + " , Error Correction: "
                + str(self.error_correction)
                + " >"
        )

    def __repr__(self) -> str:
        return self.__str__()

    def __eq__(self, other) -> bool:
        # if self.error_prob is not None and other.error_prob is not None:
        #    return self.error_prob == other.error_prob
        # else:
        return hash(self) == hash(other)

    def __lt__(self, other) -> bool:
        if self.error_prob is not None and other.error_prob is not None:
            return self.error_prob < other.error_prob
        else:
            return min(self.get_used_packets()) < min(other.get_used_packets())

    def __hash__(self):
        if self.internal_hash is None:
            self.internal_hash = hash(
                str(self.total_number_of_chunks) + str(self.id) + str(
                    self.error_prob) + self.__module__ + self.dna_data)
        return self.internal_hash


class ParallelPacket:
    def __init__(self, used_packets: typing.Set[int], total_number_of_chunks: int, p_id: int, data: bytes,
                 dna_data: str, packed, error_prob: int, packet_len_format: str, crc_len_format: str,
                 number_of_chunks_len_format: str, id_len_format: str, save_number_of_chunks_in_packet: bool,
                 org_class: str = "", prepend="", append="", org_hash=None):
        self.used_packets: typing.Set[int] = used_packets
        self.total_number_of_chunks: int = total_number_of_chunks
        self.id: int = p_id
        self.data: bytes = data
        self.dna_data: str = dna_data
        self.packed: bytes = packed
        self.error_prob: int = error_prob
        self.packet_len_format: str = packet_len_format
        self.crc_len_format: str = crc_len_format
        self.number_of_chunks_len_format: str = number_of_chunks_len_format
        self.id_len_format: str = id_len_format
        self.safe_number_of_chunks_in_packet: bool = save_number_of_chunks_in_packet
        self.org_class: str = org_class.split(".")[1]
        self.bool_arrayused_packets: typing.List[bool] = [
            x in self.used_packets for x in range(0, self.total_number_of_chunks)
        ]
        self.prepend = prepend
        self.append = append
        self.calculated_hash = org_hash

    @classmethod
    def from_packet(cls, packet):
        return ParallelPacket(packet.used_packets, packet.total_number_of_chunks, packet.id, packet.data,
                              packet.dna_data,
                              packet.packed, packet.error_prob, packet.packet_len_format, packet.crc_len_format,
                              packet.number_of_chunks_len_format, packet.id_len_format,
                              packet.save_number_of_chunks_in_packet, org_class=packet.__module__,
                              prepend=packet.prepend, append=packet.append, org_hash=hash(packet))

    def get_org_class(self) -> str:
        return self.org_class

    def __hash__(self):
        if self.calculated_hash is None:
            self.calculated_hash = hash(
                str(self.total_number_of_chunks) + str(self.id) + str(self.error_prob) + self.org_class + self.dna_data)
        return self.calculated_hash

    def __eq__(self, other) -> bool:
        # if self.error_prob is not None and other.error_prob is not None:
        #    return self.error_prob == other.error_prob
        # else:
        return hash(self) == hash(other)

    def __lt__(self, other) -> bool:
        if self.error_prob is not None and other.error_prob is not None:
            return self.error_prob < other.error_prob
        else:
            return min(self.get_used_packets()) < min(other.get_used_packets())

    def __gt__(self, other) -> bool:
        if self.error_prob is not None and other.error_prob is not None:
            return self.error_prob > other.error_prob
        else:
            return min(self.get_used_packets()) > min(other.get_used_packets())

    def get_dna_struct(self, split_to_multiple_files: bool) -> str:
        if self.dna_data is None:
            raise RuntimeError("DNA-Data should not be None!")
        return self.dna_data

    def get_used_packets(self) -> typing.Set[int]:
        return self.used_packets

    def get_struct(self, split_to_multiple_files: bool) -> bytes:
        packed = self.packed
        if not split_to_multiple_files:
            return struct.pack("<" + self.packet_len_format + str(len(packed)) + "s", len(packed), packed)
        else:
            return packed