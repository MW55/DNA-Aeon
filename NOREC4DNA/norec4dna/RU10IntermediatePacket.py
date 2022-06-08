from norec4dna.RU10Packet import RU10Packet


class RU10IntermediatePacket(RU10Packet):
    def __init__(self, data, used_packets, total_number_of_chunks, id, dist=None):
        super().__init__(data, used_packets, total_number_of_chunks, id, dist=dist)
        self.data: bytes = data
        self.id: int = id
        self.set_used_packets(used_packets)
        self.total_number_of_chunks: int = total_number_of_chunks


if __name__ == "__main__":
    print("This class must not be called by itself.")
