import socket
import multiprocessing

from multiplexer.MultiInterfaceBase import MultiInterfaceBase
from norec4dna import RU10Encoder, reed_solomon_encode
from norec4dna.distributions.RaptorDistribution import RaptorDistribution
from multiplex_OLD import Multiplexer


class MultiInterfaceFountain(MultiInterfaceBase):
    def __init__(self, port=45155):
        self.__port = port

    def partition_packets(self, packets, paranoid=False):
        """
        TODO: partition the created packets in a way that each channel wont be able to reconstruct the file
        (or even parts if paranoid is true)
        :param packets: packets to partition
        :param paranoid: if set to True, packets over each channel can not be reduced to a part of the file
        :return: list of packets for each used interface
        """
        pass

    def create_ip_socket(self, interface, broadcast=False):
        # create a UDP socket for a given interface to a destination or as a broadcast
        ip = None
        try:
            ip = self.get_ip_address(interface)
        except Exception as ex:
            print(ex)
            return None
        s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM, socket.IPPROTO_UDP)
        s.bind((ip, 0))
        s.setsockopt(socket.SOL_SOCKET, socket.SO_BROADCAST, broadcast)
        return s

    #  @staticmethod
    def send_packets_on_interface(self, packets_socket_dest=None):
        res = []
        try:
            packets, sending_socket, dest = packets_socket_dest
            if sending_socket is None:
                return res
            for packet in packets:
                print("Sending packet from %s:%s to %s:%d" % (
                    sending_socket.getsockname()[0], sending_socket.getsockname()[1], dest, self.__port))
                res.append(sending_socket.sendto(bytearray(packet), (dest, self.__port)))
        except Exception as ex:
            print(ex)
        return res

    def send_packets(self, packets, interfaces, dest=None, broadcast=False):
        """
        sends packets over the specified interface
        :param packets: list of list of packets for each interfaces
        :param interfaces: list of interfaces to use
        :param dest: list of destination addresses in the order based on the interfaces used
        :param broadcast: if true, send as a broadcast on each interface
        :return:
        """
        # TODO if broadcast is false and a list of dests is given,
        #  try to map dest-ips to the correct interfaces according to the netmask
        # TODO maybe send number of chunks as a header but as a XOR produkt of all channels?
        #  ( or use shamirs secret sharing! to get n out of m needed )
        if interfaces is None:
            interfaces = []
        if dest is None:
            dest = "255.255.255.255"
        if broadcast:
            dests = []
            for iface in interfaces:
                try:
                    dests.append(m.get_ip_address(iface, BROADCAST))
                except:
                    dests.append("255.255.255.255")
        elif isinstance(dest, str):
            dests = [dest] * len(interfaces)
        else:
            assert (len(dest) == len(interfaces))
            dests = dest
        sockets = [self.create_ip_socket(interface=iface, broadcast=broadcast) for iface in interfaces]
        p = multiprocessing.Pool(len(interfaces))
        res = p.map(self.send_packets_on_interface, [x for x in zip(packets, sockets, dests)])
        """
        # res = [self.send_packets_on_interface(x) for x in zip(packets, sockets, dest)]
        i = 0
        for iface, dest in zip(interfaces, dests):
            socket = self.create_ip_socket(interface=iface, broadcast=broadcast)
            self.send_packets_on_interface((packets[i], socket, dest))
        """
        return res


if __name__ == "__main__":
    file = "../.INFILES/logo.jpg"
    chunk_size = 100
    number_of_chunks = RU10Encoder.get_number_of_chunks_for_file_with_chunk_size(file, chunk_size)
    INSERT_HEADER = True
    NUM_IN_PACKER = True
    dist = RaptorDistribution(number_of_chunks)
    dna_rules = None  # FastDNARules()
    BROADCAST = False
    error_correction = reed_solomon_encode
    """
    encoder = RU10Encoder(file, number_of_chunks, dist, chunk_size=chunk_size, insert_header=INSERT_HEADER,
                          rules=None,
                          error_correction=error_correction, id_len_format="I", number_of_chunks_len_format="I",
                          save_number_of_chunks_in_packet=NUM_IN_PACKER, mode_1_bmp=False)
    encoder.set_overhead_limit(1.00)
    encoder.encode_to_packets()
    """
    m = MultiInterfaceFountain()
    ifaces = m.list_interfaces()
    clean_ifaces = m.filter_interfaces(ifaces)
    dests = ['192.168.0.105', '192.168.0.87', '192.168.56.1']
    print(clean_ifaces)
    channel_count = min(len(clean_ifaces), len(dests))
    mltp = Multiplexer('../.INFILES/Dorn', channel_count, factor=0.9, chunk_size=chunk_size)
    packet_list_list = mltp.do_multiplex()
    raw_packets_list_list = [[x.get_struct(True) for x in y] for y in packet_list_list]

    # "192.168.0.95"
    m.send_packets(raw_packets_list_list, clean_ifaces, dests, BROADCAST)
    # we might need a two way commuication to know when the receiver finished
