import socket
import struct
import threading
from queue import Queue

from MultiInterfaceBase import MultiInterfaceBase
from norec4dna import RU10Decoder, reed_solomon_decode
from norec4dna.helper import xor_mask


class MultiInterfaceReceiver(MultiInterfaceBase):
    def __init__(self, listen_ifaces=None, port=45155):
        self.listen_ifaces = listen_ifaces
        self.__port = port

    def get_port(self):
        return self.__port

    def create_listen_socket(self, interface, broadcast=False):
        sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        ip = self.get_ip_address(interface, broadcast)
        if ip == "0.0.0.0":
            return None
        sock.bind((ip, self.__port))
        print("Listening on %s:%s" % (ip, m_r.get_port()))
        # sock.setblocking(False)
        sock.settimeout(1)
        return sock


def listen(sock, queue, signals):
    if sock is None:
        return None
    while not signals["shutdown"]:
        try:
            data, sender = sock.recvfrom(1024)
            a, b = struct.unpack("<II", data[0:8])
            print("Packet from: %s:%s to %s - #Chunks: %s - Id: %s" % (
                sender[0], sender[1], sock.getsockname()[0], xor_mask(a), xor_mask(b)))
            queue.put(data)
        except socket.error as ex:
            # queue.put(e)
            continue


if __name__ == "__main__":
    m_r = MultiInterfaceReceiver()
    BROADCAST = False  # TODO fix broadcast not working for Windows / socket receiving broadcast on bound normal ip
    USE_HEADER_CHUNK = True
    error_correction = reed_solomon_decode
    decoder = RU10Decoder(None, use_headerchunk=USE_HEADER_CHUNK, error_correction=error_correction,
                          static_number_of_chunks=None)
    decoder.read_all_before_decode = True
    ifaces = m_r.list_interfaces()
    clean_ifaces = m_r.filter_interfaces(m_r.list_interfaces())
    socks = []
    for x in clean_ifaces:
        try:
            socks.append(m_r.create_listen_socket(x, BROADCAST))
        except Exception as e:
            if hasattr(x, 'decode'):
                x = x.decode()
            print("<%s>: %s" % (x, e))
            raise e
    pqueue = Queue()
    signals = {"shutdown": False}
    for sock in socks:
        try:
            thread = threading.Thread(target=listen, args=(sock, pqueue, signals,))
            thread.start()
        except socket.error as e:
            print(e)
    num = 0
    while True:
        try:
            packet_strs = []
            while True:
                try:
                    packet_str = pqueue.get(timeout=2)
                    packet_strs.append(packet_str)
                    if len(packet_strs) > 50 * len(socks):
                        break  # insert into decoder every 20 packet
                except Exception as ex:
                    print(ex)
                    break  # ... or if a timeout occurs
            for packet_str in packet_strs:
                pack = decoder.parse_raw_packet(packet_str, crc_len_format="L", number_of_chunks_len_format="I",
                                                packet_len_format="I", id_len_format="I")
                decoder.input_new_packet(pack)
            if decoder.GEPP is not None and decoder.solve():
                print("Success!")
                signals["shutdown"] = True
                pqueue.task_done()
                decoder.saveDecodedFile(null_is_terminator=False, print_to_output=False)
                # TODO maybe signal the sender that we are finished
                #  -> to prevent attacker from closing the connection we should use Shamirs secret sharing
                break
        except Exception as e:
            print(e)
            # raise e
