import math
import numpy as np
import matplotlib.pyplot as plt
from norec4dna.RU10Decoder import RU10Decoder
from norec4dna.OnlineDecoder import OnlineDecoder
from norec4dna.RU10Packet import RU10Packet
from norec4dna.Packet import Packet
from norec4dna.OnlinePacket import OnlinePacket
from norec4dna.helper.RU10Helper import from_true_false_list
from norec4dna.LTDecoder import LTDecoder


class AutomatedFindMinimum:

    def __init__(self, packets, error_correction, dist, use_headerchunk=False):
        self.packets = packets
        self.error_correction = error_correction
        self.dist = dist
        self.use_headerchunk = use_headerchunk
        self.org_normal_packets = [self.parallel_to_normal(packet, error_correction, dist) for packet in self.packets]
        self.normal_packets = self.org_normal_packets
        self.decoder = self.init_dec()
        self.underrep_chunks = []
        self.org_decoding_packets = self.get_decoding_packets()

    @staticmethod
    def parallel_to_normal(par_packet, error_correction, dist, quality=7, epsilon=0.06):
        """
        Converts parallel packets either to RU10Packets, normal Packets or OnlinePackets based on the original class of
        the packet.
        :param par_packet:
        :param error_correction:
        :param dist:
        :param quality:
        :param epsilon:
        :return: Converted packet
        """
        if par_packet.get_org_class() == "RU10Packet":
            packet = RU10Packet(data=par_packet.data, used_packets=par_packet.used_packets,
                                total_number_of_chunks=par_packet.total_number_of_chunks, id=par_packet.id, dist=dist,
                                error_correction=error_correction, packet_len_format=par_packet.packet_len_format,
                                crc_len_format=par_packet.crc_len_format,
                                number_of_chunks_len_format=par_packet.number_of_chunks_len_format,
                                id_len_format=par_packet.id_len_format,
                                save_number_of_chunks_in_packet=par_packet.safe_number_of_chunks_in_packet)
            packet.error_prob = par_packet.error_prob
        elif par_packet.get_org_class() == "Packet":
            packet = Packet(data=par_packet.data, used_packets=par_packet.used_packets,
                            total_number_of_chunks=par_packet.total_number_of_chunks, error_correction=error_correction,
                            packet_len_format=par_packet.packet_len_format, crc_len_format=par_packet.crc_len_format,
                            number_of_chunks_len_format=par_packet.number_of_chunks_len_format,
                            id_len_format=par_packet.id_len_format)
            packet.error_prob = par_packet.error_prob
        elif par_packet.get_org_class() == "OnlinePacket":
            packet = OnlinePacket(data=par_packet.data, used_packets=par_packet.used_packets,
                                  check_block_number=par_packet.id,
                                  total_number_of_chunks=par_packet.total_number_of_chunks,
                                  error_correction=error_correction,
                                  dist=dist, quality=quality, epsilon=epsilon, crc_len_format=par_packet.crc_len_format,
                                  number_of_chunks_len_format=par_packet.number_of_chunks_len_format,
                                  check_block_number_len_format=par_packet.id_len_format)
            packet.error_prob = par_packet.error_prob
        else:
            raise RuntimeError("Unsupported packet type")
        return packet

    def init_dec(self):
        """
        Gets the pseudodecoder based on the original class of the first packet in the packet list.
        :return: Pseudodecoder
        """
        packet = self.normal_packets[0]
        if packet.get_org_class() == "RU10Packet":
            return RU10Decoder.pseudo_decoder(packet.total_number_of_chunks)
        elif packet.get_org_class() == "OnlinePacket":
            return OnlineDecoder.pseudo_decoder(packet.total_number_of_chunks)
        elif packet.get_org_class() == "Packet":
            return LTDecoder.pseudo_decoder(packet.total_number_of_chunks)

    @staticmethod
    def get_packet_dist(packets_bool_arrays, number_of_chunks):
        """
        Takes a list of bool arrays from packets and calculates the occurence of each chunk in these packets.
        :param packets_bool_arrays: list of bool_arrays of the packets
        :param number_of_chunks:
        :return: List of chunks and usage counter
        """
        used_chunks = np.array([0 for _ in range(number_of_chunks)])
        for arr in packets_bool_arrays:
            i = 0
            for bo in arr:
                if bo == True:
                    used_chunks[i] += 1
                i += 1
        return used_chunks

    def get_decoding_packets(self):
        """
        Finds all packets needed to decode the file.
        :return: Decoding packets
        """
        for packet in self.normal_packets:
            self.decoder.input_new_packet(packet)
        self.decoder.solve()
        decoding_packet_indices = self.decoder.GEPP.packet_mapping[:self.decoder.GEPP.result_mapping.size]
        decoding_packets = []
        for x in decoding_packet_indices:
            packet = self.normal_packets[x[0]]
            decoding_packets.append(packet)
        for pack in decoding_packets:
            self.normal_packets.remove(pack)
        return decoding_packets

    def get_bool_arrays(self, packets):
        """
        Takes a list of packets and gets their bool_arrays of used chunks.
        :param packets:
        :return: List of bool_arrays for the given packets
        """
        bool_arr = []
        if isinstance(packets[0], RU10Packet):
            for pack in packets:
                bool_arr.append(self.decoder.removeAndXorAuxPackets(pack))
        else:
            for pack in packets:
                bool_arr.append(pack.get_bool_array_used_packets())
        return bool_arr

    @staticmethod
    def packet_already_exists(new_packet_bool_arr, used_packets_boolarr):
        """
        Checks if a similar packet with a different seed already exists in the used packets.
        :param new_packet_bool_arr:
        :param used_packets_boolarr:
        :return:
        """
        for arr in used_packets_boolarr:
            if np.array_equal(new_packet_bool_arr, arr):
                return True
        return False

    def set_underrep_chunks(self, used_packets_boolarr, number_of_chunks, min, prio_chunks=None):
        """
        Calculates the chunks that are underrepresented in the used packets based on a given minimum occurence. Saves a
        list of underrepresented chunk with the chunknumber and the difference of their occurence and the min for the
        occurence. If prioritized chunks are given they will be added to the list.
        :param used_packets_boolarr:
        :param number_of_chunks:
        :param min:
        :param prio_chunks:
        """
        if prio_chunks is None:
            prio_chunks = []
        cnt = 0
        used_chunks = self.get_packet_dist(used_packets_boolarr, number_of_chunks)
        for chunk in used_chunks:
            if chunk < min:
                self.underrep_chunks.append([cnt, min - chunk])
            cnt += 1
        if len(prio_chunks) > 0:
            for chunk in prio_chunks:
                for underrep_chunk in self.underrep_chunks:
                    if chunk[0] == underrep_chunk[0]:
                        underrep_chunk[1] += chunk[1]
                    else:
                        self.underrep_chunks.append([chunk[0], chunk[1]])

    def get_packets_for_chunk(self, chunk):
        """
        Searches for all packets that contain the given chunk.
        :param chunk:
        :return: All packets that contain the chunk
        """
        new_packets = set()
        for packet in self.normal_packets:
            if chunk in self.get_used_packets(packet):
                new_packets.add(packet)
        return new_packets

    def get_packets_for_all_chunks(self):
        """
        Takes a list of underrepresented chunks and gets all packets that contain these chunks.
        :return: Set of packets that contain at least one of the chunks
        """
        possible_packets = set()
        for chunk in self.underrep_chunks:
            possible_packets.update(self.get_packets_for_chunk(chunk[0]))
        return possible_packets

    def get_optimal_packets(self, decoding_packets, used_packets_boolarr, number_of_chunks, no_new_packets,
                            err_prob_fac, pack_min_occ, prio_chunks=None):
        """
        Computes the optimal packets to add to the decoding packets with all given parameters. See automated_optimization
        for further information about the parameters. Also checks the possible packets for their error_prob and similar
        existing packets.
        :param decoding_packets:
        :param used_packets_boolarr:
        :param number_of_chunks:
        :param no_new_packets:
        :param err_prob_fac:
        :param pack_min_occ:
        :param prio_chunks:
        :return:
        """
        self.set_underrep_chunks(used_packets_boolarr, number_of_chunks, pack_min_occ, prio_chunks)
        possible_packets = self.get_packets_for_all_chunks()
        avg_err_prob = sum(x.error_prob for x in decoding_packets) / len(decoding_packets)
        possible_packets = self.discard_already_exists(possible_packets, used_packets_boolarr)
        possible_packets = self.discard_packets(possible_packets, err_prob_fac, avg_err_prob)
        additional_packets = self.compute_packets(possible_packets, no_new_packets)
        return additional_packets

    def discard_already_exists(self, packets, used_packets_boolarr):
        """
        Takes a list of packets and discards all packets that already exists with different seeds/ids.
        :param packets:
        :param used_packets_boolarr:
        :return:
        """
        tmp = []
        for packet in packets:
            if isinstance(packet, RU10Packet):
                if self.packet_already_exists(self.decoder.removeAndXorAuxPackets(packet), used_packets_boolarr):
                    tmp.append(packet)
            elif self.packet_already_exists(packet.get_bool_array_used_packets(), used_packets_boolarr):
                tmp.append(packet)
        if len(tmp) > 0:
            print("-> Packets discarded since similar packets are being used: ")
            already_existing = ""
            for pack in tmp:
                packets.remove(pack)
                already_existing += str(pack.id) + " "
            print(already_existing)
        return packets

    def compute_pack(self, possible_packets, underrep_chunks):
        """
        Recursively computes packets that contain the most underrepresented chunks and returns the packet with the
        lowest error_prob that contains the most of the chunks. Takes a copy of underrep_chunks to update it in every
        iteration without manipulating the original list.
        :param possible_packets:
        :param underrep_chunks:
        :return:
        """
        if len(underrep_chunks) > 0:
            underrep_chunks.sort(key=lambda x: x[1], reverse=True)
            new_possible_packets = [x for x in possible_packets if underrep_chunks[0][0] in self.get_used_packets(x)]
            if len(new_possible_packets) == 0:
                return possible_packets
            else:
                underrep_chunks.pop(0)
                possible_packets = self.compute_pack(new_possible_packets, underrep_chunks)
                if len(underrep_chunks) > 1:
                    underrep_chunks.pop(0)
                    possible_packets = self.compute_pack(possible_packets, underrep_chunks)
                return possible_packets
        else:
            return possible_packets

    @staticmethod
    def discard_packets(possible_packets, err_prob_fac, avg_err_prob):
        """
        Discards packets based on the average error_prob of the decoding packets and an error_prob factor given by the
        user.
        :param possible_packets:
        :param err_prob_fac:
        :param avg_err_prob:
        :return: Possible packets with allowed error_probs
        """
        tmp = []
        err_prob_fac += 1
        for pack in possible_packets:
            if pack.error_prob > avg_err_prob * err_prob_fac:
                tmp.append(pack)
        if len(tmp) > 0:
            print("-> Packets discarded due to high error_prob (above " + str(
                err_prob_fac) + " times higher than the avg): ")
            ids = ""
            for pack in tmp:
                possible_packets.remove(pack)
                ids += str(pack.id) + " "
            print(ids)
        return possible_packets

    def compute_packets(self, possible_packets, no_new_packets):
        """
        Computes no_new_packets new packets. In every iteration compute_pack returns a list of "optimal" packets and the
        one with the lowest error_prob is taken. The list of underrep_chunks is updated also.
        :param possible_packets:
        :param no_new_packets:
        :return:
        """
        packets = []
        for _ in range(no_new_packets):
            if len(self.underrep_chunks) > 0:
                new_packets = self.compute_pack(possible_packets, self.underrep_chunks.copy())
                if len(new_packets) > 0:
                    new_pack = min(new_packets, key=lambda x: x.error_prob)
                    packets.append(new_pack)
                    possible_packets.remove(new_pack)
                    for pack in self.get_used_packets(new_pack):
                        rem_ch = []
                        for chunk in self.underrep_chunks:
                            if chunk[0] == pack:
                                chunk[1] -= 1
                                if chunk[1] < 1:
                                    rem_ch.append(chunk)
                        for chunk in rem_ch:
                            self.underrep_chunks.remove(chunk)
                else:
                    print(
                        "-> There are no more packets to add, but the minimum occurence for each chunk may not be reached. "
                        "Please check your limitations")
                    break
            else:
                # print("Every chunk has rechead the occurence minimum without the use of the complete overhead limit.")
                break
        return packets

    """@staticmethod
    def print_packet(packet):
        used = ""
        for x in automated_find_minimum.get_used_packets(packet):
            used += x + " "
        print("Packet-ID: " + str(packet.id) + "\nError-Prob: " + str(packet.error_prob) + "\nUsed Chunks: " + used)"""

    def get_used_packets(self, packet):
        """
        Gets the used chunks/packets of a packet and returns a list containing their numbers.
        :param packet:
        :return:
        """
        if isinstance(packet, RU10Packet):
            arr = self.decoder.removeAndXorAuxPackets(packet)
        else:
            arr = packet.get_bool_array_used_packets()
        return from_true_false_list(arr)

    @staticmethod
    def calc_used_chunks(number_of_chunks, bool_arr):
        """
        Calculates how often a chunk is used in the packets the bool_arrays are taken from.
        :param number_of_chunks:
        :param bool_arr:
        :return: A list with the occurence of every chunk in the packets. (index = chunk)
        """
        used_chunks = np.array([0 for _ in range(number_of_chunks)])
        for arr in bool_arr:
            for x in range(number_of_chunks):
                if arr[x] == True:
                    used_chunks[x] += 1
        return used_chunks

    def get_chunk_err_probs(self, packets, number_of_chunks):
        """
        Gets the err_probs of all packets containing a chunk and returns a list of lists containing these err_probs for
        every chunk.
        :param packets:
        :param number_of_chunks:
        :return: List of error_probs
        """
        err_probs = [[] for _ in range(number_of_chunks)]
        for x in range(number_of_chunks):
            err_probs[x] = []
            for packet in packets:
                if x in self.get_used_packets(packet):
                    err_probs[x].append(packet.error_prob)
        return err_probs

    @staticmethod
    def plot_with_average(used_chunks, err_probs, title):
        """
        Plots the chunks with their usage and the average err_prob of the packets containing the chunks and saves the
        plot as png.
        :param used_chunks:
        :param err_probs:
        :param title:
        """
        fig, ax1 = plt.subplots()
        fig.suptitle(title, y=1.1)
        ax1.set_xlabel('Chunk')
        ax1.set_ylabel('Occurence', color='r')
        ax1.plot(used_chunks, color='r')
        ax1.tick_params(axis='y', labelcolor='r')
        ax1.tick_params(axis='x', labelrotation=45)
        ax2 = ax1.twinx()
        ax2.set_ylabel('Avg err_prob')
        avg_err = []
        for err_prob in err_probs:
            avg_err.append(sum(err_prob) / len(err_prob))
        ax2.plot(avg_err, '--k')
        plt.grid()
        plt.xticks(rotation=45)
        plt.savefig(title + '.png', dpi=1000)
        plt.show()

    @staticmethod
    def plot_err_probs(used_chunks, err_probs, title):
        """
        Plots the chunks with their usage and boxplots for the err_probs and saves the plot as png.
        :param used_chunks:
        :param err_probs:
        :param title:
        """
        fig, ax1 = plt.subplots()
        fig.suptitle(title, y=1.1)
        plt.xticks(rotation=90)
        color_graph = "chocolate"
        color_boxplot = "lightgrey"
        ax1.set_xlabel('Chunk')
        ax1.set_ylabel('Occurence', color=color_graph)
        ax2 = ax1.twinx()
        ax1.plot(used_chunks, color=color_graph)
        bp = ax2.boxplot(err_probs, patch_artist=True, vert=True, positions=range(len(used_chunks)))
        ax1.tick_params(axis='y', labelcolor=color_graph)
        ax2.set_ylabel('Error probability')
        plt.setp(bp['medians'], color="firebrick")
        for patch in bp['boxes']:
            patch.set_facecolor(color_boxplot)
        plt.grid()
        plt.savefig(title + '.png', dpi=1000)
        ax1.set_zorder(ax2.get_zorder() + 1)
        ax1.patch.set_visible(False)
        plt.show()

    @staticmethod
    def plot_chunk_usage(bool_arrs, number_of_chunks, title, print_usage=False):
        """
        Plots the chunk usage for different variants in one plot and saves the plot as png
        :param bool_arrs:
        :param number_of_chunks:
        :param title:
        :param print_usage:
        """
        fig, ax = plt.subplots()
        fig.suptitle(title, y=1.1)
        ax.set_xlabel('Chunk')
        ax.set_ylabel('Occurence')
        for ent in bool_arrs:
            used_chunks = AutomatedFindMinimum.calc_used_chunks(number_of_chunks, ent[0])
            ax.plot(used_chunks, color=ent[1])
            if print_usage:
                print("Chunk occurence in packets:" + str(used_chunks))
        plt.grid()
        plt.xticks(rotation=45)
        plt.savefig(title + '.png', dpi=1000)
        plt.show()

    def automated_optimization(self, overhead_lim, overhead_fac, err_prob_fac, chunk_min_occ=None, plot=False,
                               plot_avg=False, prio_chunks=None):
        """
        Calculates and adds the optimal packets to meet the limitations given. The number of new packets is calculated
        with the overhead_lim. If these new packets are not enough, the overhead_fac is the maximum overhead of the
        overhead. The err_prob_fac is the factor that defines the max allowed err_prob for new packets, which is
        calculated with (avg_err_prob * (1+err_prob_fac)). The method finishes if either every chunk is in chunk_min_occ
        packets or no more packets can be added.
        :param overhead_lim:
        :param overhead_fac:
        :param err_prob_fac:
        :param chunk_min_occ:
        :param plot:
        :param plot_avg:
        :param prio_chunks:
        :return:
        """
        if prio_chunks is None:
            prio_chunks = []
        print("----- Automated optimization -----")
        # Get decoding packets and their bool_arrays
        self.normal_packets = self.org_normal_packets.copy()
        decoding_packets = self.org_decoding_packets.copy()
        number_of_chunks = decoding_packets[0].total_number_of_chunks
        decoding_packets_bool_arrays = self.get_bool_arrays(decoding_packets)
        if plot:
            plotting = [(decoding_packets_bool_arrays, 'firebrick')]
            used_chunks = self.calc_used_chunks(number_of_chunks, decoding_packets_bool_arrays)
            decoding_err_probs = self.get_chunk_err_probs(decoding_packets, number_of_chunks)
            if plot_avg:
                self.plot_with_average(used_chunks, decoding_err_probs, "opt_decoding")
            else:
                self.plot_err_probs(used_chunks, decoding_err_probs, "opt_decoding")
        # Get and add the optimal new packets without exceeding the overhead_limit
        if chunk_min_occ is None:
            packet_dist = self.get_packet_dist(decoding_packets_bool_arrays, number_of_chunks)
            _sum = 0
            for occ in packet_dist:
                _sum += occ
            chunk_min_occ = math.floor((_sum / len(packet_dist)) * 0.9)
            print("-> No minimum chunk occurence set. The calculated values is " + str(
                chunk_min_occ) + " (90% of the average occurence).")
        no_new_packets = math.floor(number_of_chunks * overhead_lim)
        add_packets = self.get_optimal_packets(decoding_packets, decoding_packets_bool_arrays,
                                               number_of_chunks,
                                               no_new_packets, err_prob_fac, chunk_min_occ, prio_chunks)
        print("-> Number of packets needed to decode: " + str(len(decoding_packets)))
        for pack in add_packets:
            decoding_packets.append(pack)
        print("-> Number of packets using the overhead to optimize: " + str(len(decoding_packets)))
        # Check if any chunks are still underrepresented and use the overhead_factor to fix it if possible
        new_decoding_packets_bool_arrays = self.get_bool_arrays(decoding_packets)
        if plot:
            plotting.append((new_decoding_packets_bool_arrays, 'orange'))
            new_decoding_err_probs = self.get_chunk_err_probs(decoding_packets, number_of_chunks)
            used_chunks = self.calc_used_chunks(number_of_chunks, new_decoding_packets_bool_arrays)
            if plot_avg:
                self.plot_with_average(used_chunks, new_decoding_err_probs, "opt_overhead")
            else:
                self.plot_err_probs(used_chunks, new_decoding_err_probs, "opt_overhead")
        if len(self.underrep_chunks) > 0:
            print("-> Some chunks are still underrepresented. Using the overhead factor now.")
            overhead_ex = math.floor(no_new_packets * overhead_fac)
            if overhead_ex == 0:
                print("-> The overhead factor doesn't allow additional packets. Please check your parameters.")
            else:
                exceeding_packets = self.get_optimal_packets(decoding_packets, new_decoding_packets_bool_arrays,
                                                             number_of_chunks, overhead_ex, err_prob_fac, chunk_min_occ)
                for pack in exceeding_packets:
                    decoding_packets.append(pack)
                print("-> Number of packets exceeding the overhead limit based on the overhead_fac: " + str(
                    len(decoding_packets)))
                exceeding_packets_bool_arrays = self.get_bool_arrays(decoding_packets)
                if plot:
                    plotting.append((exceeding_packets_bool_arrays, 'forestgreen'))
                    used_chunks = self.calc_used_chunks(number_of_chunks, exceeding_packets_bool_arrays)
                    exceeding_err_probs = self.get_chunk_err_probs(decoding_packets, number_of_chunks)
                    if plot_avg:
                        self.plot_with_average(used_chunks, exceeding_err_probs, "opt_exceeding")
                    else:
                        self.plot_err_probs(used_chunks, exceeding_err_probs, "opt_exceeding")
                if len(self.underrep_chunks) > 0:
                    print(
                        "-> There are still underrepresented chunks, but no more packets can be added with the given limitations.")
        if plot:
            self.plot_chunk_usage(plotting, number_of_chunks, "opt_all")
        print("----- Finished automated optimization -----")
        return decoding_packets

    def interactive_optimization(self):
        """
        Interactive optimization using inputs to get all params.
        :return:
        """
        print("----- Entered interactive optimization -----")
        overhead_lim = float(input("-> Set an overhead limit (0.1 = 10%): "))
        overhead_fac = float(input(
            "-> If you want to allow exceeding these limit if necessary to optimize the chunk occurence set a factor, otherwise type 0: "))
        err_prob_fac = float(
            input("-> How far above the average may the err_prob of the new packets be? (0.1 = 10%): "))
        chunk_min_occ = int(input("-> Chunk min occurence? (0 means 90% of the average): "))
        if chunk_min_occ == 0:
            chunk_min_occ = None
        optimal_packets = self.automated_optimization(overhead_lim, overhead_fac, err_prob_fac, chunk_min_occ,
                                                      plot=True)
        repeat = input("-> Do you want to repeat the calculation with new params? (y/n): ")
        if repeat == 'y':
            optimal_packets = self.interactive_optimization()
        print("----- Exited interactive optimization -----")
        return optimal_packets
