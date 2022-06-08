import multiprocessing
from functools import partial
from io import BytesIO

from norec4dna import RU10Decoder
from norec4dna.ErrorCorrection import reed_solomon_decode
from norec4dna.helper.quaternary2Bin import tranlate_quat_to_byte

STATIC_NUMBER_OF_CHUNKS = 338884
spare1core = True
PARALLEL = False


def split(a, n):
    k, m = divmod(len(a), n)
    return (a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))


def reconstruct_Packets(lst):
    """
    norepair_symbols = 3
    decoder = RU10Decoder(file=None,
                          error_correction=lambda x: reed_solomon_decode(x, norepair_symbols),
                          use_headerchunk=False, static_number_of_chunks=STATIC_NUMBER_OF_CHUNKS)
    decoder.number_of_chunks = STATIC_NUMBER_OF_CHUNKS
    """
    packet_list = []
    i = 0
    for error_prob, seed, dna_str in lst:
        packet = decoder.parse_raw_packet(BytesIO(tranlate_quat_to_byte(dna_str)).read(), crc_len_format="L",
                                          number_of_chunks_len_format="",
                                          packet_len_format="H",
                                          id_len_format="I")
        packet_list.append((decoder.removeAndXorAuxPackets(packet), packet.get_data()))
        if i % 100 == 0:
            print(str(i))
        i += 1
    return packet_list


infile = "/home/michael/Schreibtisch/MOSLA/Kuehe_WD/out.fasta"
with open(infile, "r") as in_file:
    norepair_symbols = 3
    decoder = RU10Decoder(file=None,
                          error_correction=lambda x: reed_solomon_decode(x, norepair_symbols),
                          use_headerchunk=False, static_number_of_chunks=STATIC_NUMBER_OF_CHUNKS)
    decoder.number_of_chunks = STATIC_NUMBER_OF_CHUNKS
    raw_packet_list = []
    while True:
        line = in_file.readline()
        if not line:
            break
        error_prob, seed = line[1:].replace("\n", "").split("_")
        line = in_file.readline()
        if not line:
            break
        dna_str = line.replace("\n", "")
        raw_packet_list.append((error_prob, seed, dna_str))
        if not PARALLEL:
            packet = decoder.parse_raw_packet(BytesIO(tranlate_quat_to_byte(dna_str)).read(), crc_len_format="L",
                                              number_of_chunks_len_format="",
                                              packet_len_format="H",
                                              id_len_format="I")
            res = decoder.input_new_packet(packet)
            if decoder.GEPP.n % 100 == 0:
                print("Parsed packet " + str(seed) + " - " + str(decoder.GEPP.n))
            if res:
                decoder.saveDecodedFile(last_chunk_len_format="H", null_is_terminator=False, print_to_output=False)
                break
    if PARALLEL:
        norepair_symbols = 3
        decoder = RU10Decoder(file=None,
                              error_correction=lambda x: reed_solomon_decode(x, norepair_symbols),
                              use_headerchunk=False, static_number_of_chunks=STATIC_NUMBER_OF_CHUNKS)
        decoder.number_of_chunks = STATIC_NUMBER_OF_CHUNKS
        packet_list = []
        i = 0
        error_prob, seed, dna_str = raw_packet_list[0]
        packet = decoder.parse_raw_packet(BytesIO(tranlate_quat_to_byte(dna_str)).read(), crc_len_format="L",
                                          number_of_chunks_len_format="",
                                          packet_len_format="H",
                                          id_len_format="I")
        cores = multiprocessing.cpu_count()
        if spare1core:
            cores = cores - 1
        p = multiprocessing.Pool(cores)
        list_of_lists = p.map(partial(reconstruct_Packets), split(raw_packet_list, cores))
        l = []
        map(l.extend, list_of_lists)
        print(len(l))
