import os
import json
import random
import requests

from norec4dna.helper import bin2Quaternary, quaternary2Bin

# IMPORTANT: if you plan to use this, you should change the MESA_URL to a (local) instance of your own.
MESA_URL = 'http://pc12291.mathematik.uni-marburg.de:5000/api/all'


class DNARules2:
    def __init__(self, active_rules=None):
        self.active_rules = []

    @staticmethod
    def sim_mutation(seq):
        """
        Simulates the mutation of nucleobases with a static probability for a sequence
        :param seq: Sequence to be mutated
        :return: Mutated sequence
        """
        res_seq = []
        for seq in seq:
            mod_seq = ''
            # randint(0, x) <= y defines the chance to mutate: 1 - y/x * 100% chance to mutate
            for x in seq:
                if random.randint(0, 1000) <= 1000:  # always True
                    mod_seq += x
                else:
                    mod_seq += random.choice(('A', 'T', 'G', 'C'))
            res_seq.append(mod_seq)
        return res_seq

    @staticmethod
    def apply_all_rules(packet, from_web=True):
        """
        Calculates the dropchance for either a generated packet or a list of generated packets with the option to
        simulate mutations with a static probability or getting the mutated sequences from the website. Also allows
        to import a config file with configurations downloaded from the website.
        :param packet: The packet or list of packets to calculate the dropchance for
        :param from_web: True: Get the mutated sequences from the website. False: Use the static simulation
        #:param from_file: True: Use a config file for the websiterequest. False: Enter configurations manual in @get_mutated
        :return: The dropchance for a packet or a list of packets
        """
        if type(packet) != list:
            packet = [packet]
        dna_data = []
        for pack in packet:
            dna_data.append(pack.get_dna_struct(True))
        dna_data_mutated = DNARules2.get_mutated(dna_data, from_web)
        res_err = []
        cntr = 0
        for seq in dna_data_mutated:
            seq_org = dna_data[cntr]
            cntr += 1
            dna_data_bin_enc = b""
            for i in range(0, len(seq), 4):
                try:
                    dna_data_bin_enc += quaternary2Bin.quats_to_bytes(seq[i:i + 4])
                except:
                    pass
            dna_data_bin_dec = dna_data_bin_enc
            dna_data_dec = ''
            for x in dna_data_bin_dec:
                dna_data_dec += bin2Quaternary.byte2QUATS(x)
            changes = sum([1 if dna_data_dec[x] != seq_org[x] else 0 for x in
                           range(min(len(dna_data_dec), len(seq_org)))])
            if changes < 3 - abs(len(dna_data_dec) - len(seq_org)):
                dropchance = 10 * changes
                for i in range(0, len(dna_data_bin_dec)):
                    if dna_data_bin_dec[i] != dna_data_bin_enc[i] and dropchance <= 95:
                        dropchance += 5
                res_err.append(dropchance / 100)
            else:
                res_err.append(1.0)
        if type(res_err) == list and len(res_err) == 1:
            res_err = res_err[0]
        return res_err

    # executes requests to the website with given parameters
    @staticmethod
    def get_mutated_from_web(seq, config=None, json_config=None):
        """
        Builds and executes the request for the website from given configurations, either a config file or manual
        added parameters
        :param seq: The sequence(s) to mutate
        :param config: If used, the config file downloaded from the website
        :param json_config: If used, the manually generated json_config
        :return: The response from the website as dictionary
        """
        header = {'content-type': 'application/json;charset=UTF-8'}
        if json_config:
            payload = json_config
        else:
            payload = config
        payload['sequence'] = seq
        payload['asHTML'] = False
        res = requests.post(MESA_URL, data=json.dumps(payload), headers=header)
        try:
            return res.json()
        except:
            print("Error")

    @staticmethod
    def get_mutated(seq, from_web, json_config=None):
        """
        Takes the sequences to mutate and calls either @sim_mutation or @get_mutated_from_web. If the website is used, the
        results are also processed to get only the mutated sequences for the input
        :param seq: The sequence(s) to mutate
        :param from_web: True: Use the website, False: Use the simulated mutation
        :param json_config: If used, the manually generated json_config
        :return: A list of mutated sequences for the input
        """
        if type(seq) == str:
            seq = [seq]
        if from_web:
            if json_config:
                res_all = DNARules2.get_mutated_from_web(seq=seq, json=json.dumps(json_config))
            else:
                try:
                    file = os.environ['dna_sim_config']
                except:
                    print("Could not find ENV-Var 'dna_sim_config', falling back to 'mosla.json'")
                    file = 'mosla.json'
                with open(file) as json_file:
                    config = json.load(json_file)
                    res_all = DNARules2.get_mutated_from_web(seq=seq, config=config)
            res_seq = []
            for seq in seq:
                res_seq.append(res_all[seq]['res']['modified_sequence'])
            return res_seq
        else:
            res_seq = DNARules2.sim_mutation(seq)
        return res_seq


if __name__ == "__main__":
    sequence = 'AAAACCCCGGGGTTTT'
    print('Test with the mosla.json file for: ' + sequence)
    print(DNARules2.get_mutated(sequence, True))
