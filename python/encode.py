import subprocess
import argparse
import pathlib
import json
import os
import sys

def header_crc_mapper(header_crc_conf_entry, header_bool):
    #if not header_bool and header_crc_conf_entry:
    #    raise ValueError("Header is turned off but header crc length is not 0. Please adjust.")
    #    exit(1)
    crc_choices = {0: None, 8:"B", 16:"H", 32:"I"}
    try:
        return crc_choices[header_crc_conf_entry]
    except KeyError:
        raise ValueError("Allowed choices for the header crc length are 0 bits, 8 bits, 16 bits or 32 bits. Aborting")
        exit(1)

def encode_norec_for_ac(config_data, current_path, mode="subprocess"):
    """
    outer encoder for the fountaincode-arithmetic code concatenation.
    :param file:
    :return:
    """
    input_file = config_data["encode"]["input"]
    chunk_size = config_data["NOREC4DNA"]["chunk_size"]
    packet_redundancy = config_data["NOREC4DNA"]["package_redundancy"]
    header = " --insert_header " if config_data["NOREC4DNA"]["insert_header"] else ""
    error_detection = config_data["NOREC4DNA"]["error_detection"]
    header_crc_str = "" if not header_crc_mapper(config_data["NOREC4DNA"]["header_crc_length"], config_data["NOREC4DNA"]["insert_header"]) else " --header_crc_str " + header_crc_mapper(config_data["NOREC4DNA"]["header_crc_length"], config_data["NOREC4DNA"]["insert_header"]) + "" 
    filename = input_file.split("/")[-1]
    py_command = ("{cpath}/libraries/NOREC4DNA/venv/bin/python3 {cpath}/libraries/NOREC4DNA/demo_raptor_encode.py --chunk_size {chunk_size_str} --error_correction {err_det} --save_as_zip {file}{ins_header}{crc_str} --overhead {redundancy}".format(
        cpath=current_path, 
        chunk_size_str=str(chunk_size),
        file=input_file,
        redundancy=packet_redundancy,
        ins_header=header,
        err_det=error_detection,
        crc_str=header_crc_str))
    if (mode == "subprocess"):  
        process = subprocess.Popen(py_command.split(), stdout=subprocess.PIPE)
    elif (mode == "sys"):
        process = subprocess.Popen(py_command.split(), stdout=sys.stdout, stderr=sys.stderr)
    else:
        print("Invalid mode. Exiting...")
        exit(1)
    output, error = process.communicate()
    if __debug__ and False:
        print(output)
        # it seems that output works only is the mode is set to "subprocess"
    norec_config = output.split()[-1].decode() # hardcoded
    #norec_config = output.split()[-1].decode() # this give an error is the number of Chunks is not big enough
    # print("\n" + NOREC4DNA_BASE_PATH + "/" + filename + ".ini\n")
    # os.rename(config, NOREC4DNA_BASE_PATH + "/" + filename + ".ini")
    with open(norec_config, "r", encoding="utf-8") as c_:
        line_list = c_.readlines()
        if __debug__ and False: 
            print("\n\nline list: [{cpath}/configs/data/{fname}_RU10.zip]".format(cpath=current_path, fname=filename))
        line_list[0] = "[{cpath}/configs/data/{fname}_RU10.zip]\n".format(cpath=current_path, fname=filename) #"[../data/" + filename + "_RU10.zip]\n"
        if config_data["NOREC4DNA"]["header_crc_length"] == 0:
            del(line_list[-6])
        with open("{cpath}/{ini_path}".format(cpath=current_path, ini_path=config_data["decode"]["NOREC4DNA_config"]), "w") as o_:
            o_.writelines(line_list)
    if __debug__ and False:
        print(f"norec config file: {norec_config}")
    os.remove(norec_config)
    return

def encode_ac(current_path, config, mode="subprocess"):
    """
    encode the sequences using the fountaincode-arithmetic code concatenation.

    it call the C++ binary arithmetic_modulator_error_correction with the intermediate config file.
    the intermediate config file is copy of the original config file 
    with the input file path modified to the output of the NOREC4DNA encoder.
    (in this case : ./DNA-AEON/data/{filename}_RU10.zip)
    once done is destroy the .zip file

    :param file:
    :return:
    """
    filename = config["encode"]["input"].split("/")[-1]
    if __debug__:
        pass
    config["encode"]["input"] = "{cpath}/data/{filename}_RU10.zip".format(cpath=current_path, filename=filename)
    with open("intermediate_config.json", "w") as inter:
        json.dump(config, inter)
    # Inner encoder command
    py_command = ("{cpath}/bin/arithmetic_modulator_error_correction -e {cpath}/{config_file}".format(
        cpath=current_path, config_file="intermediate_config.json"))
    if (mode == "subprocess"):
        process = subprocess.Popen(py_command.split(), stdout=subprocess.PIPE)
    elif (mode == "sys"):
        process = subprocess.Popen(py_command.split(), stdout=sys.stdout, stderr=sys.stderr)
    output, error = process.communicate()
    os.remove("intermediate_config.json")
    # Remove intermediary files (to keep to validate the datas)
    if not config["encode"]["keep_intermediary"]:
        os.remove(config["encode"]["input"])
    else :
        os.rename(config["encode"]["input"], "{cpath}/data/intermediates_results/{filename}_RU10.zip".format(cpath=current_path, filename=filename))
    #ret = parse_fasta(file + "_encoded.fasta")
    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Encode data in DNA using the concatenation of the NOREC4DNA raptor-fountain implementation and DNA-Aeon.')
    parser.add_argument('--config', '-c', dest='conf', type=str, action='store',
                        help='path to the config file.', required=True)
    parser.add_argument('--mode', '-m', dest='mode', type=str, action='store',
                        help='mode of output. "subprocess", "sys", "file"', required=False, default="subprocess")
    args = parser.parse_args()

    cpath = pathlib.Path(__file__).parent.parent.resolve()
    with open(args.conf, "r") as conf_inp:
        config_data = json.load(conf_inp)

    # Start outer encoder
    print("Starting outer encoder")
    encode_norec_for_ac(config_data, cpath) #args.mode) #modifying can cause errors for norec_config
    print("\n\nFinished outer encoding, starting inner encoder...\n")
    # Start inner encoder
    print("Starting inner encoder")
    encode_ac(cpath, config_data, args.mode)
    file_ext = "" if config_data["general"]["as_fasta"] else ".zip"
    output_path = pathlib.Path(config_data["encode"]["output"] + file_ext).resolve()
    print("Finished encoding! data can be found at {output_path}".format(output_path=output_path))
    print("The NOREC4DNA encoder config file can be found at {norec_conf}".format(norec_conf=pathlib.Path(config_data["decode"]["NOREC4DNA_config"]).resolve()))
