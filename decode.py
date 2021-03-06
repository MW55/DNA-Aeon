import subprocess
import argparse
import pathlib
import json
import os

def decode_ac(current_path, config_path):
    py_command = ("{cpath}/arithmetic_modulator_error_correction -d {conf_path}".format(cpath=current_path,
                                                                                                conf_path=config_path))
    process = subprocess.Popen(py_command.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    return

def decode_norec_for_ac(current_path, config_data):
    """
    outer decoder for the fountaincode-arithmetic code concatenation.
    :param filename:
    :return:
    """
    ini_path = pathlib.Path(config_data["decode"]["NOREC4DNA_config"]).resolve()
    with open(ini_path, "r") as c_:
        line_list = c_.readlines()
        line_list[0] = "[{cpath}/{dpath}.zip]\n".format(cpath=current_path, dpath=config_data["decode"]["output"]) #"[../data/" + filename + "_RU10.zip]\n"
        with open("{cpath}/{ini_path}".format(cpath=current_path, ini_path=config_data["decode"]["NOREC4DNA_config"]), "w") as o_:
            o_.writelines(line_list)
    pathlib.Path('data/results').mkdir(parents=True, exist_ok=True)
    os.chdir('data/results')
    py_command = ("{cpath}/NOREC4DNA/venv/bin/python3 {cpath}/NOREC4DNA/ConfigWorker.py {ipath}".format(
        cpath=current_path, ipath=ini_path))
    process = subprocess.Popen(py_command.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    os.chdir('../..')
    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Decode data that was encoded in DNA using the concatenation of the NOREC4DNA raptor-fountain implementation and DNA-Aeon.')
    parser.add_argument('--config', '-c', dest='conf', type=str, action='store',
                        help='path to the config file.', required=True)
    args = parser.parse_args()
    cpath = pathlib.Path(__file__).parent.resolve()
    conf_path = pathlib.Path(args.conf).resolve()
    with open(args.conf, "r") as conf_inp:
        config_data = json.load(conf_inp)

    # Start outer decode
    print("Starting inner decoder.")
    decode_ac(cpath, conf_path)
    print("\nFinished inner decoding, starting outer decoder...\n")
    # Start inner encoder
    decode_norec_for_ac(cpath, config_data)
    print("\nFinished outer decoding. File can be found at {cur_path}/data/results/{inp_file}".format(
        cur_path=cpath, inp_file=config_data["encode"]["input"].split("/")[-1]))
    '''
    file_ext = "" if config_data["general"]["as_fasta"] else ".zip"
    output_path = pathlib.Path(config_data["encode"]["output"] + file_ext).resolve()
    print("Finished encoding! data can be found at {output_path}".format(output_path=output_path))
    print("The NOREC4DNA encoder config file can be found at {norec_conf}".format(norec_conf=pathlib.Path(config_data["decode"]["NOREC4DNA_config"]).resolve()))
    '''