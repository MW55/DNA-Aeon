import subprocess
import argparse
import pathlib
import json
import os
import sys


# This function should printout the C++ output but it does not.
# 1) solution use sys.stdout and sys.stderr as arguments for subprocess.Popen
# @todo: I want a to print the output of the C++ program to a file.


def decode_ac(current_path, config_path, mode="subprocess"):
    
    py_command = ("{cpath}/bin/arithmetic_modulator_error_correction -d {conf_path}".format(cpath=current_path,                                                                                       conf_path=config_path))
    if (mode == "subprocess"):
        process = subprocess.Popen(py_command.split(), stdout=subprocess.PIPE)
    elif (mode == "sys"):
        process = subprocess.Popen(py_command.split(), stdout=sys.stdout, stderr=sys.stderr)
    elif (mode == "file"):
        process = subprocess.Popen(py_command.split(), stdout=subprocess.PIPE)
        with open("configs/debug/debug_output.txt", "w") as write_file:
            write_file.write(str(process.stdout.read()))
    else:
        print("Invalid mode. Exiting...")
        exit(1)
    output, error = process.communicate()
    return

def decode_norec_for_ac(current_path, config_data, mode="subprocess"):
    """
    outer decoder for the fountaincode-arithmetic code concatenation.
    :param filename:
    :return:
    :brief: 
    """
    ini_path = pathlib.Path(config_data["decode"]["NOREC4DNA_config"]).resolve()
    # print(ini_path)
    with open(ini_path, "r") as c_:
        line_list = c_.readlines()
        line_list[0] = "[{cpath}/{dpath}.zip]\n".format(cpath=current_path, dpath=config_data["decode"]["output"]) #"[../data/" + filename + "_RU10.zip]\n"
        with open("{cpath}/{ini_path}".format(cpath=current_path, ini_path=config_data["decode"]["NOREC4DNA_config"]), "w") as o_:
            o_.writelines(line_list)
    #pathlib.Path('data/results').mkdir(parents=True, exist_ok=True)
    #os.chdir('data/results')
    pathlib.Path('data').mkdir(parents=True, exist_ok=True)
    os.chdir('data')
    py_command = ("{cpath}/libraries/NOREC4DNA/venv/bin/python3 {cpath}/libraries/NOREC4DNA/ConfigWorker.py {ipath}".format(
        cpath=current_path, ipath=ini_path))
    if (mode == "subprocess"):
        process = subprocess.Popen(py_command.split(), stdout=subprocess.PIPE)
    elif (mode == "sys"):
        process = subprocess.Popen(py_command.split(), stdout=sys.stdout, stderr=sys.stderr)
    else :
        print("Invalid mode. Exiting...")
        exit(1)
    output, error = process.communicate()
    os.chdir('../..')
    return

# main should have an extra argument to specify the mode of the decode process "subprocess", "sys", "file"
# main should have an extra argument to specify debug mode of the encode process
if __name__ == "__main__":
    if __debug__: print("Debug mode is on. Asserts will be checked.")
    else: print("Debug mode is off. Asserts will not be checked.")
    parser = argparse.ArgumentParser(description='Decode data that was encoded in DNA using the concatenation of the NOREC4DNA raptor-fountain implementation and DNA-Aeon.')
    parser.add_argument('--config', '-c', dest='conf', type=str, action='store',
                        help='path to the config file.', required=True)
    parser.add_argument('--mode', '-m', dest='mode', type=str, action='store',
                        help='mode of output. "subprocess", "sys", "file"', required=False, default="subprocess")
    parser.add_argument('--skipinner', '-i', dest='skipinner', action='store_true', required=False, default=False)
    parser.add_argument('--skipouter, -o', dest='skipouter', action='store_true', required=False, default=False)
    args = parser.parse_args()
    cpath = pathlib.Path(__file__).parent.parent.resolve()
    conf_path = pathlib.Path(args.conf).resolve()
    with open(args.conf, "r") as conf_inp:
        config_data = json.load(conf_inp)
    if not args.skipinner:
        print("Starting inner decoder.")
        #args mode has been removed (observe if bug occurs)
        decode_ac(cpath, conf_path, args.mode)
        print("\nFinished inner decoding...\n")
    else:
        print("Skipping inner decoding with (-s). Starting outer decoding...\n")
    if not args.skipouter:
        print("Starting outer decoding.\n")
        #args mode has been removed (observe if bug occurs)
        decode_norec_for_ac(cpath, config_data, args.mode)
        print("\nFinished outer decoding. File can be found at {cur_path}/data/results/{inp_file}".format(
        cur_path=cpath, inp_file=config_data["encode"]["input"].split("/")[-1]))
    else:
        print("Skipping outer decoding with (-o). Exiting...")
    '''
    file_ext = "" if config_data["general"]["as_fasta"] else ".zip"
    output_path = pathlib.Path(config_data["encode"]["output"] + file_ext).resolve()
    print("Finished encoding! data can be found at {output_path}".format(output_path=output_path))
    print("The NOREC4DNA encoder config file can be found at {norec_conf}".format(norec_conf=pathlib.Path(config_data["decode"]["NOREC4DNA_config"]).resolve()))
    '''