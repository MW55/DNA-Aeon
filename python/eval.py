
import pathlib

import numpy as np
import subprocess
import random
from copy import deepcopy
import os
import shutil
import filecmp
#import yaml
from time import time
import json


def run_eval(ac_path, norec_path, ini_name, out_path, orig_path, out_zip_name, length, norec_out_path, conf_file_path):
    subprocess.call([ac_path, '-d', conf_file_path])
    print(out_zip_name)
    if(os.path.exists(out_zip_name)) :
        shutil.copy(out_zip_name, norec_path + '/' + out_zip_name)
        try:
            subprocess.call([norec_path + '/venv/bin/python', norec_path + '/ConfigWorker.py', norec_path + "/" + ini_name])
            return filecmp.cmp(orig_path, norec_out_path)
        except:
            return False
    else:
        return False


def auto_test(out_path, ac_path, norec_path, ini_name, orig_path, out_zip_name, length, iterations, results, counter, conf_file_path):
    for i in range(iterations):
        results[str(counter)] = dict()
        res = run_eval(ac_path, norec_path, ini_name, out_path, orig_path, out_zip_name, length, results, conf_file_path)
        results[str(counter)]["decoded"] = res
        counter += 1
        print(counter)

    return results

if __name__ == '__main__':
    """
    this script is used to evaluate the error correction of the arithmetic modulator
    and the NOREC4DNA tool. It uses the arithmetic modulator to encode a fasta file
    and then decodes it with the NOREC4DNA tool. The decoded file is then compared to
    the original file. The script is used to evaluate the error correction of the
    arithmetic modulator.
    
    1) out_path: the path to the encoded fasta file
    2) ac_path: the path to the arithmetic modulator executable
    3) norec_path: the path to the NOREC4DNA tool
    4) ini_name: the name of the ini file used by the NOREC4DNA tool
    5) orig_path: the path to the original fasta file
    6) norec_out_path: the path to the decoded file
    7) out_zip_name: the name of the decoded file (concatenate from norec_path and out_zip_name)
    8) conf_file_path: the path to the configuration file
    9) length: the length of the encoded sequence
    10) iterations: the number of iterations to run
    11) results: the results of the evaluation
    
    """
    project_path = str(pathlib.Path(__file__).parent.parent.resolve())
    out_path = project_path + "/data/dorn_encoded_mod.fasta"
    ac_path = project_path + "/bin/arithmetic_modulator_error_correction"
    norec_path = project_path + "/libraries/NOREC4DNA"
    ini_name = "Dorn.ini"
    #orig_path = project_path + "/arithmetic_modulator_error_correction/data/Dorn"
    orig_path = project_path + "/data/Dorn"
    norec_out_path = project_path + "/data/results/eval/Dorn"
    #out_zip_name = "dorn_decoded_auto.zip" #what is this?
    out_zip_name = "Dorn.zip"
    conf_file_path = project_path + "/configs/config-files/eval_config.json"
    length = 172
    iterations = 1
    results = dict()
    counter = 0
    config = {
        "general": {
            "sync": 4,
            "as_fasta": True
        },
        "decode": {
            "input": out_path,
            "output": out_zip_name,
            "length": length
        }
    }
    with open(conf_file_path, "w") as conf_file:
        json.dump(config, conf_file)
    if os.path.exists(norec_out_path):
        os.remove(norec_out_path)
    print("start")
    start = time()
    auto_test(out_path, ac_path, norec_path, ini_name, orig_path, out_zip_name, length, iterations, results, counter, conf_file_path)
    end = time()
    print(f"time : {(end - start)*1000:.2f} ms")
    print("done")
    #with open("test_split.yaml", "w") as o_:
    #    mod_data = yaml.load(o_, Loader=yaml.CLoader)
    #    re = z = mod_data | results
    #    yaml.dump(re, o_)
