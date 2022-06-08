import numpy as np
import subprocess
import random
from copy import deepcopy
import os
import shutil
import filecmp
import yaml
from time import time
import json


def run_eval(ac_path, norec_path, ini_name, out_path, orig_path, out_zip_name, length, norec_out_path, conf_file_path):
    subprocess.call([ac_path, '-d', conf_file_path])
    shutil.copy(out_zip_name, norec_path + '/' + out_zip_name)
    try:
        subprocess.call([norec_path + '/venv/bin/python', norec_path + '/ConfigWorker.py', norec_path + "/" + ini_name])
        return filecmp.cmp(orig_path, norec_out_path)
    except:
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

    out_path = "/home/wintermute/projects/arithmetic_modulator_error_correction/data/dorn_encoded_mod.fasta"
    ac_path = "/home/wintermute/projects/arithmetic_modulator_error_correction/cmake-build-release/arithmetic_modulator_error_correction"
    norec_path = "/home/wintermute/projects/NOREC4DNA"
    ini_name = "Dorn.ini"
    orig_path = "/home/wintermute/projects/arithmetic_modulator_error_correction/data/Dorn"
    norec_out_path = "/home/wintermute/projects/arithmetic_modulator_error_correction/Dorn"
    out_zip_name = "dorn_decoded_auto.zip"
    conf_file_path = "eval_config.json"
    length = 172
    iterations = 1
    results = dict()
    counter = 0
    config = {
        "general": {
            "sync": 4,
            "as_fasta": true
        },
        "decode": {
            "input": out_path,
            "output": out_zip_name,
            "length": length
    }
    with open(conf_file_path, "w") as conf_file:
        json.dump(config, conf_file)
    if os.path.exists(norec_out_path):
        os.remove(norec_out_path)
    print("start")
    start = time()
    auto_test(out_path, ac_path, norec_path, ini_name, orig_path, out_zip_name, length, iterations, results, counter, conf_file_path)
    end = time()
    print(end - start)
    #with open("test_split.yaml", "w") as o_:
    #    mod_data = yaml.load(o_, Loader=yaml.CLoader)
    #    re = z = mod_data | results
    #    yaml.dump(re, o_)
