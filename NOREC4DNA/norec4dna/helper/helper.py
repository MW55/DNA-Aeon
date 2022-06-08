import math
import os
import shutil

mode = "single_cpu"

if mode == "gpu":
    from norec4dna.helper.helper_cuda import *
elif mode == "gpu_simple":
    from norec4dna.helper.helper_cuda_simple import *
elif mode == "single_cpu":
    from norec4dna.helper.helper_cpu_single_core import *
else:
    from norec4dna.helper.helper_cpu import *


def split_file(in_file_name, number_of_splits):
    filesize = os.path.getsize(in_file_name)
    dirs = 'split_{}'.format(in_file_name)
    try:
        os.makedirs(dirs)
    except:
        for f in [f for f in os.listdir(dirs)]:
            os.remove(os.path.join(dirs, f))
    out_file_names = []
    float_chunk_size = filesize / number_of_splits
    chunk_size = math.ceil(float_chunk_size)
    if float_chunk_size != chunk_size:
        print("File does not perfectly split into {} Chunks. Last Chunk will be smaller!".format(number_of_splits))
    with open(in_file_name, 'rb') as in_file:
        for i in range(number_of_splits):
            out_file_name = dirs + '/{}.{}'.format(in_file_name, i)
            out_file_names.append(out_file_name)
            with open(out_file_name, 'wb') as out_file:
                tmp = in_file.read(chunk_size)
                out_file.write(tmp)
    return out_file_names


def find_ceil_power_of_four(n):
    res = 0
    while True:
        if n <= math.pow(4, res):
            return res
        else:
            res += 1


def number_to_base_str(number, target_str_length):
    digs = "ACGT"

    def int2base(x, base):
        if x < 0:
            sign = -1
        elif x == 0:
            return digs[0]
        else:
            sign = 1

        x *= sign
        digits = []

        while x:
            digits.append(digs[int(x % base)])
            x = int(x / base)

        if sign < 0:
            digits.append('-')

        digits.reverse()

        return ''.join(digits)

    res = int2base(number, 4)
    res = res.rjust(target_str_length, 'A')
    print(base_str_to_int(res))
    return res


def base_str_to_int(base_str):
    mapping_dict = {'A': "0", 'C': "1", 'G': "2", 'T': "3"}
    return int("".join([mapping_dict[base] for base in base_str]), 4)


def cluster_and_remove_index(split_index_position: str, split_index_length: int, folder: str) -> typing.Tuple[
    typing.List[str], str]:
    number_folder_mapping = dict()
    if os.path.exists("cluster_out"):
        shutil.rmtree('cluster_out')
    for file in os.listdir(folder):
        file = os.path.join(folder, file)
        if file.endswith("DNA") and not os.path.isdir(file):
            with open(file, "r") as f:
                content = f.read()
            if split_index_position == "start":
                bin_number, base_str = base_str_to_int(content[:split_index_length]), content[split_index_length:]
            else:
                bin_number, base_str = base_str_to_int(content[-split_index_length:]), content[:-split_index_length]
            if bin_number not in number_folder_mapping:
                created_folder = os.path.join("cluster_out", str(bin_number))
                os.makedirs(created_folder)
                number_folder_mapping[bin_number] = created_folder
                # os.mkdir(os.path.join(number_folder_mapping[bin_number]))
            with open(os.path.join(number_folder_mapping[bin_number], os.path.basename(file)), "w") as out_f:
                out_f.write(base_str)
    return [x for x in number_folder_mapping.values()], number_folder_mapping[max(number_folder_mapping.keys())]


def fasta_cluster_and_remove_index(split_index_position: str, split_index_length: int, file: str) -> typing.Tuple[
    typing.List[str], str]:
    number_file_mapping = dict()
    if os.path.exists("cluster_out"):
        shutil.rmtree('cluster_out')
    os.makedirs("cluster_out")
    with open(file, "r") as in_file:
        while True:
            first_line = in_file.readline()
            if not first_line:
                break
            line = in_file.readline()
            if not line:
                break
            dna_str = line.replace("\n", "")
            if split_index_position == "start":
                bin_number, base_str = base_str_to_int(dna_str[:split_index_length]), dna_str[split_index_length:]
            else:
                bin_number, base_str = base_str_to_int(dna_str[-split_index_length:]), dna_str[:-split_index_length]
            if bin_number not in number_file_mapping:
                bin_out_fasta_file = os.path.join("cluster_out", str(bin_number) + ".fasta")
                number_file_mapping[bin_number] = bin_out_fasta_file
            with open(number_file_mapping[bin_number], "a+") as out_f:
                out_f.write(first_line + base_str + "\n")
    return [x for x in number_file_mapping.values()], number_file_mapping[max(number_file_mapping.keys())]


def merge_folder_content(src_folder_of_folders, dest_folder, append_folder_name=True, clear_dest_folder=False):
    print(src_folder_of_folders)
    print(dest_folder)
    if clear_dest_folder and os.path.exists(dest_folder):
        try:
            shutil.rmtree(dest_folder)
        except:
            print("Could not delete folder")
    try:
        os.mkdir(dest_folder)
    except Exception:
        print("Folder already exists.")
    if len(os.listdir(dest_folder)) > 0:
        raise FileExistsError("dest_folder was not empty!")

    for folder in os.listdir(src_folder_of_folders):
        folder = os.path.join(os.path.abspath(src_folder_of_folders), folder)
        if os.path.isdir(folder):
            for file in os.listdir(folder):
                if os.path.isfile(os.path.join(folder, file)):
                    if append_folder_name:
                        dest_file = os.path.join(dest_folder, os.path.basename(folder) + "_" + file)
                    else:
                        dest_file = dest_folder + "_" + file
                    shutil.copy(os.path.join(folder, file), dest_file)


def split_first(x):
    split_text = x.rsplit('.', 1)
    if split_text[1].lower() == "fasta":
        return split_text[0].split("_")[1]
    else:
        return split_text[1]


def merge_parts(filenames, remove_tmp_on_success=False):
    numbers = [int(split_first(x)) for x in filenames]
    max_num = max(numbers)
    assert len(filenames) == max_num + 1, "Some parts were not decoded, try manual merge."
    base_name = filenames[0].rsplit('.', 1)[0]
    try:
        os.remove(base_name)
    except:
        print("Error while removing file: {}, new content will be appended!".format(base_name))
    with open(base_name, "ab") as out_f:
        for i in range(max_num + 1):
            with open(base_name + "." + str(i), "rb") as in_f:
                while True:
                    buffer = in_f.read(65536)
                    if buffer:
                        out_f.write(buffer)
                    else:
                        break
    print("Final file saved as {}".format(base_name))
    if remove_tmp_on_success:
        for i in range(max_num + 1):
            os.remove(base_name + "." + str(i))


if __name__ == "__main__":
    print(os.listdir(os.path.curdir))
    merge_folder_content("../../split_Dorn", "../../Dorn_combined_split_output", True, True)
