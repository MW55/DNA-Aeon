from subprocess import run
import pathlib
from venv import create
from os.path import abspath
import sys

#todo: activate an option to just run cmake and just compile the project


def create_norec_env(current_path):
    venv_dir = "{cpath}/NOREC4DNA/venv".format(cpath=current_path)
    create(venv_dir, with_pip=True)

def install_norec_packages(current_path):
    run(["venv/bin/pip", "install", "wheel"], cwd="{cpath}/NOREC4DNA".format(cpath=current_path))
    run(["venv/bin/pip", "install", "-r", abspath("./NOREC_requirements.txt")], cwd="{cpath}/NOREC4DNA".format(cpath=current_path)) #"./NOREC_requirements.txt"

def compile_dna_aeon(current_path, debug_mode):
    if debug_mode == True :
        run(["cmake", "-DENABLE_DEBUG_MACRO=ON", current_path])
    else :
        run(["cmake", current_path])
    run(["cmake", "--build", current_path]) 

# Main function should offers debug options for cmake
if __name__ == "__main__":
    args = sys.argv[1:]
    if args[0] == "debug":
        debug_mode = True
        print("Debug mode enabled")
    else:
        debug_mode = False
        print("Release mode enabled")
    cpath = pathlib.Path(__file__).parent.resolve()
    if len(args) > 1:
        if args[1] == "cmake":
            compile_dna_aeon(cpath, debug_mode)
        else:
            print("Setting up NOREC4DNA virtual environment.")
            create_norec_env(cpath)
            print("Installing packages required for NOREC4DNA.")
            install_norec_packages(cpath)
            print("Compiling DNA Aeon")
            compile_dna_aeon(cpath, debug_mode)
            print("Installation finished!")

