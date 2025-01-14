from subprocess import run
import pathlib
from venv import create
from os.path import abspath


def create_norec_env(current_path, python_executable="python3.9"):
    venv_dir = f"{current_path}/NOREC4DNA/venv"
    run([python_executable, "-m", "venv", venv_dir])


def install_norec_packages(current_path):
    run(["venv/bin/pip", "install", "wheel"], cwd="{cpath}/NOREC4DNA".format(cpath=current_path))
    run(["venv/bin/pip", "install", "-r", abspath("./NOREC_requirements.txt")], cwd="{cpath}/NOREC4DNA".format(cpath=current_path)) #"./NOREC_requirements.txt"

def compile_dna_aeon(current_path):
    run(["cmake", current_path])
    run(["cmake", "--build", current_path]) 

if __name__ == "__main__":
    cpath = pathlib.Path(__file__).parent.resolve()
    print("Setting up NOREC4DNA virtual environment.")
    create_norec_env(cpath)
    print("Installing packages required for NOREC4DNA.")
    install_norec_packages(cpath)
    print("Compiling DNA Aeon")
    compile_dna_aeon(cpath)
    print("Installation finished!")

