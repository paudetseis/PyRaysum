import os
import subprocess
import pyraysum
from pathlib import Path
from pyraysum import src


def main():

    # Find out if conda is installed
    conda_exists = os.getenv("CONDA_PREFIX")

    # Check out where binaries will be located
    if conda_exists:
        prs_path = Path(os.getenv("CONDA_PREFIX")) / "bin"
    else:
        prs_path = Path(os.getenv("PATH").split(":")[0])

    # Print message to screen where binaries are installed
    print()
    print("*"*len(str(prs_path)))
    print()
    print("Raysum will be installed in : \n"+str(prs_path))
    print()
    print("*"*len(str(prs_path)))
    print()

    # Change directory to where the Fortran source files are located
    src_path = src.__path__[0]
    os.chdir(src_path)

    # Compile Raysum using the Makefile: clean, compile, copy binaries, clean
    subprocess.call(["make", "clean"])
    subprocess.call(["make", "all"])
    subprocess.call(["cp", "-f", "seis-spread", str(prs_path)])
    subprocess.call(["make", "clean"])


if __name__ == "__main__":

    # Run main program
    main()
