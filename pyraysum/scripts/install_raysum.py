import os
import subprocess
from pathlib import Path
from pyraysum import src


def main():

    # Check out where binaries will be located
    prs_path = Path(os.getenv("CONDA_PREFIX")) / "bin"

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
    subprocess.call(["cp", "seis-spread", prs_path.name])
    subprocess.call(["make", "clean"])


if __name__ == "__main__":

    # Run main program
    main()
