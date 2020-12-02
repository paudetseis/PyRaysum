import os
import subprocess
import pyraysum
from pathlib import Path
from pyraysum import src
from argparse import ArgumentParser


def get_install_args(argv=None):

    parser = ArgumentParser(
        usage="%(prog)s [options]",
        description="Script used to compile the Raysum Fortran modules " +
        "and install them locally on the system. The script will attempt " +
        "to use default paths, or you can specify the path where you want " +
        "the binaries to be installed. You can also specify the Fortran " +
        "compiler of your choice that exists on your system. If PyRaysum " +
        "was installed using conda with the conda-provided fortran " +
        "compiler, you can safely ignore these options. However, specifying " +
        "these options will take precedence over the conda defaults.")
    parser.add_argument(
        "--path",
        action="store",
        type=str,
        dest="install_path",
        default="",
        help="Destination where the Fortran binaries will be installed. " +
        "Use absolute path only - the install may fail otherwise. Check " +
        "your $PATH environment variable to find out what these are. " +
        "If your path is in the root directory, you may have to run this " +
        "script with super user privileges.")
    parser.add_argument(
        "--fcompiler",
        action="store",
        type=str,
        dest="fcompiler",
        default="gfortran",
        help="Fortran compiler used to produce binaries for Raysum. " +
        "For example, --fcompiler=gfortran or --fcompiler=ifort. Check " +
        "your $FC environment variable to find out if you have a Fortran " +
        "compiler installed.")

    args = parser.parse_args(argv)

    return args


def main(args=None):

    if args is None:
        args = get_install_args()

    conda_exists = os.getenv("CONDA_PREFIX")

    # If user specifies path, use it
    if args.install_path:
        prs_path = Path(args.install_path)

    # Otherwise, attempt to find suitable path
    else:
        # First find out if conda is installed
        if conda_exists:
            prs_path = Path(os.getenv("CONDA_PREFIX")) / "bin"
        # Conda not installed - use first path in $PATH environment
        else:
            prs_path = Path(os.getenv("PATH").split(":")[0])

    # If user specifies fortran compiler, use it
    if args.fcompiler:
        os.environ["FC"] = args.fcompiler

    # Print message to screen where binaries are installed
    print()
    print("*"*len(str(prs_path)))
    print()
    print("Fortran compiler is : \n"+str(os.getenv("FC")))
    print()
    print("Raysum will be installed in : \n"+str(prs_path))
    print()
    print("*"*len(str(prs_path)))
    print()

    # Change directory to where the Fortran source files are located
    src_path = src.__path__[0]
    os.chdir(src_path)

    # Set FC environment variable

    # Compile Raysum using the Makefile: clean, compile, copy binaries, clean
    subprocess.call(["make", "clean"])
    subprocess.call(["make", "all"])
    subprocess.call(["cp", "-f", "seis-spread", str(prs_path)])
    subprocess.call(["make", "clean"])


if __name__ == "__main__":

    # Run main program
    main()
