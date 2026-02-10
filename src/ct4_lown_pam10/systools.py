"""This module contains different utility making use of linux programs like awk, grep etc."""

import logging
import os
import subprocess
import sys
from pathlib import Path
from typing import Any, TextIO

from numpy import uint32

__module_name__ = "System Tools"
__source__ = "sys_tools.py"
__author__ = "Salvatore Cosentino"
# __copyright__ = ""
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Cosentino Salvatore"
__email__ = "salvo981@gmail.com"


# Logger that will be used in this module
# It is child of the root logger and
# should be initialiazied using the function set_logger()
logger: logging.Logger = logging.getLogger()



def info() -> None:
    """This module contains different utility making use of linux programs like awk, grep etc."""
    print(f"MODULE NAME:\t{__module_name__}")
    print(f"SOURCE FILE NAME:\t{__source__}")
    print(f"MODULE VERSION:\t{__version__}")
    print(f"LICENSE:\t{__license__}")
    print(f"AUTHOR:\t{__author__}")
    print(f"EMAIL:\t{__email__}")



def bzip2(path: str, outdir: str = os.getcwd(), level: int = 9, overwrite: bool = False, keep: bool = True) -> str:
    """Use bzip2 to compress the input archive to the output directory."""
    import bz2

    # debugStr: str = f"""bzip2 :: START
    # Input file:\t{path}
    # Output dir:\t{outdir}
    # Compression level:\t{level:d}
    # Overwrite existing compressed file:\t{overwrite}
    # Keep original file:\t{keep}
    # """
    # logger.debug(debugStr)

    logger.debug(f"""bzip2 :: START
    Input file:\t{path}
    Output dir:\t{outdir}
    Compression level:\t{level:d}
    Overwrite existing compressed file:\t{overwrite}
    Keep original file:\t{keep}
    """)

    # check if the zip file is valid
    if not os.path.isfile(path):
        sys.stderr.write(f"\nERROR: {path} is not a valid file.")
        sys.exit(-2)
    # create the output directory
    makedir(outdir)
    # create the new path
    newname: str = f"{os.path.basename(path)}.bz2"
    outpath: str = os.path.join(outdir, newname)
    # check that the output file does not already exist
    if os.path.isfile(outpath):
        if not overwrite:
            sys.stderr.write(
                f"\nERROR: the archive {outpath} already exists.\nSet overwrite=True to overwrite."
            )
            sys.exit(-2)
    # open output file
    # note that this would completely load the input file in memory
    new_file = bz2.BZ2File(outpath, mode="wb", compresslevel=level)
    # now write the data
    new_file.write(bz2.compress(open(path, "rb").read(), compresslevel=9))
    logger.debug(f"\nThe archive\n{path}\nwas compressed to\n{outpath}")
    # keep/remove original
    if not keep:
        os.remove(path)
        logger.debug(f"\nThe original raw file \n{path}\n has been removed.")
    # return the path to the extracted archive
    return outpath



def copy(src: str, dst: str, metadata: bool = False) -> bool:
    """Copy src file/dir to dst."""

    logger.debug(f"""copy :: START
    Src file: {src}
    Dest file: {dst}
    Metadata:\t{metadata}
    """)

    # check the existence of the input file
    if not os.path.isfile(src):
        sys.stderr.write(f"The file {src} was not found, please provide a valid file path")
        sys.exit(-2)
    # if src and dst are same, do nothing...
    if src == dst:
        sys.stderr.write("\nWARNING: Source and destination files are the same, nothing will be done.\n")
        return False
    import shutil

    # let's execute commands
    if metadata:  # then also copy the metadata
        try:
            shutil.copy2(src, dst)
        # eg. src and dest are the same file
        except shutil.Error as e:
            print(f"Error: {e}")
        # eg. source or destination doesn't exist
        except IOError as e:
            print(f"Error: {e.strerror}")
    else:
        try:
            shutil.copy(src, dst)
        # eg. src and dest are the same file
        except shutil.Error as e:
            print(f"shutil.Error: {e}")
        # eg. source or destination doesn't exist
        except IOError as e:
            print(f"IOError: {e.strerror}")
    return True



def count_lines_wc(infile: str) -> uint32:
    """Takes in input a text file and uses WC to count the number of lines."""
    logger.debug(f"""count_line_wc :: START
    Input file: {infile}
    """)

    # check the existence of the input file
    if not os.path.isfile(infile):
        sys.stderr.write(f"The file {infile} was not found, please provide a input path")
        sys.exit(-2)

    # EXAMPLES of command to be executed
    # wc -l <infile.fna>
    cmd: str = f"wc -l {infile}"
    #execute the system call
    completed: subprocess.CompletedProcess = subprocess.run(cmd, shell=True, capture_output=True)
    if completed.returncode != 0:
        logger.error(f"wc command failed.\nExit code:\t{completed.returncode}\nSTDERR: {completed.stderr}")
        sys.exit(completed.returncode)

    inlines: uint32 = uint32(completed.stdout.split(b" ")[0])

    logger.debug(f"""Count lines CMD:\t{cmd}
    Counted lines:\t{inlines:d}
    """)

    return inlines



def chop_string(s: str, n: int):
    """
    Chop strings to a given size.
    Produce (yield) \'n\'-character chunks from \'s\'.
    """

    logger.debug(f"""chop_string :: START
    Input: {s}
    Chunk size:\t{n:d}
    Input length:\t{len(s):d}
    """)

    for start in range(0, len(s), n):
        yield s[start : start + n]



def create_flogger(logpath: str, loggername: str, lev: int = 10, mode: str = "a", propagate: bool = False) -> logging.Logger:
    """Create a logger that writes into a file"""
    # THIS CREATES THE GENERAL LOGGER
    logger: logging.Logger = logging.getLogger(loggername)
    # set same level as the root
    logger.setLevel(lev)
    logger.propagate = propagate
    logfh: logging.FileHandler = logging.FileHandler(logpath, mode=mode)
    # This makes sure that the log file is created even if not in debug mode
    logfh.setLevel(lev)
    logfh.setFormatter(fmt=logging.Formatter("%(message)s"))
    logger.addHandler(logfh)

    return logger



def get_cpu_count() -> int:
    """Get the number of cpu available in the system."""
    from multiprocessing import cpu_count
    logger.debug("get_cpu_count :: START")
    return cpu_count()



def get_sys_info() -> dict[str, Any]:
    """Obtain system information."""
    logger.debug("get_sys_info :: START")
    sysdict: dict[str, Any] = {}
    # Indentify the operative system
    from platform import uname

    sysdict["os"] = uname().system
    # Check if it is a MacOS
    if sysdict["os"] == "Darwin":
        sysdict["is_darwin"] = True
    else:
        sysdict["is_darwin"] = False
    # Cpu architecture
    sysdict["architecture"] = uname().machine
    from psutil import virtual_memory

    # now compute the memory per thread
    availMem: float = round(virtual_memory().total / 1073741824.0, 2)
    sysdict["mem"] = str(availMem)
    # CPU count
    sysdict["cpu"] = str(get_cpu_count())
    # Get info about the python installation
    sysdict["py_ver"] = sys.version
    # Python bin path
    sysdict["py_path"] = sys.executable

    return sysdict



def load_input_paths(dir: Path, fmt: str = "fasta") -> list[Path]:
    """Load the input paths of files in a specified directory."""

    logger.debug(f"""load_input_paths :: START
    Directory: {dir}
    Load files with format: {fmt}
    """)

    # Check input dir and requested format
    if not dir.is_dir():
        logger.error(f"The provided path is not a valid directory:\n{dir}")
        sys.exit(-2)
    validfmts: list[str] = ["fasta", "fna", "faa", "tsv"]
    if fmt not in validfmts:
        logger.error(f"{fmt} is not a valid format.")
        sys.exit(-10)

    tmppath: Path = Path()
    fname: str = ""
    paths: list[Path] = []
    fnames: list[str] = os.listdir(dir)
    fnames.sort()
    for fname in fnames:
        if fname.endswith(f".{fmt}"):
            tmppath = dir.joinpath(fname)
            if tmppath.is_file():
                paths.append(tmppath)
            else:
                logger.error(f"{fname} is not a valid {fmt} file.\n{tmppath}")
                sys.exit(-2)

    logger.debug(f"{fmt} files loaded:\t{len(paths)}")

    return paths



def makedir(path: Path) -> None:
    """Create a directory including the intermediate directories in the path if not existing."""
    # check the file or dir does not already exist
    if path.is_file():
        sys.stderr.write(f"\nWARNING: {path}\nalready exists as a file, and the directory cannot be created.\n")
    try:
        os.makedirs(path)
    except OSError:
        if not path.is_dir():
            raise



def move(src: str, dst: str, debug: bool = False) -> None:
    """Recursively moves src to dst."""
    if debug:
        print("move :: START")
        print(f"SRC:\n{src}")
        print(f"DEST:\n{dst}")
    # check the existence of the input file
    if not os.path.exists(src):
        sys.stderr.write(f"{src} was not found, please provide a valid path")
        sys.exit(-2)
    import shutil

    # let's execute command
    if os.path.exists(dst):  # then we should use copy
        # copy and remove the source file
        copy(src, dst, True)
        os.remove(src)
    else:
        shutil.move(src, dst)



def set_logger(
    logger_name: str,
    lev: int,
    propagate: bool,
    custom_fmt: logging.Formatter = logging.Formatter(fmt=None),
) -> None:
    """Set the global logger for this module"""
    global logger
    logger = logging.getLogger(logger_name)
    logger.setLevel(lev)
    logger.propagate = propagate
    # Create the handler and
    clsLogger: logging.StreamHandler = logging.StreamHandler(stream=sys.stdout)
    # This makes sure that the log file is created even if not in debug mode
    clsLogger.setLevel(logger.level)
    # Set the formatter
    if custom_fmt is not None:
        clsLogger.setFormatter(custom_fmt)
    logger.addHandler(clsLogger)
    # write some log about it!
    logger.debug(f"General logger for {logger_name} loaded!")
