"""Annotate input proteomes using Kofamscan.
"""
import logging
import os
import sys
from pathlib import Path
import pandas as pd

from numpy import float64, uint32
from rich.logging import RichHandler
from typing_extensions import TextIO
import argparse

import genome_data_analysis as gda

# Set the logger
logger = logging.getLogger(__name__)



### FUNCTIONS ###
def get_params() -> argparse.Namespace:
    """Parse and analyse command line parameters."""

    # define the parameter list
    parser = argparse.ArgumentParser(description="Annotate input proteomes using kofamscan", usage='%(prog)s -i <PROTEOME_SET> -o <OUTPUT_DIRECTORY>[options]', prog="annotate-proteomes")
    # Mandatory arguments
    parser.add_argument("-i", "--input-proteomes", type=str, required=True, help="Directory containing the proteomes to be analyzed.")
    parser.add_argument("-o", "--output-dir", type=str, required=True, help="The directory in which the output will be stored.", default="")
    parser.add_argument("-k", "--kofamscan-path", type=str, required=True, help="Path to the main Kofamscan script (exec_annotate).", default="")
    # General run options
    parser.add_argument("-p", "--profiles", type=str, required=True, help="Path to the file directory, or profile list file for Kofamscan.", default="")
    parser.add_argument("-ko", "--ko-list", type=str, required=True, help="Path to the file file with a list of KO ids for Kofamscan.", default="")
    parser.add_argument("-e", "--file-extension", type=str, required=False, help="Extension of the input FASTA files containing the proteins sequences. Default=faa", default="faa")
    parser.add_argument("-t", "--threads", type=int, required=False, help="Maximum number of CPUs to be used. Default=8", default=8)
    parser.add_argument("-wl", "--write-logs", required=False, help="Write logs to file.", default=False, action="store_true")
    parser.add_argument("-d", "--debug", required=False, help="Show debug lines.", default=False, action="store_true")

    # parse the arguments
    args = parser.parse_args()

    return args




def extract_brite_hierarchies(annotated_file_paths: list[Path], kegg_raw_dir: Path, kegg_json_dir: Path, outdir: Path, kofam_fmt: str = "detail-tsv", best_hit_only: bool = True, exclude_uncomplete_brite: bool = True, overwrite_kegg_obj: bool = False, waittime: float64 = float64(0.5), concatenate: bool=True):
    """Create a single file containing all BRITE hierachies for K-matches from KofamScan for all species"""

    dbgstr: str = f"""extract_brite_hierarchies :: START
    Annotated files:\t{len(annotated_file_paths)}
    Directory with KEGG object: {kegg_raw_dir}
    Directory with KEGG object in JSON: {kegg_json_dir}
    Output directory: {outdir}
    Kofamscan output format:\t{kofam_fmt}
    Only use the first significant hit:\t{best_hit_only}
    Exclude some BRITE fields with poor support:\t{exclude_uncomplete_brite}
    Download KEGG objects:\t{overwrite_kegg_obj}
    Time between requests to KEGG:\t{waittime}
    Concatenate into a single file:\t{concatenate}
    """
    logger.debug(dbgstr)

    prefix: str = ""
    brite_info_files: list[Path] = []
    brite_list_path: Path = Path()

    for p in annotated_file_paths:
        # print(p)
        # base names have the following format:
        # <bacteria>.<id>.kofam.detail-tsv.tsv
        prefix = os.path.basename(p).split(".kofam.", maxsplit=1)[0]
        brite_list_path = gda.annotation.extract_brite_info(kofam_raw=p, kegg_raw_dir=kegg_raw_dir, kegg_json_dir=kegg_json_dir, outdir=outdir, kegg_db_abbrev="ko", kofam_fmt=kofam_fmt, best_hit_only=best_hit_only, exclude_uncomplete_brite=exclude_uncomplete_brite, overwrite_kegg_obj=overwrite_kegg_obj, waittime=waittime)

        # add the path to the list
        brite_info_files.append(brite_list_path)
        # print(prefix)

    # This file would contain all the extracted brite information
    # for all the annotated proteomes
    # including a column for describing the species
    master_brite_file: Path = Path()
    # says if the header has been already written to the file
    has_hdr: bool = False

    # Concatenate all brite info into a single file
    # adding also the species name extracting it from the original file me
    if concatenate:
        master_brite_file = outdir.joinpath("master.brite.list.tsv")

        ofd: TextIO = open(master_brite_file, "wt")
        # Iterate through the files
        for p in brite_info_files:
            prefix = os.path.basename(p).replace(".brite.list.tsv", "")
            with open(p, "rt") as ifd:
                for line in ifd:
                    # print(line)
                    if line.endswith("sequence\n"):
                        if not has_hdr:
                            ofd.write(f"Species\t{line}")
                            has_hdr = True
                        # continue
                    else:
                        ofd.write(f"{prefix}\t{line}")

        ofd.close()

    return master_brite_file



#####  MAIN  #####
def main():
    """Main function executing SonicParanoid"""

    # start setting the needed variables and logger
    args: argparse.Namespace = get_params()
    debug: bool = args.debug
    ## Configure rich logging
    if debug:
        logging.basicConfig(
            level="DEBUG",
            format="%(message)s",
            datefmt="[%X]",
            handlers=[RichHandler()])
    else:
        logging.basicConfig(
            level="INFO",
            format="%(message)s",
            datefmt="[%X]",
            handlers=[RichHandler()])

    # Set main parameters
    proteomes_dir: Path = Path(args.input_proteomes)
    outdir: Path = Path(args.output_dir)
    kofam_path: Path = Path(args.kofamscan_path)
    kofam_profiles: Path = Path(args.profiles)
    ko_list: Path = Path(args.ko_list)
    ext: str = args.file_extension
    threads: uint32 = uint32(args.threads)
    writelog: bool = args.write_logs

    # Show some info
    logger.info(f"""The annotation will be perfomed using the following parameters:
    Proteomes directory: {proteomes_dir}
    Output directory: {outdir}
    Kofamscan program: {kofam_path}
    Kofamscan hmm profiles:\t{kofam_profiles}
    Kofamscan KO list:\t{ko_list}
    FASTA file extension:\t{ext}
    Write logs:\t{writelog}
    CPUs to be used:\t{threads}
    Debug mode:\t{debug}
    """)

    logger.info("\nPerforming proteomes annotation using Kofamscan...\nThis might take a long time.")
    # Load proteome paths
    annotation_dir: Path = outdir.joinpath("kofam-annotation")
    # kofamdb_dir: Path = Path("/home/salvocos/Downloads/kofamscan/kofamdb")
    # kolist: Path = kofamdb_dir.joinpath("ko_list")
    kofam_outfmt: str = "detail-tsv"
    # prokaryote_profiles: Path = kofamdb_dir.joinpath("profiles/prokaryote.hal")
    proteome_paths: list[Path] = gda.systools.load_input_paths(proteomes_dir, fmt=ext)
    # KOfamScan ANNOTATION
    for p in proteome_paths:
        # print(p)
        bname = os.path.basename(p).replace(".faa", "")
        bname = f"{bname}.kofam.{kofam_outfmt}"
        # print(bname)
        gda.annotation.kofamscan_annotate(p, outdir=annotation_dir, kofam_path=kofam_path, profiles=kofam_profiles, kolist=ko_list, outprefix=bname, outfmt=kofam_outfmt, threads=threads, writelog=writelog)
        # break
    print("Kofam protein annotation completed successfully.")

    # Filter to only keep the significant hits
    kdecoder_input_dir: Path = annotation_dir.joinpath("kegg_decoder")
    gda.systools.makedir(kdecoder_input_dir)
    annotated_file_paths: list[Path] = gda.systools.load_input_paths(annotation_dir, fmt="tsv")
    prefix: str = ""
    kdecoder_input: Path = kdecoder_input_dir.joinpath("kdecoder.input.tsv")
    # recreate the file if needed
    if kdecoder_input.is_file():
        kdecoder_input.unlink(missing_ok=False)
        kdecoder_input.touch()

    # Extract K-object ids from kofam annotation files
    for p in annotated_file_paths:
        # base names have the following format:
        # <bacteria>.<id>.kofam.detail-tsv.tsv
        prefix = os.path.basename(p).split(".kofam.", maxsplit=1)[0]
        # Extract kegg ids
        gda.annotation.kofamscan2keggdecoder(p, outpath = kdecoder_input, kofam_fmt = "detail-tsv", prefix = prefix, best_hit_only = True, append = True)
        # Run KEGG-decoder
        # NOTE: KEGG-decoder requires its own environment be loaded in advance
        # The example command would be the following:
        # KEGG-decoder -i <kdecoder_input> -o kdecoder.output.list -v static

    # Retrieve KEGG files and convert them to JSON
    klist_file: Path = kdecoder_input_dir.joinpath("kegg_ids.txt")
    outdir_kegg: Path = kdecoder_input_dir.joinpath("raw")
    gda.systools.makedir(outdir_kegg)
    outdir_json: Path = kdecoder_input_dir.joinpath("json")
    gda.systools.makedir(outdir_json)

    # Create file containing only the ko_ids
    df: pd.DataFrame = pd.read_csv(kdecoder_input, sep="\t", header=None)
    df.sort_values(by=[1], ascending=[True], inplace=True)
    df.to_csv(klist_file, columns=[1], index=False, header=False)
    del df

    # Extract BRITE info from KofamScan output files
    brite_cnt_dir: Path = annotation_dir.joinpath("brite_counts")

    master: Path = extract_brite_hierarchies(annotated_file_paths, kegg_raw_dir=outdir_kegg, kegg_json_dir=outdir_json, outdir=brite_cnt_dir, kofam_fmt="detail-tsv", best_hit_only=True, exclude_uncomplete_brite=False, overwrite_kegg_obj=False, waittime=float64(1.50), concatenate=True)
    print(master)
    sys.exit("DEBUG: main")




if __name__ == "__main__":
    main()
