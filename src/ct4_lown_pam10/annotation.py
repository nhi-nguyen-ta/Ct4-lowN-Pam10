"""Utility methods for sequence annotation."""

import logging
import os
import subprocess
import sys
from collections import Counter
from json import dump, load
from pathlib import Path
from time import sleep
from typing import Any, TextIO
from shutil import which

# import bioservices
from bioservices import KEGG
from numpy import float64, uint32

from genome_data_analysis import systools

# Logger that will be used in this module
# It is child of the root logger and
# should be initialiazied using the function set_logger()
logger: logging.Logger = logging.getLogger()



def get_kegg_db_mapping(kresolver: KEGG, org: str = "") -> dict[str, str]:
    """Return a dictionary that maps KEGG DB names to their corresponding abbrev names used in the KEGG API"""

    dbgstr: str = f"""get_kegg_db_mapping :: START
    KEGG resolver:\t{kresolver}
    KEGG organisms ID:\t{org}
    """
    logger.debug(dbgstr)

    # Mapping of DB names to abbrev
    # The list of DB names and abbrev can be found at
    # https://www.kegg.jp/kegg/rest/keggapi.html
    KEGG_DB_NAMES2ABBREV: dict[str, str] = {
        "pathway":"path",
        "brite":"br",
        "module":"md",
        "orthology":"ko",
        "genome":"gn",
        "compound":"cpd",
        "glycan":"gl",
        "reaction":"rn",
        "rclass":"rc",
        "enzyme":"ec",
        "network":"ne",
        "variant":"hsa_var",
        "disease":"ds",
        "drug":"dr",
        "dgroup":"dg",
        "genes_vg":"vg",
        "genes_vp":"vp",
        "genes_ag":"ag"
        }

    # Retrieve KEGG organisms
    if len(org) > 0:
        if org in kresolver.organismIds:
            logging.info(f"genes_{org} added to the possible mappings with {org} as abbreviation.")
            # This mapping is only used for genes
            KEGG_DB_NAMES2ABBREV[f"genes_{org}"] = org
        else:
            logging.warning(f"The organism in input ({org}) is not valid. Use a valid KEGG organism code.")

    return KEGG_DB_NAMES2ABBREV




def extract_brite_info(kofam_raw: Path, kegg_raw_dir: Path, kegg_json_dir: Path, outdir: Path, kegg_db_abbrev: str, kofam_fmt: str = "detail-tsv", best_hit_only: bool = True, exclude_uncomplete_brite: bool = True, overwrite_kegg_obj: bool = False, waittime: float64 = float64(0.5)) -> Path:
    """Given an output file from KofamScan for each hit extra the BRITE ontologies, and write inot an output file."""

    dbgstr: str = f"""extract_brite_info :: START
    KofamScan output file: {kofam_raw}
    Directory with KEGG object: {kegg_raw_dir}
    Directory with KEGG object in JSON: {kegg_json_dir}
    Output directory: {outdir}
    KEGG DB abbreviation:\t{kegg_db_abbrev}
    Kofamscan output format:\t{kofam_fmt}
    Only use the first significant hit:\t{best_hit_only}
    Exclude some BRITE fields with poor support:\t{exclude_uncomplete_brite}
    Download again KEGG objects:\t{overwrite_kegg_obj}
    Time between requests to KEGG:\t{waittime}
    """
    logger.debug(dbgstr)

    # Check that the input file is valid
    if not kofam_raw.is_file():
        logger.error(f"The KofamScan result file is not valid:\n{kofam_raw}")
        sys.exit(-2)

    # Validate output formatting
    if kofam_fmt not in ["detail", "detail-tsv", "mapper", "mapper-oneline"]:
        logger.warning(f"Invalid output formating for kofamscan ({kofam_fmt}).")

    # Create output directory
    systools.makedir(outdir)

    logger.warning("For now only the \"detail-tsv\" is supported as kofamscan output format.")

    # Associate significant KO hits to contigs
    proteinseq2ko: dict[str, list[str]] = {}
    # Count occorrencies for each K number with a significant hit
    kocnt: Counter = Counter()
    proteinseq: str = ""
    koid: str = ""
    tmp_str_list: list[str] = []

    # Obtain KEGG resolver
    kresolver: KEGG = initialize_kegg_resolver()

    # This dictionary contains information about KEGG entry (e.g. K00500)
    # Specifically it contains:
    # SYMBOL: e.g. phhA, PAH (for K00500)
    # NAME: e.g. phenylalanine-4-hydroxylase [EC:1.14.16.1] (for K00500)
    koinfo: dict[str, tuple[str, str]] = dict()
    kegg_name: str = ""
    kegg_symbol: str = ""
    # List of BRITE lev 1 entries that could be excluded from processing
    brite_exclude_list: list[str] = ["Brite Hierarchies", "Not Included in Pathway or Brite"]

    # NOTE: this is hardcoded and is not a great solution
    bname: str = os.path.basename(kofam_raw)
    bname = bname.rsplit(".kofam", maxsplit=1)[0]
    # finalize the output file name
    bname = f"{bname}.brite.list.tsv"
    brite_list_path: Path = outdir.joinpath(bname)

    # Paths to json and raw KEGG object file
    jsonkpath: Path = Path()
    rawkpath: Path = Path()
    # Create the output with a entry for each BRITE entry at Lev 3
    # For example:
    # - Propanoate metabolism is lev3
    # - Carbohydrate metabolism is lev 2
    # - Metabolism is lev 3
    # Each entry could be as follows:
    # Brite-lev3	Brite-lev2	Brite-lev1	KEGG-Name	KEGG-Symbol	KEGG-Entry	Contig
    # Phenylalanine, tyrosine and tryptophan biosynthesis	Amino acid metabolism	Metabolism	phenylalanine-4-hydroxylase [EC:1.14.16.1]	phhA, PAH	K00500	Plasmid2_150
    # Note: Protein_sequence is the protein for which the given KEGG object was matched
    # Create the output file
    ofd: TextIO = open(brite_list_path, "wt")
    ofd.write("Brite-lev3\tBrite-lev2\tBrite-lev1\tKEGG-Name\tKEGG-Symbol\tKEGG-Entry\tProtein_sequence\n")
    # open the input file
    ifd: TextIO = open(kofam_raw, "rt")

    if kofam_fmt == "detail-tsv":
        # Lines in the kofam scan output would look as follows
        # #	gene name	KO	thrshld	score	E-value	"KO definition"
        # #	---------	------	-------	------	---------	-------------
        # *	Chromosome1_33	K09470	557.37	634.9	1.4e-191	"gamma-glutamylputrescine synthase [EC:6.3.1.11]"
        # *	Chromosome1_33	K01915	39.20	187.3	1.5e-55	"glutamine synthetase [EC:6.3.1.2]"
        # Chromosome1_33	K01949	828.73	125.9	3.7e-37	"glutamate---methylamine ligase [EC:6.3.4.12]"
        # Where lines starting with '*' represent significant hits

        # skip the header lines
        ifd.readline()
        ifd.readline()

        # process the hits
        for ln in ifd:
            if ln[0] == "*":
                tmp_str_list = ln.rstrip("\n").split("\t", maxsplit=6)
                proteinseq = tmp_str_list[1]
                koid = tmp_str_list[2]
                # print(f"Processing:\t{koid}")
                # Skip the hit if it is not the most significant
                if proteinseq in proteinseq2ko:
                    if best_hit_only:
                        logger.debug(f"Contig {proteinseq} already have an associated K object, hence the current one will be ignored.")
                        continue
                else:
                    proteinseq2ko[proteinseq] = [koid]
                # Increment the counter for the current K object
                kocnt[koid] += 1
                # Paths to the JSON and RAW KEGG object files
                rawkpath = kegg_raw_dir.joinpath(f"{koid}.txt")
                jsonkpath = kegg_json_dir.joinpath(f"{koid}.json")

                # Retrieve records if one of the two files is missing
                if (not os.path.isfile(rawkpath)) or (not os.path.isfile(jsonkpath)) or overwrite_kegg_obj:
                    # Retrieve the file again
                    get_keggobj(kresolver=kresolver, kid=koid, kegg_db_abbrev=kegg_db_abbrev, outdir_txt=kegg_raw_dir, outdir_json=kegg_json_dir, store_json=True, ignore_txt=False, overwrite=overwrite_kegg_obj)
                    # get_keggobj(kresolver=kresolver, kid=koid, outdir_kegg=kegg_raw_dir, outdir_json=kegg_json_dir, store_json=True, ignore_txt=False, overwrite=overwrite_kegg_obj)
                    sleep(waittime)

                # Extract the BRITE info and update the counters accordingly
                kegg_json_info = load(open(jsonkpath, "rt"))
                kegg_symbol = kegg_json_info["SYMBOL"]
                kegg_name = kegg_json_info["NAME"][0]
                # add the KEGG object info in dict
                if koid not in koinfo:
                    koinfo[koid] = (kegg_symbol, kegg_name)

                # Extract BRITE information
                # Level 1
                for kl1, vl1 in kegg_json_info["BRITE"]["root"].items():
                    # print("\nLev 1:")
                    # print(kl1, vl1)
                    if exclude_uncomplete_brite:
                        if kl1 in brite_exclude_list:
                            logger.debug(f"Excluding \"{kl1}\" section.")
                            continue
                    # Level 2
                    for kl2, vl2 in vl1.items():
                        # print("\nLev 2:")
                        # print(kl2, vl2)
                        # Level 3
                        for kl3, vl3 in vl2.items():
                            # print("\nLev 3:")
                            # print(kl3, vl3)
                            ofd.write(f"{kl3}\t{kl2}\t{kl1}\t{kegg_name}\t{kegg_symbol}\t{koid}\t{proteinseq}\n")
    else:
        logger.error("For now only the \"detail-tsv\" is supported as kofamscan output format.")
        sys.exit(-100)

    # Close the file with the Brite entry list
    ofd.close()
    ifd.close()

    # Compute total hits
    total_hits: uint32 = uint32(kocnt.total())
    logger.debug(f"Total KO hits for {os.path.basename(kofam_raw)}:\t{total_hits}")
    logger.debug(f"Uniq K objects:\t{len(kocnt.keys())}")
    logger.debug(f"Total hits:\t{total_hits}")
    logger.debug(f"KO with most hits:\t{kocnt.most_common(n=5)}")

    # Create the output with counts for KEGG objects
    bname = os.path.basename(kofam_raw)
    bname = bname.rsplit(".", maxsplit=1)[0]
    # finalize the output file name
    bname = f"{bname}.kcnt.tsv"
    # create the output file
    outpath_kcnt: str = os.path.join(outdir, bname)
    # The file KO counts will have the following format
    # KEGG-Entry	Hits	Hits-percentage	KEGG-Symbol	KEGG-Name
    # K00500 50 1.15	phhA, PAH	phenylalanine-4-hydroxylase [EC:1.14.16.1]
    # The percentage is computed on the total number of hits in the KofamScan output file
    ofd = open(outpath_kcnt, "tw")
    ofd.write("KEGG-Entry\tHits\tHits-percentage\tKEGG-Symbol\tKEGG-Name\n")
    # Write KEGG object counts in the output file
    for kid, cnt in sorted(kocnt.items()):
        # ofd.write("KEGG-Entry\tHits\tHits-percentage\tKEGG-Symbol\tKEGG-Name\n")
        ofd.write(f"{kid}\t{cnt:d}\t{(cnt/total_hits)*100.}\t{koinfo[kid][0]}\t{koinfo[kid][1]}\n")
    ofd.close()

    # Return the path to the file with BRITE information
    return brite_list_path



def download_kegg_objects(klist_file: str, outdir_txt: str, outdir_json: str, store_json: bool = False, ignore_txt: bool = False, overwrite: bool = False, waittime: float64 = float64(0.50)) -> tuple[uint32, uint32, uint32, uint32, uint32]:
    """Download KEGG entries in raw TXT and/or JSON."""

    dbgstr: str = f"""download_kegg_objects :: START
    List with K IDs: {klist_file}
    Output directory: {outdir_txt}
    Output directory for JSON files: {outdir_json}
    Store a JSON file:\t{store_json}
    Do not create a TXT file:\t{ignore_txt}
    Overwrite existing file:\t{overwrite}
    Time between requests to KEGG:\t{waittime}
    """
    logger.debug(dbgstr)

    # Load the KEGG ids from the input file
    # Each KEGG object identifier starts with a 'K' and is followed by 5 digits (6 chars in total)
    # Following is an example: K00500
    klist: list[str] = []
    kid: str = ""
    with open(klist_file, "rt") as ifd:
        for ln in ifd:
            # The following conditions must apply for the ID to be valid:
            # String length must be 6
            # The first character must be a 'k' or 'K'
            # The chars from position 1 to 5 must be digits
            ln = ln.rstrip("\n")
            if (len(ln) != 6) or (ln[0] not in ["k", "K"]) or (not ln[1:6].isnumeric()):
                logger.error(f"The ID {repr(ln)} is not valid. Valid KEGG objects have the following format:\nK00000")
                sys.exit(-10)
            else:
                klist.append(ln)

    # sort the list
    klist.sort()
    logger.debug(f"Loaded KEGG IDs:\t{len(klist)}")
    # Initialize some counters
    downloaded_raw: uint32 = uint32(0)
    downloaded_json: uint32 = uint32(0)
    skipped_raw: uint32 = uint32(0)
    skipped_json: uint32 = uint32(0)
    total_ids: uint32 = uint32(len(klist))
    cnt: uint32 = uint32(0)

    # Initialize the KEGG resolver
    kservice: KEGG = initialize_kegg_resolver()

    # Tuple in the format returned from get_keggobj()
    restpl: tuple[str, str, bool, bool] = ("", "", False, False)
    for kid in klist:
        cnt = cnt + 1
        print(f"Processing ID {kid}\t{cnt}/{total_ids}")
        restpl = get_keggobj(kresolver=kservice, kid=kid, outdir_txt=outdir_txt, outdir_json=outdir_json, store_json=store_json, ignore_txt=ignore_txt, overwrite=overwrite)
        # update the counters for raw KEGG objets
        if restpl[2] is True:
            downloaded_raw = downloaded_raw + 1
        else:
            skipped_raw = skipped_raw + 1
        # and for JSON ones
        if restpl[2] is True:
            downloaded_json = downloaded_json + 1
        else:
            skipped_json = skipped_json + 1
        # Wait for between requests
        if restpl[2] or restpl[3]:
            sleep(waittime)

    logger.debug(total_ids, downloaded_raw, skipped_raw, downloaded_json, skipped_json)
    return (total_ids, downloaded_raw, skipped_raw, downloaded_json, skipped_json)



def get_keggobj(kresolver: KEGG, kid: str, kegg_db_abbrev: str = "ko", outdir_txt: Path = Path(), outdir_json: Path = Path(), store_json: bool = False, ignore_txt: bool = False, overwrite: bool = False) -> tuple[Path, Path, bool, bool]:
    """Obtain the object using a KEGG object ID"""
    # The information are retrieved using bioservices
    # https://bioservices.readthedocs.io/en/main/index.html
    # Alternatively the KEGG REST API could be directly used
    # https://www.kegg.jp/kegg/rest/keggapi.html

    # Add valid DB list

    dbgstr: str = f"""get_keggobj :: START
    KEGG Resolver:\t{kresolver}
    KEGG object ID:\t{kid}
    KEGG Database abbrevation:\t{kegg_db_abbrev}
    Output directory: {outdir_txt}
    Output directory for files JSON: {outdir_json}
    Store a JSON file:\t{store_json}
    Do not create a TXT file:\t{ignore_txt}
    Overwrite existing file:\t{overwrite}
    """
    logger.debug(dbgstr)

    # Check consistency among input parameters
    if store_json:
        if not outdir_json.is_dir():
            logger.error("You must specify the output directory for the JSON files.")
            sys.exit(-2)
        else:
            systools.makedir(outdir_json)
    if not ignore_txt:
        if not outdir_txt.is_dir():
            logger.error("You must specify the output directory for the RAW KEGG files.")
            sys.exit(-2)
        else:
            systools.makedir(outdir_txt)

    # Create the output paths
    outtxt: Path = outdir_txt.joinpath(f"{kid}.txt")
    outjson: Path = outdir_json.joinpath(f"{kid}.json")
    downloaded_raw: bool = False
    downloaded_json: bool = False

    # Avoid downloading the file if not required
    if not overwrite:
        if outtxt.is_file() and outjson.is_file():
            return (outtxt, outjson, False, False)

    # API get string
    kegg_get_str: str = f"{kegg_db_abbrev}:{kid}"
    # rawdata: Any = kresolver.get(f"ko:{kid}", parse=False)
    rawdata: Any = kresolver.get(kegg_get_str, parse=False)

    # NOTE: we cannot directy give the type "str" when we define rawdata, because the HTML "get" could return something that is not str

    """
    # Write the TXT file if required
    # Skip creation of the textfile
    if ignore_txt:
        outtxt = ""
    else:
        # Do not download unleass required
        if os.path.isfile(outtxt) and (not overwrite):
            logger.info(f"The raw file for  {kid} was previously downloaded, and will be skipped.")
        else:
            with open(outtxt, "wt") as ofdtxt:
                ofdtxt.write(rawdata)
                downloaded_raw = True
    """

    # Write the TXT file if required
    # Skip creation of the textfile
    if not ignore_txt:
        # Do not download unleass required
        if outtxt.is_file() and (not overwrite):
            logger.info(f"The raw file for  {kid} was previously downloaded, and will be skipped.")
        else:
            with open(outtxt, "wt") as ofdtxt:
                ofdtxt.write(rawdata)
                downloaded_raw = True
    else:
        # Set it to a directory
        outtxt = outdir_txt


    # Parse the file and store in JSON format if requested
    '''
    TODO: this includes the paring of the BRITE sections which is not necessarily required
    An extra parameter to parse BRITE might be included.
    '''

    if store_json:
        parsed_brite: dict[str, dict[str, dict[str, dict[str, str]]]] = dict()
        # Do not download unleass required
        if outjson.is_file() and (not overwrite):
            logger.info(f"The JSON file for {kegg_db_abbrev}:{kid} was previously downloaded, and will be skipped.")
        else:
            parsed_kegg_obj: dict[str, Any] = kresolver.parse(rawdata)

            '''
            This is not a good solution,
            a function parameter (e.g. parse_brite) could be used to decide if
            BRITE should be parsed or not
            '''
            if "BRITE" in parsed_kegg_obj.keys():
                if kegg_db_abbrev == "ko":
                    # The function from bioservice does not parse the BRITE section
                    # which is return as raw string as value for the "BRITE" key in the dictionay
                    # We will parse this text using an in-house parser and update the value for the BRITE key
                    parsed_brite = parse_brite_lines(brite_txt=parsed_kegg_obj["BRITE"], kid=kid)
                    # update the value in the dictionary before dumping it
                    parsed_kegg_obj["BRITE"] = parsed_brite
            with open(outjson, "w") as ofdjson:
                dump(parsed_kegg_obj, ofdjson, indent=4)
                downloaded_json = True
    else:
        # Set it to a directory
        outjson = Path()

    # return the paths to the generated files
    return (outtxt, outjson, downloaded_raw, downloaded_json)



def filter_kofamscan_hits(kofam_raw: str, outdir: str, kofam_fmt: str = "mapper", best_hit_only: bool = True) -> str:
    """Filter the output file from Kofamscan and keep only the hit with significant scores"""

    dbgstr: str = f"""filter_kofamscan_hits :: START
    Kofamscan output file: {kofam_raw}
    Output dir: {outdir}
    Kofamscan output format: {kofam_fmt}
    Only keep the best hit:\t{best_hit_only}
    """
    logger.debug(dbgstr)

    # Validate output formatting
    if kofam_fmt not in ["detail", "detail-tsv", "mapper", "mapper-oneline"]:
        logger.warning(f"Invalid output formating for kofamscan ({kofam_fmt}).")

    logger.warning("For now only the \"detail-tsv\" is supported as kofamscan output format.")

    # Associate significant KO hits to contigs
    ctg2ko: dict[str, list[str]] = {}
    ctg: str = ""
    koid: str = ""
    bname: str = os.path.basename(kofam_raw)
    tmp_str_list: list[str] = []
    outpath: str = ""

    # open the input and output files
    ifd: TextIO = open(kofam_raw, "rt")

    if kofam_fmt == "detail-tsv":
        # Lines in the kofam scan output would look as follows
        # #	gene name	KO	thrshld	score	E-value	"KO definition"
        # #	---------	------	-------	------	---------	-------------
        # *	Chromosome1_33	K09470	557.37	634.9	1.4e-191	"gamma-glutamylputrescine synthase [EC:6.3.1.11]"
        # *	Chromosome1_33	K01915	39.20	187.3	1.5e-55	"glutamine synthetase [EC:6.3.1.2]"
        # Chromosome1_33	K01949	828.73	125.9	3.7e-37	"glutamate---methylamine ligase [EC:6.3.4.12]"
        # Where lines starting with '*' represent significant hits

        # print(f"bname.before:\t{bn}")
        # remove the file extension
        bname = bname.rsplit(".", maxsplit=1)[0]
        # finalize the output file name
        bname = f"{bname}.filtered.tsv"
        # create the output file
        outpath = os.path.join(outdir, bname)
        ofd: TextIO = open(outpath, "tw")

        # skip the header lines
        ifd.readline()
        ifd.readline()

        # process the hits
        for ln in ifd:
            if ln[0] == "*":
                # print(ln)
                if not best_hit_only:
                    ofd.write(ln)
                else:
                    # only allow the first signicant hit to be written
                    tmp_str_list = ln.split("\t", maxsplit=3)[1:3]
                    ctg = tmp_str_list[0]
                    koid = tmp_str_list[1]
                    # Insert the KO element only if never seen before for the current contig
                    if ctg not in ctg2ko:
                        ctg2ko[ctg] = [koid]
                        ofd.write(ln)
                    else:
                        logger.warning(f"The KO number {koid} was previously found for contig {ctg} and will be ignored.")
    else:
        logger.error("For now only the \"detail-tsv\" is supported as kofamscan output format.")
        sys.exit(-100)

    ifd.close()
    ofd.close()

    return outpath



def kofamscan2keggdecoder(kofam_raw: Path, outpath: Path, kofam_fmt: str = "mapper", prefix: str = "", best_hit_only: bool = True, append: bool = False) -> Path:
    """Conver Kofamscan output to Keggdecoder input file"""

    dbgstr: str = f"""kofamscan2keggdecoder :: START
    Kofamscan output file: {kofam_raw}
    Output file: {outpath}
    Kofamscan output format: {kofam_fmt}
    Prefix: {prefix}
    Only keep the best hit:\t{best_hit_only}
    Write the file in append mode:\t{append}
    """
    logger.debug(dbgstr)

    # Validate output formatting
    if kofam_fmt not in ["detail", "detail-tsv", "mapper", "mapper-oneline"]:
        logger.warning(f"Invalid output formating for kofamscan ({kofam_fmt}).")

    logger.warning("For now only the \"detail-tsv\" is supported as kofamscan output format.")

    # Associate significant KO hits to contigs
    ctg2ko: dict[str, list[str]] = {}
    ctg: str = ""
    koid: str = ""
    tmp_str_list: list[str] = []

    # open the input and output files
    ifd: TextIO = open(kofam_raw, "rt")

    if kofam_fmt == "detail-tsv":
        # Lines in the kofam scan output would look as follows
        # #	gene name	KO	thrshld	score	E-value	"KO definition"
        # #	---------	------	-------	------	---------	-------------
        # *	Chromosome1_33	K09470	557.37	634.9	1.4e-191	"gamma-glutamylputrescine synthase [EC:6.3.1.11]"
        # *	Chromosome1_33	K01915	39.20	187.3	1.5e-55	"glutamine synthetase [EC:6.3.1.2]"
        # Chromosome1_33	K01949	828.73	125.9	3.7e-37	"glutamate---methylamine ligase [EC:6.3.4.12]"
        # Where lines starting with '*' represent significant hits

        # Open the file in append mode
        ofd: TextIO = open(outpath, "ta")
        if not append:
            ofd.close()
            open(outpath, "tw")

        # skip the header lines
        ifd.readline()
        ifd.readline()

        # process the hits
        for ln in ifd:
            if ln[0] == "*":
                # extract the required info
                tmp_str_list = ln.split("\t", maxsplit=3)[1:3]
                ctg = tmp_str_list[0]
                koid = tmp_str_list[1]
                # Insert the KO element only if never seen before for the current contig
                if ctg not in ctg2ko:
                    ctg2ko[ctg] = [koid]
                else:
                    if best_hit_only:
                        logger.warning(f"The KO number {koid} was previously found for contig {ctg} and will be ignored.")
                        continue
                    else:
                        ctg2ko[ctg].append(koid)
                # write the info to the output file
                if len(prefix) > 0:
                    ofd.write(f"{prefix}_{ctg}\t{koid}\n")
                else:
                    ofd.write(f"{ctg}\t{koid}\n")
    else:
        logger.error("For now only the \"detail-tsv\" is supported as kofamscan output format.")
        sys.exit(-100)

    ifd.close()
    ofd.close()

    return outpath



def kofamscan_annotate(inseq: Path, outdir: Path, kofam_path: Path, profiles: Path, kolist: Path, outprefix: str = "kofamscan", outfmt: str = "detail", threads: uint32 = uint32(8), writelog: bool = True) -> Path:
    """ Annotate input proteins with KEGG information using KofamScan.
        ftp://ftp.genome.jp/pub/tools/kofam_scan/
    """
    # Set the log file path
    logpath: str = os.path.join(outdir, f"log.{outprefix}.txt")
    currentlev: int = logger.level

    dbgstr: str = f"""kofamscan_annotate :: START
    Protein file: {inseq}
    Output directory: {outdir}
    Path to the Kofamscan executables: {kofam_path}
    Profiles: {profiles}
    KO list: {kolist}
    Output name:\t{outprefix}
    Output type:\t{outfmt}
    Write log file:\t{writelog}
    Threads:\t{threads}
    """
    logger.debug(dbgstr)

    # create the output directory if required
    systools.makedir(outdir)
    # Directory for kofamscam temporary files
    tmpdir: Path = outdir.joinpath(f"tmp.{outprefix}")

    flogger: logging.Logger = logging.Logger("")
    if writelog:
        flogger = systools.create_flogger(logpath, loggername = f"{__name__}.file_logger", lev = currentlev, mode="a", propagate=False)
        flogger.log(20, dbgstr)

    if not inseq.is_file():
        logger.error(f"The file with the input protein sequences was not found:\n{inseq}")
        sys.exit(-2)

    if not kofam_path.is_file():
        logger.error(f"The path to the kofamscan executables (exec_annotate) is not valid:\n{kofam_path}")
        sys.exit(-2)

    # Validate output formatting
    if outfmt not in ["detail", "detail-tsv", "mapper", "mapper-oneline"]:
        logger.warning(f"Invalid output formating for kofamscan ({outfmt}).\nThe default formating will be used: 'detail'")
        outfmt = "detail"

    # check that profiles and KO list file are valid
    if not profiles.exists():
        logger.error(f"The file or dir with HMM profiles is not valid:\n{profiles}")
        sys.exit(-2)
    if not kolist.exists():
        logger.error(f"The file with the KO list is not valid:\n{kolist}")
        sys.exit(-2)

    # set output paths
    outpath: Path = outdir.joinpath(f"{outprefix}.txt")
    # set to tsv format if not 'detail'
    if outfmt.endswith("tsv"):
        outpath = outdir.joinpath(f"{outprefix}.tsv")

    # kofamscan example
    # exec_annotation <input.faa> -o <outfile.tsv> --format mapper --cpu <threads> -k <path/to/ko_list> -p <path/to/profiles/prokaryote.hal>
    cmd: str = f"{kofam_path} {inseq} -o {outpath} --format {outfmt} --cpu {threads} -k {kolist} -p {profiles} --tmp-dir {tmpdir}"

    logger.debug(f"\nKofamScan CMD:\n{cmd}")
    if writelog:
        # Write also in the log file
        flogger.log(currentlev, f"\nKofamScan CMD:\n{cmd}")

    #execute the system call
    completed: subprocess.CompletedProcess = subprocess.run(cmd, shell=True, capture_output=True)

    if writelog:
        # Write also in the log file
        flogger.log(currentlev, f"\nEXIT CODE:\t{completed.returncode}\n\
        STDOUT:\n{completed.stdout.decode()}\n\
        STDERR:\n{completed.stderr.decode()}\
        ")

    if completed.returncode != 0:
        logger.error(f"KofamScan failed with exit code {completed.returncode}\n{completed.stderr}")
        sys.exit(completed.returncode)

    # Return the path to the file with proteins
    return outpath



def initialize_kegg_resolver() -> KEGG:
    """Initialize the KEGG resolver"""
    logger.debug("initialize_kegg_resolver :: START")

    # Initialize the resolver for KEGG
    return KEGG(verbose=False, cache=True)



def map_kegg_db2abbrev(kdb: str, org: str, kresolver: KEGG) -> tuple[str, list[str]]:
    """
    Map a KEGG database name to its abbrevation
    A list of KEGG DBs can found at
    https://www.kegg.jp/kegg/rest/keggapi.html
    BioService might use different db names
    https://bioservices.readthedocs.io/en/main/index.html

        Parameters:
                db (str): KEGG DB name
                org (str): KEGG organism ID
                kresolver (str): KEGG connector
        Returns:
                mapping_info tuple[str, list[str]]: abbrev of DB name, and other information
                    mapping_info[0] : str
                        KEGG API abbrev for DB name (e.g. path for pathway)
                    mapping_info[1] : list[str]
                        Examples of kid prefix to be used in API calls (e.g. map03010 for pathaways)
    """

    dbgstr: str = f"""map_kegg_db2abbrev :: START
    Database name:\t{kdb}
    Organism:\t{org}
    KEGG connection:\t{kresolver}
    """
    logger.debug(dbgstr)
    # Mapping dictionary containing the kegg abbrev
    kegg_db_name2abbrev: dict[str, str] = get_kegg_db_mapping(kresolver, org)

    # Used to identify the special case
    # in which the DB is a Organisms ID
    org_as_db: bool = False

    # Make sure that the name of the KEGG DB is valid
    if kdb == "genes":
        gene_abbrevs: list[str] = ["genes_vg", "genes_vp", "genes_ag", "gene_<org>"]
        logging.error(f"\nPlease us one of the following DB names instead:\n{gene_abbrevs}")
        sys.exit(-10)
    elif(kdb not in kegg_db_name2abbrev.keys()):
        print(f"genes_{org}" not in kegg_db_name2abbrev.keys())
        # If the special case in which a organisms ID was used do nothing
        if f"genes_{org}" not in kegg_db_name2abbrev.keys():
            logging.error(f"\n{kdb} is not a valid KEGG DB name, use one of the following instead:\n{list(kegg_db_name2abbrev.keys())}")
            sys.exit(-10)
        else:
            org_as_db = True

    # print(f"Abbrev:\t{abbrev}")
    # Map each KEGG abbrev to examples of entries for the requested DB
    # This is corrent only if we use an organisms as the DB name
    # hence org contains a valid KEGG organism ID
    abbrev: str = org
    if not org_as_db:
        abbrev = kegg_db_name2abbrev[kdb]
    # examples of kid to be used
    kid_prefixes: list[str] = []
    if abbrev == "path":
        # KEGG pathway maps
        # path:map03010
        # path:hsa04930 (not in bioservices)
        kid_prefixes.append("map03010")
        kid_prefixes.append("hsa04930")
    elif abbrev == "br":
        # BRITE functional hierarchies
        # brite br (br, jp, ko, <org>)
        # API call examples
        # br:br08303
        # br:jp08303 (only for Japanese)
        # br:ko00003 (not in bioservices)
        kid_prefixes.append("br08303")
        kid_prefixes.append("jp08303")
        kid_prefixes.append("ko00003")
    elif abbrev == "md":
        logging.warning("...and organisms ID could be provided (e.g. hsa)")
        # KEGG modules
        # module md (M number) | <org>_M
        # API call examples
        # md:M00020
        # the presence of the module for a given species
        # can be verified from the module file
        # Which lists all the species for which the module is available
        # under the COMPLETE section
        # md:hsa_M00020
        kid_prefixes.append("M00020")
        kid_prefixes.append("hsa_M00020")
    elif abbrev == "ko":
        # KEGG orthology
        # orthology ko (K number)
        # NOTE: koIds in Bioservices returns a list with the valid K-numbers
        # API call examples
        # ko:K001234
        kid_prefixes.append("K001234")
    elif abbrev == "cpd":
        # Small molecules (compounds)
        # compound cpd (C number)
        # API call examples
        # cpd:C00031
        kid_prefixes.append("C00031")
    elif abbrev == "gl":
        # Glycans
        # glycan gl (G number)
        # API call examples
        # gl:G00109
        kid_prefixes.append("G00109")
    elif abbrev == "rn":
        # Biochemical reactions
        # reaction rn (R number)
        # API call examples
        # rn:R00259
        kid_prefixes.append("R00259")
    elif abbrev == "rc":
        # Reaction class
        # rclass rc  RC
        # API call examples
        # rc:RC00046
        kid_prefixes.append("RC00046")
    elif abbrev == "ec":
        # Enzyme nomenclature
        # Enzyme ec
        # API call examples
        # ec:2.7.10.1
        kid_prefixes.append("2.7.10.1")
    elif abbrev == "gn":
        # KEGG Organisms (Genomes)
        # Mostly information about the organism
        # genome size, taxonomy etc.
        # genome gn (T number)
        # API call examples
        # gn:T01001
        kid_prefixes.append("T01001")
    elif abbrev == "ne":
        # Network elements
        # network ne (N number | nt00001)
        # NOTE: this could require a separate methods
        # API call examples
        # ne:N00002 (network element)
        # ne:nt06210 (network variation map)
        kid_prefixes.append("N00002")
        kid_prefixes.append("nt06210")
    elif abbrev == "hsa_var":
        # Human gene variants
        # hsa_var variant id
        # API call examples
        # hsa_var:118v1
        kid_prefixes.append("118v1")
    elif abbrev == "ds":
        # Diseases (Human diseases)
        # disease ds (H number)
        # API call examples
        # ds:H00004
        kid_prefixes.append("H00004")
    elif abbrev == "dr":
        # Drug
        # drug dr (D number)
        # API call examples
        # dr:D01441
        kid_prefixes.append("D01441")
    elif abbrev == "dg":
        # Drug groups
        # dgroup dg (DG number)
        # API call examples
        # dg:DG00710
        kid_prefixes.append("DG00710")
    elif abbrev == "vg":
        # viral gene or proteins
        # API call examples
        # vg:155971
        kid_prefixes.append("155971")
    elif abbrev == "vp":
        # mature peptides in viruses
        # API call examples
        # vp:155971-1
        kid_prefixes.append("155971-1")
    elif abbrev == "ag":
        # Functionally characterized proteins
        # API call examples
        # ag:CAA76703
        kid_prefixes.append("CAA76703")
    elif abbrev == org:
        # Gene for a valid organism
        # <org>:<gene_id>
        # API call examples
        # hsa:3643
        kid_prefixes.append("3643")

    # Return abbrevation and examples of prefixes
    return (abbrev, kid_prefixes)



def parse_brite_lines(brite_txt: str, kid: str) -> dict[str, dict[str, dict[str, dict[str, str]]]]:
    """Parse lines containing BRITE levels.
    brite_str: the content of the BRITE section of a KEGG object.

    Following is an example of the BRITE content for the KEGG object K00001:

    BRITE       KEGG Orthology (KO) [BR:ko00001]
                 09100 Metabolism
                  09101 Carbohydrate metabolism
                   00010 Glycolysis / Gluconeogenesis
                    K00001  E1.1.1.1, adh; alcohol dehydrogenase
                   00620 Pyruvate metabolism
                    K00001  E1.1.1.1, adh; alcohol dehydrogenase
                  09103 Lipid metabolism
                   00071 Fatty acid degradation
                    K00001  E1.1.1.1, adh; alcohol dehydrogenase
                  09105 Amino acid metabolism
                   00350 Tyrosine metabolism
                    K00001  E1.1.1.1, adh; alcohol dehydrogenase
                  09108 Metabolism of cofactors and vitamins
                   00830 Retinol metabolism
                    K00001  E1.1.1.1, adh; alcohol dehydrogenase
                  09111 Xenobiotics biodegradation and metabolism
                   00625 Chloroalkane and chloroalkene degradation
                    K00001  E1.1.1.1, adh; alcohol dehydrogenase
                   00626 Naphthalene degradation
                    K00001  E1.1.1.1, adh; alcohol dehydrogenase
                   00980 Metabolism of xenobiotics by cytochrome P450
                    K00001  E1.1.1.1, adh; alcohol dehydrogenase
                   00982 Drug metabolism - cytochrome P450
                    K00001  E1.1.1.1, adh; alcohol dehydrogenase
                Enzymes [BR:ko01000]
                 1. Oxidoreductases
                  1.1  Acting on the CH-OH group of donors
                   1.1.1  With NAD+ or NADP+ as acceptor
                    1.1.1.1  alcohol dehydrogenase
                     K00001  E1.1.1.1, adh; alcohol dehydrogenase
    """

    dbgstr: str = f"""parse_brite_lines :: START
    KEGG ID:\t{kid}\n
    BRITE content: \n{brite_txt}\n
    """
    logger.debug(dbgstr)

    # set the string of white spaces used
    # to identify the BRITE level
    lev1start: str = "             "
    lev2start: str = "              "
    lev3start: str = "               "
    lev4start: str = "                "
    # This is only present in extra fields as described below
    lev5start: str = "                 "

    # Aside from the main brite levels also introduce other info
    # such us Enzyme, Ion channels and so on
    # Such levels start with a 12-chars long string
    # An example of an extra section is the following:
    # Enzymes [BR:ko01000]
    # The internal entries have the length levels 1 to 4 above, so such variables can be reused for parsing
    extra_entity_start: str = "            "
    tmpstr: str = ""
    # The dictionary will contain all the information
    # From level 0 (e.g. the line contaning BRITE keyword)
    # To level 3
    # Level 4 is ignored because it is repeated and contains info already present in the Kegg object
    brite_dict: dict[str, dict[str, dict[str, dict[str, str]]]] = dict()
    # Variables used to keep track of the BRITE level during parsing
    current_lev0: str = "root"
    current_lev1: str = ""
    current_lev2: str = ""
    current_lev3: str = ""
    tmpsplit: list[str] = []
    # Check if we are in a extra section
    extra_section: bool = False

    # Start parsing the BRITE info
    for ln in brite_txt.split("\n"):
        if ln.startswith(lev5start):
            # This section of the BRITE hierarchy is ignored i nthis version of the parser
            tmpstr = ln.lstrip(" ")
            if not extra_section:
                logger.error(f"This string is expected to the value of a 'extra' BRITE section.\n{tmpstr}")
                sys.exit(-10)
            # pass
        elif ln.startswith(lev4start):
            # level 4 example: K00500  phhA, PAH; phenylalanine-4-hydroxylase
            tmpstr = ln.lstrip(" ")
            logger.debug(f"LEV 4:\n\t{tmpstr}")
            if tmpstr.startswith(kid):
                # print(f"{current_lev0}\n\t{current_lev1}\n\t\t{current_lev2}\n\t\t{current_lev3}")
                # The for Lev 2 could be overwritten hence an extra check on the dict is required
                if current_lev0 in brite_dict:
                    # Avoid lev 1 to be overwritten
                    if current_lev1 in brite_dict[current_lev0]:
                        logger.warning(f"The value for Lev 1 ({current_lev1}) could be overwritten")
                        # check lev 2
                        if current_lev2 in brite_dict[current_lev0][current_lev1]:
                            logger.warning(f"The value for Lev 2 ({current_lev2}) could be overwritten")
                            # check lev 3
                            if current_lev3 in brite_dict[current_lev0][current_lev1][current_lev2]:
                                logger.warning(f"The value for Lev 3 ({current_lev3}) could be overwritten")
                            else:
                                brite_dict[current_lev0][current_lev1][current_lev2][current_lev3] = tmpstr
                        else:
                            brite_dict[current_lev0][current_lev1][current_lev2] = {current_lev3:tmpstr}
                    else:
                        brite_dict[current_lev0][current_lev1] = {current_lev2:{current_lev3:tmpstr}}
                # if current_lev1 in brite_dict[current_lev0]:
                else:
                    brite_dict[current_lev0] = {current_lev1:{current_lev2:{current_lev3:tmpstr}}}
            else:
                if extra_section:
                    logger.warning("This is an extra field in the BRITE hierarchy and is ignored in this version of the parser.")
                    continue
                else:
                    logger.error(f"LEV 4 is expected to start with a KEGG ID {kid}")
                    print(tmpstr)
                    sys.exit(-10)
        elif ln.startswith(lev3start):
            # lev 3 example: 00010 Glycolysis / Gluconeogenesis
            tmpstr = ln.lstrip(" ")
            logger.debug(f"LEV 3:\n\t{tmpstr}")
            tmpsplit = tmpstr.split(" ", maxsplit=1)
            if len(tmpsplit) == 1:
                logger.error(f"The field at level 3 does not contain spaces\n{tmpstr}\nthe parsing will be stopped at this level.")
            else:
                current_lev3 = tmpsplit[1]
                # Add to the dictionary without the numeric node
                # brite_dict[current_lev0][current_lev1][current_lev2] = tmpstr
                # brite_dict[current_lev0][current_lev1][current_lev2] = current_lev3
                # brite_dict[current_lev0] = {current_lev1:{current_lev2:current_lev3}}
            # break
        elif ln.startswith(lev2start):
            # lev 2 example: 09101 Carbohydrate metabolism
            tmpstr = ln.lstrip(" ")
            logger.debug(f"LEV 2:\n\t{tmpstr}")
            tmpsplit = tmpstr.split(" ", maxsplit=1)
            if len(tmpsplit) == 1:
                if not extra_section:
                    logger.error(f"The field at level 2 does not contain spaces\n{tmpstr}\nthe parsing will be stopped at this level.")
            else:
                current_lev2 = tmpsplit[1]
            # break
        elif ln.startswith(lev1start):
            # lev 1 example: 09100 Metabolism
            tmpstr = ln.lstrip(" ")
            logger.debug(f"LEV 1:\n\t{tmpstr}")
            tmpsplit = tmpstr.split(" ", maxsplit=1)
            if len(tmpsplit) == 1:
                logger.error(f"The field at level 1 does not contain spaces\n{tmpstr}\nthe parsing will be stopped at this level.")
            else:
                current_lev1 = tmpsplit[1]
            # current_lev1 = tmpsplit[1]
            # break
        elif ln.startswith(extra_entity_start):
            # Example fo extra level start: Exosome [BR:ko04147]
            tmpstr = ln.lstrip(" ")
            if tmpstr[-1] == "]":
                logger.debug(f"LEV 0: EXTRA SECTION\n\t{tmpstr}")
                extra_section = True
                current_lev0 = tmpstr
            else:
                logger.error("Extra sections start levels must and with a ']' symbol.")
                sys.exit("DEBUG")
            # break
        # break

    return brite_dict



def prodigal(inseq: str, outdir: str = os.getcwd(), outname: str = "prodigal_prediction", mode: str = "single", quiet: bool = False, writelog: bool = False) -> str:
    """Predict proteins from input DNA sequences using Prodigal.
    mode: [single|meta] says if the DNA sequences come from a metagenomic or single genome study
    """
    # Set the log file path
    logpath: str = os.path.join(outdir, f"log.{outname}.txt")
    currentlev: int = logger.level

    dbgstr: str = f"""prodigal :: START
    Protein file: {inseq}
    Output directory: {outdir}
    Output name:\t{outname}
    Prodigal mode:\t{mode}
    Quite mode:\t{quiet}
    Write log file:\t{writelog}
    """
    logger.debug(dbgstr)

    # create the output directory if required
    systools.makedir(outdir)

    flogger: logging.Logger = logging.Logger("")
    if writelog:
        flogger = systools.create_flogger(logpath, loggername = f"{__name__}.file_logger", lev = currentlev, mode="a", propagate=False)
        flogger.log(20, dbgstr)

    if not os.path.isfile(inseq):
        logger.error(f"The file with the input protein sequences was not found:\n{inseq}")
        sys.exit(-2)

    proteins_outpath: str = os.path.join(outdir, f"{outname}.faa")
    gbk_outpath: str = os.path.join(outdir, f"{outname}.gbk")
    options: list[str] = ["-p", mode]
    if quiet:
        options.append("-q")

    # EXAMPLES of command to be executed
    # prodigal -i <input.fna> -a <predicted_proteins.faa> -p single -o <main_output.gbk> -q
    cmd: str = f"prodigal -i {inseq} -a {proteins_outpath} -o {gbk_outpath} {' '.join(options)}"
    logger.debug(f"\nProdigal CMD:\n{cmd}")
    if writelog:
        # Write also in the log file
        flogger.log(currentlev, f"\nProdigal CMD:\n{cmd}")

    #execute the system call
    completed: subprocess.CompletedProcess = subprocess.run(cmd, shell=True, capture_output=True)

    if writelog:
        # Write also in the log file
        flogger.log(currentlev, f"\nCMD:\t{completed.args}\n\
        EXIT CODE:\t{completed.returncode}\n\
        STDOUT:\n{completed.stdout.decode()}\n\
        STDERR:\n{completed.stderr.decode()}\
        ")

    if completed.returncode != 0:
        logger.error(f"Prodigal failed with exit code {completed.returncode}\n{completed.stderr}")
        sys.exit(completed.returncode)
    else:
        print("Protein prediction completed successfully.")

    # Return the path to the file with proteins
    return proteins_outpath
