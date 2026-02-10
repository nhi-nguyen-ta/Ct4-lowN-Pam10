"""Script to analyze Hiruma's bacterial genomes"""
from polars.testing.parametric import columns
import logging
import os
import sys
from os.path import sep
from pathlib import Path

import matplotlib.pyplot as plt
import polars as pl
import seaborn as sns
from numpy import float64, uint32
from rich.logging import RichHandler
from typing_extensions import TextIO
import argparse

import ct4_lown_pam10 as ct

# Set the logger
logger = logging.getLogger(__name__)



### FUNCTIONS ###
def get_params() -> argparse.Namespace:
    """Parse and analyse command line parameters."""

    # define the parameter list
    parser = argparse.ArgumentParser(description="Plot KEGG annotations", usage='%(prog)s -i <input_SET> -o <OUTPUT_DIRECTORY>[options]', prog="plot-kegg-annotations")
    # Mandatory arguments
    parser.add_argument("-i", "--brite-list", type=str, required=True, help="File containing extrated KEGG Brite information for each protein annotated with Kofamscan. It is commonly called 'master.brite.list.tsv'")
    parser.add_argument("-o", "--output-dir", type=str, required=True, help="The directory in which the output will be stored.", default="")
    parser.add_argument("-m", "--strain-names-mapping", type=str, required=False, help="TSV table file to map the original proteome file names to strain names used in the manuscript.", default="")
    parser.add_argument("-e", "--file-extension", type=str, required=False, help="Extension used for the output plots. Default=svg", default="svg")
    parser.add_argument("-t", "--threads", type=int, required=False, help="Maximum number of CPUs to be used. Default=8", default=8)
    parser.add_argument("-wl", "--write-logs", required=False, help="Write logs to file.", default=False, action="store_true")
    parser.add_argument("-d", "--debug", required=False, help="Show debug lines.", default=False, action="store_true")

    # parse the arguments
    args = parser.parse_args()

    return args



def count_hits_for_kentries(df: pl.DataFrame, outpath: Path) -> None:
    """Group and count hits for each KEGG object and for each species in the main brite table"""

    dbgstr: str = f"""count_hits_for_kentries :: START
    Columns: {df.columns}
    Shape: {df.shape}
    Output table:\t{outpath}
    """
    logger.debug(dbgstr)

    # Groups by 4 columns
    grp_results: pl.GroupBy = df.group_by(pl.col("Species"), pl.col("Brite-lev3"), pl.col("Brite-lev2"), pl.col("KEGG-Entry"), maintain_order=True)
    # Count the entries in the groups
    cntdf: pl.DataFrame = grp_results.agg(pl.len().alias("KEGG-entry-hits"))
    del grp_results
    # Sort columns
    cntdf = cntdf.sort(["Species", "KEGG-Entry"])

    print(cntdf.head(n=5))
    print(cntdf.tail(n=3))

    # Compute the total number of hits for a given KObject
    cntdf_unique: pl.DataFrame = cntdf.unique(subset=["Species", "KEGG-Entry"], keep="any", maintain_order=True)
    print(cntdf_unique.head(n=5))
    print(cntdf["KEGG-Entry"].sort().unique().len())
    # This dictionary contains:
    # total number of k-objects with at least one hits
    # Cumulative sum of k-obj hits (each k-object could have multiple hits)
    perspecies_total_hits: dict[str, tuple[int, int]] = {}
    sp_subset: pl.DataFrame = pl.DataFrame([], schema=["a"])
    # Count uniq hits and totals
    for sp in list(cntdf.select(pl.col("Species")).unique().to_series()):
        sp_subset = cntdf_unique.filter(pl.col("Species") == sp)
        perspecies_total_hits[sp] = (sp_subset["KEGG-entry-hits"].len(), int(sp_subset["KEGG-entry-hits"].sum()))

    del cntdf_unique

    # Create the final df which includes the percentages
    df_final: pl.DataFrame = pl.DataFrame([], schema=cntdf.schema)
    # add the column with percentage
    df_final.insert_column(df_final.shape[1], pl.Series("KEGG-entry-hits-pct", [], dtype=pl.Float64))
    print(df_final.head())

    tmp_total: int = 0

    # Compute the percentage of hits for each species
    for sp, v in perspecies_total_hits.items():
        sp_subset = cntdf.filter(pl.col("Species") == sp)
        tmp_total = v[1]
        sp_subset = sp_subset.with_columns(((pl.col("KEGG-entry-hits") / tmp_total) * 100.).alias("KEGG-entry-hits-pct"))
        # print(sp_subset.head(n=5))
        # print(sp_subset.sample(n=10))
        # Add records to the df
        df_final = pl.concat([df_final, sp_subset])

        # break
    del sp_subset

    # Write dataframe to TSV file
    df_final = df_final.sort(["Species", "KEGG-Entry", "KEGG-entry-hits", "Brite-lev2", "Brite-lev3"])
    df_final.write_csv(outpath, separator="\t")



def generate_tables_for_heatmaps(df: pl.DataFrame, outdir: Path, brite_l2_hierarchies: list[str] = ["Amino acid metabolism"], val_cols: list[str] = ["KEGG-entry-hits-pct", "KEGG-entry-hits"]) -> list[Path]:
    """Generate the datapoint tables by pivoting the inpu dataframe"""

    dbgstr: str = f"""generate_tables_for_heatmaps :: START
    Columns: {df.columns}
    Shape: {df.shape}
    Output table:\t{outdir}
    Brite lev2 hirarchies: {brite_l2_hierarchies}
    Value columns: {val_cols}
    """
    logger.debug(dbgstr)

    # Temporary variables
    df_brite_l2: pl.DataFrame = pl.DataFrame([], schema=[("dummy", pl.UInt32)])
    df_brite_l2_pivoted: pl.DataFrame = pl.DataFrame([], schema=[("dummy_pivoted", pl.UInt32)])
    outname_p1: str = ""
    outname_p2: str = ""
    outpath: Path = Path()
    outpaths: list[Path] = []

    # Iterate through brite_l2_hierarchies
    # And generate a Pivoted table for each pair (brite_l2_hierarchies[i], val_cols[j])
    for bl2 in brite_l2_hierarchies:
        # df_brite_l2 = df.filter(pl.col("Brite-lev2") == "Amino acid metabolism")
        df_brite_l2 = df.filter(pl.col("Brite-lev2") == bl2)
        # print(df_brite_l2.shape)
        outname_p1 = bl2.replace(" ", ".")

        # Now iterate through the value columns
        for val in val_cols:
            # Pivot the dataframe (we aggregate using the sum because for a given Brite-lev2 hierachy there could be multiple entries with different K-objects)
            df_brite_l2_pivoted = df_brite_l2.pivot("Species", index="Brite-lev3", values=val, aggregate_function="sum")
            df_brite_l2_pivoted = df_brite_l2_pivoted.sort(["Brite-lev3"])
            # Fill cells that have no value with zeros
            df_brite_l2_pivoted = df_brite_l2_pivoted.fill_null(strategy="zero")
            # print(df_brite_l2_pivoted.head(10))

            # Add the values column name to the output name
            outname_p2 = f"_{val}.tsv"
            outpath = outdir.joinpath(f"{outname_p1}{outname_p2}")
            # print(outname_p2)
            # print(outpath)
            # Write the output TSV file
            df_brite_l2_pivoted.write_csv(outpath, separator="\t")
            # Insert the file path to the output list
            outpaths.append(outpath)

    return outpaths



def load_mapping_names(mapping_tbl: Path, old_names_col: str, new_names_col: str) -> dict[str, str]:
    """Generate the datapoint tables by pivoting the inpu dataframe"""

    dbgstr: str = f"""load_mapping_names :: START
    Mapping file: {mapping_tbl}
    Mapping file:\t{old_names_col}
    Mapping file:\t{new_names_col}
    """
    logger.debug(dbgstr)

    # The mapping file should have a format similar to the one below
    # ID      Taxonomy        RAST_name
    # Ac8     Acinetobacter guillouiae        RAST_ATAC8
    # Br13    Brevundimonas sp.       RAST_ATBR13

    # load the tables
    df: pl.DataFrame = pl.read_csv(source=mapping_tbl, columns=[old_names_col, new_names_col], separator="\t", has_header=True)

    return dict(zip(df[old_names_col].to_list(), df[new_names_col].to_list()))



def plot_heatmap(tbl_file: Path, plots_dir: Path, excluded_sp: list[str], rename_sp: dict[str, str], fmt: str) -> Path:

    dbgstr: str = f"""plot_heatmap :: START
    Datapoints file: {tbl_file}
    Output directory: {plots_dir}
    Excluded species: {excluded_sp}
    Rename species: {rename_sp}
    Output format:\t{fmt}
    """
    logger.debug(dbgstr)

    df: pl.DataFrame = pl.read_csv(tbl_file, separator="\t")
    # Identify the indexes of the columns with the numberic values
    datapoints: pl.DataFrame = df.select(pl.nth(range(1, (df.shape[1]))))
    # remove species if needed
    if len(excluded_sp) > 0:
        datapoints = datapoints.drop(excluded_sp, strict=True)
    # Rename the column names if required
    if len(rename_sp) > 0:
        datapoints = datapoints.rename(rename_sp)
    # sns.color_palette("mako", as_cmap=True)

    # Set the output file name using the input datapoint file name
    # Input files have the following pattern: Amino.acid.metabolism_KEGG-entry-hits-pct.tsv
    basename: str = os.path.basename(tbl_file)
    outfname: str = basename.rsplit("_KEGG-", maxsplit=1)[0]
    plot_title: str = outfname.replace(".", " ")
    # Now add the Hits or '%' to the file names and titles
    if basename.endswith("pct.tsv"):
        plot_title = f"{plot_title} (%)"
        outfname = f"{outfname}.pct.{fmt}"
    else:
        plot_title = f"{plot_title} (hits)"
        outfname = f"{outfname}.hits.{fmt}"
    outpath: Path = plots_dir.joinpath(outfname)
    # Generate the plot
    heatmap = sns.heatmap(datapoints, xticklabels=datapoints.columns, yticklabels=df["Brite-lev3"], annot=False, vmin=0, cmap=sns.cm.rocket_r)
    # Set title, x and y labels
    heatmap.set_title(plot_title)
    plt.xlabel("Species")
    plt.ylabel("Pathway (Brite Level 3)")
    # Save the plot
    plt.savefig(outpath, bbox_inches="tight", transparent=True, dpi=300, format=fmt)
    # This is to avoid that plots get overwritten
    plt.close()
    # sys.exit("DEBUG :: plot_heatmaps")

    return outpath




#####  MAIN  #####
def main():
    """Main function executing SonicParanoid"""

    # start setting the needed variables
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

    # file with BRITE information
    # This file contains all information extrated from KEGG files
    # for each annotated input proteins.
    # Its input should follow the following structure
    # Species Brite-lev3 Brite-lev2 Brite-lev1 KEGG-Name KEGG-Symbol KEGG-Entry Protein_sequence
    # RAST_ATAC8.fasta Two-component system Signal transduction Environmental Information Processing chromosomal replication initiator protein dnaA K02313 ATAC8_0001
    # RAST_ATAC8.fasta Cell cycle - Caulobacter Cell growth and death Cellular Processes chromosomal replication initiator protein dnaA K02313 ATAC8_0001
    # RAST_ATAC8.fasta DNA replication Replication and repair Genetic Information Processing DNA polymerase III subunit beta [EC:2.7.7.7] dnaN K02338 ATAC8_0002
    # NOTE: this file is usually called 'master.brite.list.tsv'
    brite_master_file: Path = Path(args.brite_list)
    # proteomes_dir: Path = Path(args.input_proteomes)
    # Create the output directories
    outdir: Path = Path(args.output_dir)
    ct.systools.makedir(outdir)
    datapoints_dir: Path = outdir.joinpath("datapoints")
    ct.systools.makedir(datapoints_dir)
    plots_dir: Path = outdir.joinpath("plots")
    ct.systools.makedir(plots_dir)
    strain_names_mapping_tbl: Path = Path(args.strain_names_mapping)
    ext: str = args.file_extension
    threads: uint32 = uint32(args.threads)
    writelog: bool = args.write_logs

    # Show some info
    logger.info(f"""The plots will be created using the following parameters:
    Master BRITE list file: {brite_master_file}
    Output directory: {outdir}
    Datapoints directory: {datapoints_dir}
    Plots directory: {plots_dir}
    Strain names mapping file: {strain_names_mapping_tbl}
    File extension for plots:\t{ext}
    Write logs:\t{writelog}
    CPUs to be used:\t{threads}
    Debug mode:\t{debug}
    """)

    # load the dataframe with the master BRITE list
    df: pl.DataFrame = pl.read_csv(brite_master_file, separator="\t")
    # Output path for the final table (could be used to generate the plots)
    brite_counts_tsv: Path = datapoints_dir.joinpath("master.brite.counts.tsv")
    # Count hits for K-entries
    print(brite_counts_tsv)
    count_hits_for_kentries(df, outpath=brite_counts_tsv)
    # Load df with counts
    df: pl.DataFrame = pl.read_csv(brite_counts_tsv, separator="\t")

    # Generate tables for heatmaps
    # Brite lev2 hierachies and value columns on which to perform the Pivots
    brite_l2_hierachies: list[str] = ["Amino acid metabolism", "Biosynthesis of other secondary metabolites", "Carbohydrate metabolism", "Cell motility", "Cellular community - prokaryotes", "Energy metabolism", "Xenobiotics biodegradation and metabolism"]
    val_cols: list[str] = ["KEGG-entry-hits-pct", "KEGG-entry-hits"]

    # generate heatmap datapoints
    heatmap_datapoint_files: list[Path] = generate_tables_for_heatmaps(df=df, outdir=datapoints_dir, brite_l2_hierarchies=brite_l2_hierachies, val_cols=val_cols)

    '''
    # Original Names based on prodigal predictions
    # Species were renamed under Nhi's suggestion as follows:
    # Poryzihabitans.1 -> Ps3
    # Shigellasp.2 -> Sh6
    # Aguillouiae.4 -> Ac8
    # Paraburkholderia.5 -> Pam10
    # Kaureofaciens.7 -> Ki12
    # Brevundimonassp.8 -> Br13
    # Paeruginosa.11 -> Psa14
    # The following 3 species will not be inlcuded in the manuscript
    # Aguillouiae.3
    # Brevundimonassp.9
    # Rpickettii.6
    ################

    sp_mapping: dict[str, str] = {
        "Poryzihabitans.1" : "Ps3",
        "Shigellasp.2" : "Sh6",
        "Aguillouiae.4" : "Ac8",
        "Paraburkholderia.5" : "Pam10",
        "Kaureofaciens.7" : "Ki12",
        "Brevundimonassp.8" : "Br13",
        "Paeruginosa.11" : "Psa14",
    }
    '''


    # Name mapping for annotation based on proteins obtained from RAST predictions
    # sp_mapping: dict[str, str] = {
    #     "GCF_001580545.1_ASM158054v1_protein" : "Parabulk. Ref. (Rast)",
    #     "PAM10_aa" : "Pam10 (Rast)",
    # }

    # strain mapping
    sp_mapping: dict[str, str] = {}
    # Load the mapping file if any
    if len(args.strain_names_mapping) > 0:
        sp_mapping = load_mapping_names(strain_names_mapping_tbl, old_names_col="RAST_file_name", new_names_col="ID")
    else:
        print("no mapping file!")

    # The following species will not be included in the manuscript
    # excluded_sp: list[str] = ["RAST_ATAC8", "RAST_ATSH6"]
    excluded_sp: list[str] = []

    for x in heatmap_datapoint_files:
        # print(x)
        plot_heatmap(tbl_file=x, plots_dir=plots_dir, excluded_sp=excluded_sp, rename_sp=sp_mapping, fmt=ext)

    # sys.exit("DEBUG :: main")


if __name__ == "__main__":
    main()
