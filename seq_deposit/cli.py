"""
command line interface for seq_deposit script
"""

import os
import shutil
import cloup
import pandas as pd
import sys
import glob
from pathlib import Path
from typing import List, Dict, Any
from tabulate import tabulate

from seq_tools import to_fasta, to_dna, to_rna, add

from gsheets.sheet import get_next_code, get_sequence_sheet, get_oligo_sheet

from seq_deposit.management import (
    get_rna_dataframe_from_row,
    get_dna_dataframe_from_row,
    fasta_to_dataframe,
)
from seq_deposit.logger import setup_logging, get_logger
from seq_deposit.prepare import generate_rna_dataframe
from seq_deposit.process import process_constructs
from seq_deposit.util import create_directories

log = get_logger("cli")


# setup functions ############################################################


def log_initial_info(args: Dict[str, Any]) -> None:
    """Log the initial configuration details.

    Args:
        args (Dict[str, Any]): Dictionary containing setup arguments.

    Returns:
        None
    """
    log.info("Initial configuration:")
    for key, value in args.items():
        log.info("%s: %s", key, value)
    log.info("Ran at commandline as: %s", " ".join(sys.argv))


def handle_missing_sequences(args: Dict[str, Any]) -> None:
    """Handle the logic for missing T7 and RT sequences.

    Args:
        args (Dict[str, Any]): Dictionary containing setup arguments.

    Returns:
        None
    """
    if args.get("ignore_missing_t7"):
        log.warning("Ignore missing T7 promoter. This is not recommended!!")
    else:
        log.info("T7 promoter sequence: %s", args.get("t7_seq"))
    if args.get("ignore_missing_rt_seq"):
        log.warning("Ignore missing RT sequence. This is not recommended!!")
    else:
        log.info("RT sequence: %s", args.get("rt_seq"))


def setup(
    args: Dict[str, Any],
) -> None:
    """Set up the necessary configurations and directories for seq-deposit.

    Args:
        args (Dict[str, bool or str]): Dictionary containing setup arguments.

    Returns:
        None
    """
    if os.path.exists("seq-deposit-output") and not args["overwrite"]:
        print("seq-deposit-output directory exists")
        exit()
    setup_logging(file_name="seq-deposit-output/seq-deposit.log")
    log_initial_info(args)
    create_directories()
    handle_missing_sequences(args)


# cli functions ###############################################################
def common_options() -> cloup.OptionGroup:
    """
    Returns a Cloup option group with common options.

    Returns:
        cloup.OptionGroup: A Cloup option group containing the common options.

    Options:
        --ignore-missing-t7: Ignore the t7 promoter. (default: False)
        --ignore-missing-rt-seq: Ignore the rt sequence. (default: False)
        --overwrite: Overwrite existing files. (default: False)
        --t7-seq: The t7 promoter sequence. (default: "TTCTAATACGACTCACTATA")
        --rt-seq: The rt sequence. (default: "AAAGAAACAACAACAACAAC")
        --threads: The number of threads to use. (default: 1)
    """
    return cloup.option_group(
        "common options",
        cloup.option("--ignore-missing-t7", is_flag=True, help="ignore t7 promoter"),
        cloup.option("--ignore-missing-rt-seq", is_flag=True, help="ignore rt seq"),
        cloup.option("--overwrite", is_flag=True, help="overwrite existing files"),
        cloup.option(
            "--t7-seq", default="TTCTAATACGACTCACTATA", help="t7 promoter sequence"
        ),
        cloup.option("--rt-seq", default="AAAGAAACAACAACAACAAC", help="rt sequence"),
        cloup.option("--threads", default=1, help="number of threads to use"),
    )


@cloup.group()
def cli():
    """
    command line interface for seq_deposit script
    """
    pass


@cli.command(help="")
@common_options()
def from_final_dir(**kwargs):
    setup(kwargs)
    if os.path.isdir("final"):
        log.info("final directory exists")
    else:
        log.error("final directory does not exist")
        exit()
    df = pd.read_json("final/summary_jsons/without_codes.json")
    process_constructs(df)


@cli.command(help="generate required files for opools")
@cloup.option("--ntype", default="DNA", help="type of nucleic acid")
@common_options()
@cloup.argument("csvs", nargs=-1)
def opools(
    csvs: List[str],
    ntype: str,
    **kwargs,
) -> None:
    """
    Generate required files for opools.

    Args:
        csvs (List[str]): List of CSV file paths.
        ntype (str): Type of nucleic acid. Default is "DNA".
    Returns:
        None
    """
    setup(kwargs)
    log.info("processing %d csvs", len(csvs))
    log.info("ntype: %s", ntype)
    del kwargs["overwrite"]
    del kwargs["rt_seq"]
    if len(csvs) == 0:
        log.error("no csvs supplied")
        return
    dfs = []
    for csv in csvs:
        df = pd.read_csv(csv)
        df["construct"] = Path(csv).stem
        df["type"] = "OPOOL"
        dfs.append(df)
    df = pd.concat(dfs)
    if ntype == "DNA":
        df_rna = generate_rna_dataframe(df, ntype, **kwargs)
        if df_rna is None:
            log.error("rna generation failed")
            return
        df_rna["dna_sequence"] = df["sequence"]
        df = df_rna
    else:
        df_dna = to_dna(df)
        df_dna = add(df_dna, kwargs["t7_seq"], "")
        if "structure" or "ens_defect" not in df.columns:
            df = generate_rna_dataframe(df, ntype, **kwargs)
        df["dna_sequence"] = df_dna["sequence"]
    process_constructs(df)


@cli.command(help="generate required files for primer assemblies")
@cloup.argument("results_json", type=cloup.Path(exists=True))
@common_options()
def assembly(
    results_json: str,
    **kwargs,
):
    setup(kwargs)
    del kwargs["overwrite"]
    del kwargs["rt_seq"]
    df = pd.read_json(results_json)
    df_rna = generate_rna_dataframe(df, "DNA", **kwargs)
    if df_rna is None:
        log.error("rna generation failed")
        return
    df_rna["dna_sequence"] = df["sequence"]
    df_rna["construct"] = df["name"]
    df_rna["type"] = "ASSEMBLY"
    process_constructs(df_rna)


# TODO search by name size etc return csv with results
@cli.command()
@cloup.option("--code", default=None, help="code of the sequence")
def get_sequence_info(code):
    setup_logging()
    df = get_sequence_sheet()
    if code is not None:
        df = df[df["code"] == code]
        if len(df) == 0:
            log.error(f"code {code} not found")
            return
    log.info("num of constructs %d", len(df))
    seq_path = os.environ["SEQPATH"]
    for i, row in df.iterrows():
        if row["size"] == 1:
            print("name: ", row["name"])
            print("sequence: ", row["rna_sequence"])
            print("structure: ", row["rna_structure"])
        else:
            df_rna = pd.read_csv(f"{seq_path}/rna/{row['code']}.csv")
            rna_row = df_rna.iloc[0]
            print("example sequence ->")
            print("name: ", rna_row["name"])
            print("rna sequence: ", rna_row["sequence"])
            print("rna structure: ", rna_row["structure"])
        print("code: ", row["code"])
        print("type: ", row["type"])
        print("-" * 40)


@cli.command()
@cloup.option("--deposit-path", default=None, help="deposit path")
@cloup.option("--force", is_flag=True, help="force copy of files")
def deposit(deposit_path, force):
    # check if the deposit path exists if not use supplied option
    setup_logging()
    output = input(
        "are you sure you want to deposit all files (i.e. did you check they are corect)? [y/n]"
    )
    if output == "y":
        print("depositing files")
    else:
        return
    if not os.path.exists("seq-deposit-output"):
        log.error("seq-deposit-output does not exist")
        return
    if deposit_path is None:
        deposit_path = os.environ.get("SEQPATH")
    if deposit_path is None:
        log.error("no deposit path supplied")
        return
    dna_csvs = glob.glob("seq-deposit-output/dna/*.csv")
    for csv in dna_csvs:
        if os.path.isfile(f"{deposit_path}/dna/{csv}") and not force:
            log.error(f"file exists: {csv} in DNA use --force to overwrite")
            return
        shutil.copy(csv, f"{deposit_path}/dna/")
    rna_csvs = glob.glob("seq-deposit-output/rna/*.csv")
    for csv in rna_csvs:
        if os.path.isfile(f"{deposit_path}/rna/{csv}") and not force:
            log.error(f"file exists: {csv} in RNA use --force to overwrite")
            return
        shutil.copy(csv, f"{deposit_path}/rna/")
    fasta_files = glob.glob("seq-deposit-output/fastas/*.fasta")
    for fasta in fasta_files:
        if os.path.isfile(f"{deposit_path}/rna/{fasta}") and not force:
            log.error(f"file exists: {csv} in fasta use --force to overwrite")
            return
        shutil.copy(fasta, f"{deposit_path}/fastas/")


@cli.command()
@cloup.argument("code")
@cloup.option("--deposit-path", default=None, help="deposit path")
def generate_files(code, deposit_path):
    setup_logging()
    if deposit_path is None:
        deposit_path = os.environ.get("SEQPATH")
    if deposit_path is None:
        log.error("no deposit path supplied")
        return
    df = get_sequence_sheet()
    df = df[df["code"] == code]
    if len(df) == 0:
        log.error(f"code {code} not found")
        return
    if len(df) > 1:
        log.error("more than one row found")
        return
    if os.path.isfile(f"{deposit_path}/rna/{code}.csv"):
        log.error(f"file exists: {code}.csv in RNA")
        return
    if os.path.isfile(f"{deposit_path}/dna/{code}.csv"):
        log.error(f"file exists: {code}.csv in DNA")
        return
    if os.path.isfile(f"{deposit_path}/fastas/{code}.fasta"):
        log.error(f"file exists: {code}.fasta in fastas")
        return
    row = df.iloc[0]
    rna_df = get_rna_dataframe_from_row(row, df.columns)
    dna_df = get_dna_dataframe_from_row(row, df.columns)
    fasta_df = to_dna(rna_df)
    rna_df.to_csv(f"{deposit_path}/rna/{code}.csv", index=False)
    dna_df.to_csv(f"{deposit_path}/dna/{code}.csv", index=False)
    to_fasta(fasta_df, f"{deposit_path}/fastas/{code}.fasta")


if __name__ == "__main__":
    cli()
