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

from seq_tools import to_fasta, to_dna

from gsheets.sheet import get_next_code, get_sequence_sheet, get_oligo_sheet

from seq_deposit.logger import setup_logging, get_logger
from seq_deposit.deposit import deposit_files, deposit_rna_csv
from seq_deposit.prepare import generate_dna_dataframe, generate_rna_dataframe
from seq_deposit.construct_info import get_construct_entry


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


def check_existing_output_directory(overwrite: bool) -> None:
    """Check if the output directory already exists and handle the overwrite condition.

    Args:
        overwrite (bool): Flag indicating whether to overwrite existing output directory.

    Returns:
        None
    """
    output_dir = "seq-deposit-output"
    if os.path.exists(output_dir) and not overwrite:
        log.error("%s already exists. Use --overwrite.", output_dir)
        sys.exit(1)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    log.info("Outputs will be written in %s", output_dir)


def create_directories() -> None:
    """Create the necessary directories for the output.

    Returns:
        None
    """
    directories = [
        "seq-deposit-output/dna",
        "seq-deposit-output/rna",
        "seq-deposit-output/fastas",
    ]
    for directory in directories:
        os.makedirs(directory, exist_ok=True)
    log.info("Directories created: %s", directories)


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
    setup_logging()
    log_initial_info(args)
    check_existing_output_directory(args["overwrite"])
    create_directories()
    handle_missing_sequences(args)


# helper functions ############################################################


def process_primers(df):
    last_code, last_primer_code = get_last_codes()
    df["primer_codes"] = [[] for _ in range(len(df))]
    primers = []
    for construct, g in df.groupby("name"):
        for i, row in g.iterrows():
            pcodes = []
            pos = 1
            for pname, pseq in zip(row["primer_names"], row["primers"]):
                pcode = get_next_code(last_primer_code)
                primers.append(
                    [
                        pname,
                        pcode,
                        construct,
                        g.iloc[0]["code"],
                        "UNK",
                        "N/A",
                        pseq,
                        "ASSEMBLY",
                        pos,
                    ]
                )
                pcodes.append(pcode)
                last_primer_code = pcode
                pos += 1
            df.at[i, "primer_codes"] = pcodes
    df_primers = pd.DataFrame(
        primers,
        columns="p_name p_code name code useable a_temp sequence type p_num".split(),
    )
    df_assembly = df_primers[["name", "code", "p_num", "p_name", "p_code"]]
    df_assembly.to_csv("seq-deposit-output/assemblies.csv", index=False)
    df_primers.drop(columns=["p_num", "name"], inplace=True)
    df_primers.rename(columns={"p_name": "name", "code": "c_code"}, inplace=True)
    df_primers.rename(columns={"p_code": "code"}, inplace=True)
    df_primers = df_primers[
        ["name", "code", "c_code", "useable", "a_temp", "sequence", "type"]
    ]
    df_primers.drop_duplicates(subset="name", inplace=True)
    df_primers.to_csv("seq-deposit-output/primers.csv", index=False)
    df_primers = pd.DataFrame(
        primers,
        columns=[
            "p_name",
            "p_code",
            "construct",
            "code",
            "useable",
            "a_temp",
            "sequence",
            "type",
            "p_num",
        ],
    )


def process_constructs(df):
    last_code, last_primer_code = get_last_codes()
    df_seqs = get_sequence_sheet()
    constructs = []
    df["code"] = ""
    code = get_next_code(last_code)
    for construct, group in df.groupby("construct"):
        df_dna = group[["name", "dna_sequence"]].rename(
            columns={"dna_sequence": "sequence"}
        )
        df_rna = group[["name", "sequence", "structure", "mfe", "ens_defect"]].copy()
        log.info("processing %s it has %d sequence(s)", construct, len(df_dna))
        df_dna.to_csv(f"seq-deposit-output/dna/{code}.csv", index=False)
        df_rna.to_csv(f"seq-deposit-output/rna/{code}.csv", index=False)
        construct_entry = get_construct_entry(
            df_dna, df_rna, construct, group.iloc[0]["type"], code
        )
        constructs.append(construct_entry)
        df_fasta = df_rna[["name", "sequence"]]
        df_fasta = to_dna(df_fasta)
        to_fasta(df_fasta, f"seq-deposit-output/fastas/{code}.fasta")
        df.loc[group.index, "code"] = code
        last_code = code
        code = get_next_code(last_code)
    df_constructs = pd.DataFrame(constructs, columns=df_seqs.columns)
    log.info(
        "\n"
        + tabulate(
            df_constructs[["name", "code", "type", "size", "dna_len"]],
            headers="keys",
            tablefmt="psql",
        )
    )
    df_sub = df.query("type == 'ASSEMBLY'")
    if len(df_sub) > 0:
        process_primers(df_sub)
    df_constructs.to_csv("seq-deposit-output/constructs.csv", index=False)


def get_last_codes():
    """
    get the last codes from the google sheet
    :param params: module parameters
    """
    log = get_logger("get_last_codes")
    df_seqs = get_sequence_sheet()
    df_primers = get_oligo_sheet()
    last_code = df_seqs["code"].loc[df_seqs["code"].last_valid_index()]
    last_primer_code = df_primers["code"].loc[df_primers["code"].last_valid_index()]
    log.info("last construct code on sheet: %s", last_code)
    log.info("last primer code on sheet: %s", last_primer_code)
    return last_code, last_primer_code


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
        ignore_missing_t7 (bool): Flag to ignore t7 promoter. Default is False.
        ignore_missing_rt_seq (bool): Flag to ignore rt seq. Default is False.
        t7_seq (str): t7 promoter sequence. Default is "TTCTAATACGACTCACTATA".
        threads (int): Number of threads to use. Default is 1.
        overwrite (bool): Flag to overwrite existing files. Default is False.

    Returns:
        None
    """
    setup(kwargs)
    dfs = []
    for csv in csvs:
        df = pd.read_csv(csv)
        df["construct"] = Path(csv).stem
        df["type"] = "OPOOL"
        dfs.append(df)
    df = pd.concat(dfs)
    if ntype == "DNA":
        df.rename(columns={"sequence": "dna_sequence"}, inplace=True)


@cli.command(help="generate required files for primer assemblies")
@cloup.argument("results_csv")
@common_options()
def assembly(
    construct_csv,
    primer_csv,
    ignore_missing_t7,
    ignore_missing_rt_seq,
    t7_seq,
    overwrite,
):
    setup(ignore_missing_t7, ignore_missing_rt_seq, t7_seq, overwrite)
    last_code, last_primer_code = get_last_codes()
    df_constructs = pd.read_csv(construct_csv)
    df_primers = pd.read_csv(primer_csv)
    log.info("processing %s it has %d sequences", construct_csv, len(df_constructs))
    constructs = []
    p_codes = {}
    codes = {}
    all_primers = []
    for _, row in df_constructs.iterrows():
        df_row = pd.DataFrame([row], columns=df_constructs.columns, index=[0])
        df_primers_sub = df_primers[df_primers["name"] == row["name"]]
        code = get_next_code(last_code)
        codes[row["name"]] = code
        for _, p_row in df_primers_sub.iterrows():
            if p_row["primer_name"] not in p_codes:
                p_codes[p_row["primer_name"]] = get_next_code(last_primer_code)
                last_primer_code = p_codes[p_row["primer_name"]]
                all_primers.append(
                    [
                        p_row["primer_name"],
                        last_primer_code,
                        code,
                        "UNK",
                        "N/A",
                        p_row["primer_sequence"],
                        "ASSEMBLY",
                    ]
                )
        deposit_files(df_row, code, params["deposit_path"], dry_run, ignore_missing_t7)
        centry = get_construct_entry(df_row, row["name"], code, ignore_missing_t7)
        last_code = code
        constructs.append(centry)
    df_primers["primer_code"] = df_primers["primer_name"].apply(lambda x: p_codes[x])
    df_primers["code"] = df_primers["name"].apply(lambda x: codes[x])
    df_primers = df_primers[
        ["name", "code", "primer_num", "primer_name", "primer_code"]
    ]
    df_primers.to_csv("seq-deposit-output/assemblies.csv", index=False)


@cli.command(help="get last code")
def get_last_code():
    params = get_params()
    df_seqs = fetch_sequence_gsheet(params)
    df_primers = fetch_primers_gsheet(params)
    print(df_seqs.iloc[-1])
    print(df_primers.iloc[-1])


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


# TODO takes a seq-desposit directory and actually transfers files!
@cli.command()
@cloup.option("--deposit-path", default=None, help="deposit path")
@cloup.option("--force", is_flag=True, help="force copy of files")
def deposit(deposit_path, force):
    # check if the deposit path exists if not use supplied option
    # deposit_files(df, code, params["deposit_path"], dry_run, ignore_missing_t7)
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


if __name__ == "__main__":
    cli()
