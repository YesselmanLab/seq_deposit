"""
command line interface for seq_deposit script
"""

import os
import shutil
import click
import pandas as pd
import yaml
from pathlib import Path
from typing import List

from seq_deposit.logger import setup_logging, get_logger
from seq_deposit.deposit import deposit_files

from seq_deposit.gsheets import (
    get_construct_entry,
    get_last_codes,
)

from gsheets.sheet import get_next_code

log = get_logger(__name__)


# helper functions ############################################################


def setup(ignore_missing_t7: bool, overwrite: bool) -> None:
    """
    Set up the necessary configurations for seq-deposit.

    Args:
        ignore_missing_t7 (bool): Flag indicating whether to ignore missing t7 promoter.
        overwrite (bool): Flag indicating whether to overwrite existing output directory.

    Returns:
        None
    """
    # setup logger
    setup_logging()
    # create output directory
    if os.path.exists("seq-deposit-output") and not overwrite:
        log.error("seq-deposit-output already exists use --overwrite")
        exit(1)
    log.info("outputs will be written in seq-deposit-output")
    if not os.path.exists("seq-deposit-output"):
        os.makedirs("seq-deposit-output")
    if ignore_missing_t7:
        log.warning("ignore missing t7 promoter. This is not recommended!!")


# logging #####################################################################
# these functions are to create informationa about the new sequences


def log_new_libaries(params, new_libs, codes):
    """
    logs the new libraries
    :param params: module params
    :param new_libs: the new libraries
    :return: None
    """
    path = os.path.join(params["deposit_path"], "org_data_csvs", "log.csv")
    df_log = pd.read_csv(path)
    # create backup in case something goes wrong
    df_log.to_csv(
        os.path.join(params["deposit_path"], "org_data_csvs", "log_old.csv"),
        index=False,
    )
    pos = len(df_log) - 1
    for lib, code in zip(new_libs, codes):
        df_log.loc[pos] = [
            Path(lib).stem,
            code,
            os.path.abspath(lib),
        ]
        shutil.copy(
            lib,
            os.path.join(
                params["deposit_path"], "org_data_csvs", "csvs", f"{code}.csv"
            ),
        )
        pos += 1
    df_log.to_csv(path, index=False)


def log_constructs(constructs):
    """
    logs the new constructs
    :param constructs: the new constructs
    :return: None
    """
    log.info("Add seq-deposit-output/constructs.csv to sequences.gsheet on all sheet")
    with open("seq-deposit-output/constructs.csv", "w") as f:
        f.write(constructs[0].to_csv_header() + "\n")
        for c in constructs:
            f.write(c.to_csv_str() + "\n")


def log_primers(primers):
    """
    logs the new primers
    :param primers: the new primers
    :return: None
    """
    log.info("seq-deposit-output/primers.csv to the primers.gsheet on oligo sheet")
    with open("seq-deposit-output/primers.csv", "w") as f:
        f.write("name,code,c_code,useable,a_temp,sequence,type\n")
        for p in primers:
            f.write(",".join([str(x) for x in p]) + "\n")


# cli functions ###############################################################


@click.group()
def cli():
    """
    command line interface for seq_deposit script
    """


@cli.command(help="generate required files for opools")
@click.option("--ignore-missing-t7", is_flag=True, help="ignore t7 promoter")
@click.option("--ignore-missing-rt-seq", is_flag=True, help="ignore rt seq")
@click.option("--overwrite", is_flag=True, help="overwrite existing files")
@click.option(
    "--ntype", default="DNA", help="type of nucleic acid", allowed=["DNA", "RNA"]
)
@click.option("--t7-seq", default="TTCTAATACGACTCACTATA", help="t7 promoter sequence")
@click.argument("csvs", nargs=-1)
def opools(
    csvs: List[str],
    ntype: str,
    ignore_missing_t7: bool,
    ignore_missing_rt_seq: bool,
    t7_seq: str,
    overwrite: bool,
) -> None:
    """
    Generate required files for opools.

    Args:
        csvs (List[str]): List of CSV file paths.
        ignore_missing_t7 (bool): Flag to ignore missing t7 promoter.
        ignore_missing_rt_seq (bool): Flag to ignore rt seq.
        t7_seq (str): T7 promoter sequence.
        overwrite (bool): Flag to overwrite existing files.
    """
    setup(ignore_missing_t7, overwrite)
    last_code, _ = get_last_codes()
    constructs = []
    codes = []
    for csv in csvs:
        if ntype == "DNA":
            df_dna = pd.read_csv(csv)

            pass
        log.info("processing %s it has %d sequences", csv, len(df_dna))
        exit()
        code = get_next_code(last_code)
        centry = get_construct_entry(df, Path(csv).stem, code, ignore_missing_t7)
        last_code = code

    # log_constructs(constructs)
    # if not dry_run:
    #    log.info("updating log csv ")
    #    log_new_libaries(params, csvs, codes)


@cli.command(help="generate required files for primer assemblies")
@click.argument("construct_csv")
@click.argument("primer_csv")
@click.option("--ignore-missing-t7", is_flag=True, help="ignore t7 promoter")
@click.option("--dry-run", is_flag=True, help="do not write any files")
@click.option("--overwrite", is_flag=True, help="overwrite existing files")
def assembly(construct_csv, primer_csv, ignore_missing_t7, dry_run, overwrite):
    """
    generate required files for primer assemblies
    :param construct_csv: path to construct csv
    :param primer_csv: path to primer csv
    :param ignore_missing_t7: ignore missing t7 promoter
    :param dry_run: do not write any files
    :param overwrite: overwrite existing files
    """
    setup(ignore_missing_t7, dry_run, overwrite)
    params = get_params()
    last_code, last_primer_code = get_last_codes(params)
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
    log_constructs(constructs)
    log_primers(all_primers)
    df_primers["primer_code"] = df_primers["primer_name"].apply(lambda x: p_codes[x])
    df_primers["code"] = df_primers["name"].apply(lambda x: codes[x])
    df_primers = df_primers[
        ["name", "code", "primer_num", "primer_name", "primer_code"]
    ]
    df_primers.to_csv("seq-deposit-output/assemblies.csv", index=False)
    if not dry_run:
        log.info("updating log csv ")
        log_new_libaries(params, [construct_csv], [";".join(codes)])


@cli.command(help="update deposit path")
@click.argument("path")
def set_deposit_path(path):
    """
    update deposit path
    :param path: path to deposit
    """
    params = get_params()
    if not os.path.exists(path):
        raise ValueError(f"path {path} does not exist")
    params["deposit_path"] = path
    path = os.path.join(os.path.dirname(__file__), "resources", "params.yml")
    with open(path, "w", encoding="utf-8") as f:
        yaml.dump(params, f, default_flow_style=False, allow_unicode=True)


@cli.command(help="update libraries")
@click.argument("construct_csv")
@click.argument("csvs", nargs=-1)
def update_libraries(construct_csv, csvs):
    params = get_params()
    df_construct = pd.read_csv(construct_csv)
    for i, row in df_construct.iterrows():
        name = row["name"]
        current_csv = None
        for csv in csvs:
            if name in csv:
                current_csv = csv
                break
        df = pd.read_csv(current_csv)
        log.info("processing %s it has %d sequences", csv, len(df))
        if not has_t7_promoter(df):
            raise ValueError(f"no t7 promoter found in sequences in {csv}")
        code = row["code"]
        deposit_files(df, code, params["deposit_path"], False, False)


@cli.command(help="get last code")
def get_last_code():
    params = get_params()
    df_seqs = fetch_sequence_gsheet(params)
    df_primers = fetch_primers_gsheet(params)
    print(df_seqs.iloc[-1])
    print(df_primers.iloc[-1])


# TODO takes a seq-desposit directory and actually transfers files!
@cli.command()
def deposit():
    # deposit_files(df, code, params["deposit_path"], dry_run, ignore_missing_t7)
    pass


if __name__ == "__main__":
    cli()
