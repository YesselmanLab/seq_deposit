"""
command line interface for seq_deposit script
"""

import os
import shutil
import click
import pandas as pd
import yaml
from pathlib import Path

from seq_tools import has_t7_promoter

from seq_deposit.logger import setup_applevel_logger, get_logger
from seq_deposit.deposit import deposit_files

from seq_deposit.gsheets import (
    fetch_sequence_gsheet,
    get_construct_entry,
    get_last_codes,
)
from seq_deposit.settings import LIB_PATH


def get_next_code(code: str) -> str:
    """
    generates the next construct code
    :param code: the previous code
    :return: next construct code
    """
    order = "0123456789ABCDEFGHIJKLMNPQRSTUVWXYZ"
    comps = list(code)
    pos = len(comps) - 1
    while pos > 1:
        order_pos = order.index(comps[pos])
        if order_pos == len(order) - 1:
            pos -= 1
            continue
        comps[pos] = order[order_pos + 1]
        for i in range(pos + 1, len(comps)):
            comps[i] = "1"
        break
    return "".join(comps)


def get_params():
    """
    load parameters from yaml file
    :return:
    """
    # get path relative to this file + resources/params.yaml
    path = os.path.join(os.path.dirname(__file__), "resources", "params.yml")
    with open(path, "r", encoding="utf-8") as f:
        params = yaml.safe_load(f)
    params["deposit_path"] = os.path.abspath(params["deposit_path"])
    return params


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
    log = get_logger("log_constructs")
    log.info(
        "Add seq-deposit-output/constructs.csv to sequences.gsheet on all sheet"
    )
    with open("seq-deposit-output/constructs.csv", "w") as f:
        f.write(constructs[0].to_csv_header() + "\n")
        for c in constructs:
            f.write(c.to_csv_str() + "\n")


def log_primers(primers):
    """
    logs the new primers
    :param constructs: the new primers
    :return: None
    """
    log = get_logger("log_primers")
    log.info(
        "seq-deposit-output/primers.csv to the primers.gsheet on oligo sheet"
    )
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
@click.argument("csvs", nargs=-1)
def opools(csvs):
    """
    generate required files for opools
    :param csvs: list of csv files
    """
    setup_applevel_logger()
    log = get_logger("opools")
    params = get_params()
    df_seqs = fetch_sequence_gsheet(params)
    last_code = df_seqs.iloc[-1]["code"]
    log.info("last code on sheet: %s", last_code)
    data = []
    codes = []
    for csv in csvs:
        df = pd.read_csv(csv)
        log.info("processing %s it has %d sequences", csv, len(df))
        if not has_t7_promoter(df):
            raise ValueError(f"no t7 promoter found in sequences in {csv}")
        code = get_next_code(last_code)
        last_code = code
    for row in data:
        print(",".join([str(x) for x in row]))
    log_new_libaries(params, csvs, codes)


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
    """
    setup_applevel_logger()
    log = get_logger("assembly")
    log.info("outputs will be written in seq-deposit-output")
    if os.path.exists("seq-deposit-output") and not overwrite:
        raise ValueError("seq-deposit-output already exists use --overwrite")
    os.makedirs("seq-deposit-output", exist_ok=True)
    params = get_params()
    last_code, last_primer_code = get_last_codes(params)
    df_constructs = pd.read_csv(construct_csv)
    df_primers = pd.read_csv(primer_csv)
    log.info(
        "processing %s it has %d sequences", construct_csv, len(df_constructs)
    )
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
        centry = get_construct_entry(
            df_row, row["name"], code, ignore_missing_t7
        )
        last_code = code
        constructs.append(centry)
    log_constructs(constructs)
    log_primers(all_primers)
    df_primers["primer_code"] = df_primers["primer_name"].apply(
        lambda x: p_codes[x]
    )
    df_primers["code"] = df_primers["name"].apply(lambda x: codes[x])
    df_primers = df_primers[
        ["name", "code", "primer_num", "primer_name", "primer_code"]
    ]
    df_primers.to_csv("seq-deposit-output/assemblies.csv", index=False)


if __name__ == "__main__":
    cli()
