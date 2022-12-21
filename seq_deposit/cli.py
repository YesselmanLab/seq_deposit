"""
command line interface for seq_deposit script
"""

import os
import shutil
import subprocess
import click
import pandas as pd
import yaml
from pathlib import Path

from seq_tools import (
    get_length,
    get_molecular_weight,
    get_default_names,
    get_extinction_coeff,
    has_t7_promoter,
    to_dna,
    to_fasta,
    transcribe,
)
from seq_deposit.logger import setup_applevel_logger, get_logger
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


def get_seq_fwd_primer_code(df: pd.DataFrame) -> str:
    """
    gets the sequence forward primer code
    :param df: the dataframe with sequences
    """
    df = df.copy()
    df = to_dna(df)
    path = os.path.join(LIB_PATH, "resources", "p5_sequences.csv")
    df_p5 = pd.read_csv(path)
    for _, row in df_p5.iterrows():
        # if all sequences in df start with the p5 sequence then return the p5 code
        if all(df["sequence"].str.startswith(row["sequence"])):  # type: ignore
            return row["code"]
    return ""


def generate_dna_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """
    generates the dna dataframe
    :param df:
    :return: the dna dataframe
    """
    df_dna = df.copy()
    if "name" not in df_dna.columns:
        df_dna = get_default_names(df_dna)
    df_dna = to_dna(df_dna)
    df_dna = df_dna[["name", "sequence"]]
    df_dna = get_length(df_dna)
    df_dna = get_molecular_weight(df_dna, "DNA", True)
    df_dna = get_extinction_coeff(df_dna, "DNA", True)
    return df_dna


def generate_rna_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """
    generates the rna dataframe
    :param df:
    :return: the rna dataframe
    """
    df_rna = df.copy()
    df_rna = transcribe(df_rna)
    df_rna = df_rna[["name", "sequence"]]
    df_rna = get_length(df_rna)
    df_rna = get_molecular_weight(df_rna, "RNA", False)
    df_rna = get_extinction_coeff(df_rna, "RNA", False)
    return df_rna


def generate_fasta(df: pd.DataFrame, code, params) -> None:
    """
    generates the fasta file
    :param df: the dataframe with sequences
    :param code: the construct code
    :param params: module params
    :return: None
    """
    df = df.copy()
    df = to_dna(df)
    df = df[["name", "sequence"]]
    log = get_logger("generate_fasta")
    path = os.path.join(params["deposit_path"], "fastas", f"{code}.fasta")
    log.info(f"writing fasta file to {path}")
    to_fasta(df, path)


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


def fetch_sequence_gsheet(params) -> pd.DataFrame:
    """
    fetches the sequence google sheet
    :param params: module params
    :return: dataframe of sequence sheet
    """
    gsheet_path = params["sequence_gsheet_url"]
    subprocess.call(
        f'wget --output-file="logs.csv" "{gsheet_path}" -O "temp.csv"',
        shell=True,
    )
    df = pd.read_csv("temp.csv")
    os.remove("temp.csv")
    os.remove("logs.csv")
    cols = (
        "name,code,type,arrived,usuable,size,dna_len,dna_sequence,rna_len,"
        "rna_sequence,rna_structure,fwd_p,rev_p,rt_p,seq_fwd_p,project,"
        "comment".split(",")
    )
    df = df[cols]
    path = os.path.join(LIB_PATH, "resources/all.csv")
    df.to_csv(path, index=False)
    return df


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
        df_dna = generate_dna_dataframe(df)
        log.info(f"writing dna csv to {params['deposit_path']}/dna/")
        path = os.path.join(params["deposit_path"], "dna", f"{code}.csv")
        df_dna.to_csv(path, index=False)
        df_rna = generate_rna_dataframe(df_dna)
        log.info(f"writing rna csv to {params['deposit_path']}/rna/")
        path = os.path.join(params["deposit_path"], "rna", f"{code}.csv")
        df_rna.to_csv(path, index=False)
        generate_fasta(df_rna, code, params)
        dtype = "OPOOL"
        if len(df) > 100:
            dtype = "AGILENT"
        dna_length = round(df_dna["length"].mean())
        seq_fwd_primer = get_seq_fwd_primer_code(df_rna)
        row = [
            Path(csv).stem,
            code,
            dtype,
            "NO",
            "",
            len(df),
            dna_length,
            "LIBARY",
            dna_length - 20,
            "LIBARY",
            "LIBARY",
            "P001E",
            "P001F",
            "RTB",
            "P000R",
            seq_fwd_primer,
        ]
        codes.append(code)
        last_code = code
        data.append(row)
    for row in data:
        print(",".join([str(x) for x in row]))
    log_new_libaries(params, csvs, codes)


"""
@cli.command()
@click.argument("csv")
@click.argument("assembly-csv")
@click.argument("primers-csv")
@click.option("-sc", "--start-code", required=True)
@click.option("-p", "--path", required=True, default=None, help="path")
def primers(csv, assembly_csv, primers_csv, start_code, path):
    assembly = {}
    all_primers = {}
    df = pd.read_csv(primers_csv)
    for i, row in df.iterrows():
        all_primers[row["sequence"]] = row["name"]
    for l in open(assembly_csv).readlines():
        name, primers_str = l.strip().split(",")
        primers = primers_str.split(";")

    df = pd.read_csv(csv)
"""

if __name__ == "__main__":
    cli()
