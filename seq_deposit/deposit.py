"""
module for generating the differnt csvs that will be stored for future use
"""

import os
import pandas as pd

from seq_tools import (
    get_length,
    get_molecular_weight,
    get_default_names,
    get_extinction_coeff,
    to_dna,
    to_fasta,
    trim,
    transcribe,
)

from seq_deposit.logger import get_logger

log = get_logger("deposit")


def deposit_dna_csv(df: pd.DataFrame, code: str, deposit_path: str) -> None:
    """
    generates the dna dataframe
    :param df: the dataframe with sequences
    :param code: the construct code
    :param deposit_path: where to deposit the csv file
    :param dry_run: if true, don't actually deposit the file
    """
    log.info(f"writing dna csv to {deposit_path}/dna/")
    path = os.path.join(deposit_path, "dna", f"{code}.csv")
    df.to_csv(path, index=False)


def deposit_rna_csv(
    df: pd.DataFrame,
    code: str,
    deposit_path: str,
) -> None:
    log.info(f"writing rna csv to {deposit_path}/rna/")
    path = os.path.join(deposit_path, "rna", f"{code}.csv")
    df.to_csv(path, index=False)


def deposit_fasta_file(df: pd.DataFrame, code, deposit_path, dry_run=False) -> None:
    """
    generates the fasta file
    :param df: the dataframe with sequences
    :param code: the construct code
    :param deposit_path: where to deposit the fasta file
    :param dry_run: if true, don't actually deposit the file
    :return: None
    """
    df = df.copy()
    df = df[["name", "sequence"]]
    df = trim(df, 20, 0)
    if not dry_run:
        path = os.path.join(deposit_path, "fastas", f"{code}.fasta")
        log.info(f"writing fasta file to {path}")
        to_fasta(df, path)


def deposit_files(
    df: pd.DataFrame, code, deposit_path, dry_run=False, ignore_missing_t7=False
) -> None:
    """
    deposits the files
    :param df: the dataframe with sequences
    :param code: the construct code
    :param deposit_path: where to deposit the files
    :param dry_run: if true, don't actually deposit the file
    """
    deposit_dna_csv(df, code, deposit_path, dry_run)
    deposit_rna_csv(df, code, deposit_path, dry_run, ignore_missing_t7)
    deposit_fasta_file(df, code, deposit_path, dry_run)
