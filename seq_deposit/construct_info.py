"""
handles processing of google sheets data
"""

import os
import pandas as pd
from dataclasses import dataclass
from typing import Dict, Any

from seq_tools import to_dna, get_length, trim, has_5p_sequence
import vienna

from gsheets.sheet import get_sequence_sheet, get_oligo_sheet

from seq_deposit.logger import get_logger
from seq_deposit.settings import LIB_PATH


log = get_logger(__name__)


def get_seq_fwd_primer_code(df: pd.DataFrame) -> str:
    """
    Returns the code for the forward primer based on the given DataFrame.

    Args:
        df (pd.DataFrame): The DataFrame containing the sequences.

    Returns:
        str: The code for the forward primer, or an empty string if no match is found.
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


def get_default_construct_entry(name, code, ctype, size):
    """
    Create a dictionary representing a construct entry.

    Args:
        name (str): The name of the construct.
        code (str): The code of the construct.
        ctype (str): The type of the construct.
        size (int): The size of the construct.

    Returns:
        dict: A dictionary representing the construct entry with default values for other fields.
    """
    return {
        "name": name,
        "code": code,
        "type": ctype,
        "size": size,
        "arrived": "NO",
        "usuable": "UNK",
        "dna_len": -1,
        "dna_sequence": "LIBRARY",
        "rna_len": -1,
        "rna_sequence": "LIBRARY",
        "rna_structure": "LIBRARY",
        "fwd_p": "NONE",
        "rev_p": "NONE",
        "rt_p": "RTB",
        "seq_fwd_p": "P000R",
        "seq_rev_p": "P000Y",
        "dir": "",
        "designed_by": "",
        "project": "",
        "comment": "",
    }


def get_construct_entry(
    df_dna: pd.DataFrame,
    df_rna: pd.DataFrame,
    name: str,
    ctype: str,
    code: str,
) -> Dict[str, Any]:
    """
    Get the construct entry for a given DNA and RNA DataFrame.

    Args:
        df_dna (pd.DataFrame): DataFrame containing DNA sequences.
        df_rna (pd.DataFrame): DataFrame containing RNA sequences.
        name (str): Name of the construct.
        ctype (str): Type of the construct.
        code (str): Code of the construct.

    Returns:
        Dict[str, Any]: Construct information dictionary.

    """
    construct_info = get_default_construct_entry(name, ctype, code, len(df_dna))
    df_dna = get_length(df_dna)
    df_rna = get_length(df_rna)
    if len(df_dna) == 1:
        construct_info["dna_sequence"] = df_dna.iloc[0]["sequence"]
        construct_info["rna_sequence"] = df_rna.iloc[0]["sequence"]
    else:
        construct_info["fwd_p"] = "P001E"
        construct_info["rev_p"] = "P001F"
    construct_info["dna_len"] = round(df_dna["length"].mean())
    construct_info["rna_len"] = round(df_rna["length"].mean())
    construct_info["seq_rev_p"] = get_seq_fwd_primer_code(df_rna)
    return construct_info


def get_last_codes():
    """
    Retrieves the last construct code and last primer code from the sequence and oligo
    sheets respectively.

    Returns:
        tuple: A tuple containing the last construct code and last primer code.
    """
    df_seqs = get_sequence_sheet()
    df_primers = get_oligo_sheet()
    last_code = df_seqs["code"].loc[df_seqs["code"].last_valid_index()]
    last_primer_code = df_primers["code"].loc[df_primers["code"].last_valid_index()]
    log.info("last construct code on sheet: %s", last_code)
    log.info("last primer code on sheet: %s", last_primer_code)
    return last_code, last_primer_code
