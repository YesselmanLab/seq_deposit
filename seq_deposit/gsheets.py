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


def get_construct_entry(name, code, ctype, size):
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
        "size": -1,
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
        "project": "",
        "comment": "",
    }


def generate_rna_dataframe(df, ignore_missing_t7):
    pass


def get_construct_entry_from_dna(
    df: pd.DataFrame,
    name: str,
    ctype: str,
    code: str,
    no_t7: bool = False,
    t7_seq="TTCTAATACGACTCACTATA",
) -> Dict[str, Any]:
    construct_info = get_construct_entry(name, ctype, code, len(df))
    df = df.copy()
    df = get_length(df)
    if len(df) == 1:
        dna_seq = df.iloc[0]["sequence"]
        construct_info["dna_sequence"] = dna_seq
        if not no_t7:  # need to remove t7 promoter
            if not dna_seq.startswith(t7_seq):
                log.error(f"{name} does not contain a t7 promoter")
                exit(1)
            dna_seq = dna_seq[20:]
        construct_info["rna_sequence"] = dna_seq.replace("T", "U")
    else:
        construct_info["fwd_p"] = "P001E"
        construct_info["rev_p"] = "P001F"
    construct_info["dna_len"] = df["length"].mean()
    if not no_t7:
        if not dna_seq.startswith(t7_seq):
            log.error(f"{name} does not contain a t7 promoter")
            exit(1)
        construct_info["rna_len"] = construct_info["dna_len"] - 20
        df = trim(df, 20, 0)
    else:
        construct_info["rna_len"] = construct_info["dna_len"]


def get_construct_entry(
    df: pd.DataFrame, name: str, code: str, no_t7: bool = False
) -> ConstructEntry:
    """
    gets the construct entry from dataframe
    :param df: the dataframe with sequences
    :param code: the code of the construct
    """

    def _assign_ctype(df: pd.DataFrame) -> str:
        """
        assigns the construct type
        """
        if len(df) == 1:
            return "ASSEMBLY"
        elif len(df) < 100:
            return "OPOOL"
        else:
            return "AGILENT"

    df = df.copy()
    df = get_length(df)
    centry = ConstructEntry(name, code, _assign_ctype(df))
    centry.size = len(df)
    if centry.ctype == "ASSEMBLY":
        centry.dna_len = df.iloc[0]["length"]
        centry.dna_sequence = df.iloc[0]["sequence"]
        if no_t7:
            centry.rna_sequence = df.iloc[0]["sequence"].replace("T", "U")
        else:
            centry.rna_sequence = df.iloc[0]["sequence"][20:].replace("T", "U")
        centry.rna_len = len(centry.rna_sequence)
        centry.rna_structure = vienna.folded_structure(centry.rna_sequence)
    else:
        centry.rna_len = df["length"].mean()
        centry.dna_len = df["length"].mean() + 20
        if no_t7:
            centry.dna_len -= 20
        centry.fwd_p = "P001E"
        centry.rev_p = "P001F"
    if not no_t7:
        df = trim(df, 20, 0)
    centry.seq_rev_p = get_seq_fwd_primer_code(df)
    if centry.seq_rev_p == "":
        log.warning(f"forward primer cannot be determined for {code}")
    return centry


def get_last_codes():
    """
    get the last codes from the google sheet
    :param params: module parameters
    """
    df_seqs = get_sequence_sheet()
    df_primers = get_oligo_sheet()
    last_code = df_seqs["code"].loc[df_seqs["code"].last_valid_index()]
    last_primer_code = df_primers["code"].loc[df_primers["code"].last_valid_index()]
    log.info("last construct code on sheet: %s", last_code)
    log.info("last primer code on sheet: %s", last_primer_code)
    return last_code, last_primer_code
