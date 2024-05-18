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
from seq_deposit.settings import LIB_PATH

log = get_logger(__file__)


def check_rt_seq(df: pd.DataFrame) -> str:
    """
    Check if all sequences in the given DataFrame end with a specific reverse
    transcription (RT) primer sequence.

    Args:
        df (pd.DataFrame): The DataFrame containing the sequences to check.

    Returns:
        str: The name of the reverse transcription sequence if all sequences in the DataFrame
            end with it, otherwise returns None.
    """
    path = os.path.join(LIB_PATH, "resources", "rt_seqs.csv")
    df_p5 = pd.read_csv(path)
    for _, row in df_p5.iterrows():
        # if all sequences in df start with the p5 sequence then return the p5 code
        if all(df["sequence"].str.endswith(row["sequence"])):  # type: ignore
            return row["name"]
    return None


def generate_dna_dataframe(
    df: pd.DataFrame,
    ntype: str,
    ignore_missing_t7: bool = False,
    ignore_missing_rt_seq: bool = False,
    t7_seq: str = "TTCTAATACGACTCACTATA",
) -> pd.DataFrame:
    df_dna = df.copy()
    if "name" not in df_dna.columns:
        df_dna = get_default_names(df_dna)
    df_dna = to_dna(df_dna)
    df_dna = df_dna[["name", "sequence"]]
    df_dna = get_length(df_dna)
    df_dna = get_molecular_weight(df_dna, "DNA", True)
    df_dna = get_extinction_coeff(df_dna, "DNA", True)
    return df_dna


def generate_rna_dataframe(df: pd.DataFrame, ignore_missing_t7: bool) -> pd.DataFrame:
    """
    generates the rna dataframe
    :param df: the dataframe with sequences
    :param ignore_missing_t7: if true, ignore missing t7 promoter
    :return: the rna dataframe
    """
    df_rna = df.copy()
    df_rna = transcribe(df_rna, ignore_missing_t7)
    df_rna = df_rna[["name", "sequence", "structure"]]
    df_rna = get_length(df_rna)
    df_rna = get_molecular_weight(df_rna, "RNA", False)
    df_rna = get_extinction_coeff(df_rna, "RNA", False)
    return df_rna
