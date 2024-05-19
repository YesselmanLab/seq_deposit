import os
import pandas as pd
import multiprocessing
from typing import List

from seq_tools import (
    get_length,
    get_molecular_weight,
    get_default_names,
    get_extinction_coeff,
    to_dna,
    to_fasta,
    trim,
    transcribe,
    has_5p_sequence,
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


def split_dataframe(df, chunk_size=10000) -> List[pd.DataFrame]:
    """
    Splits a dataframe into smaller chunks.

    Args:
        df (pandas.DataFrame): The dataframe to be split.
        chunk_size (int, optional): The size of each chunk. Defaults to 10000.

    Returns:
        list: A list of dataframes, each representing a chunk of the original dataframe.
    """
    chunks = list()
    num_chunks = len(df) // chunk_size + 1
    for i in range(num_chunks):
        chunks.append(df[i * chunk_size : (i + 1) * chunk_size])
    return chunks


def run_func_in_parallel(
    func, df: pd.DataFrame, p: multiprocessing.Pool
) -> pd.DataFrame:
    """
    Runs a given function in parallel on a DataFrame using a multiprocessing Pool.

    Args:
        func: The function to be executed in parallel.
        df (pd.DataFrame): The DataFrame to be processed.
        p (multiprocessing.Pool): The multiprocessing Pool object.

    Returns:
        pd.DataFrame: The concatenated DataFrame after executing the function in parallel.
    """
    inputs = split_dataframe(df, int(len(df) / len(p._pool)) + 1)
    dfs = p.starmap(func, zip(inputs))
    df = pd.concat(dfs)
    return df


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
    if ntype == "DNA":
        df_dna = to_dna(df_dna)
        if not ignore_missing_t7:
            if not has_5p_sequence(df_dna, t7_seq):
                log.error("Missing T7 promoter sequence")
                exit()
    elif ntype == "RNA":
        pass

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
