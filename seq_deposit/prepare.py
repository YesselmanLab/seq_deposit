import os
import pandas as pd
import multiprocessing
from typing import List

from seq_tools import (
    get_length,
    get_default_names,
    fold,
    to_dna,
    to_rna,
    to_fasta,
    trim,
    transcribe,
    has_5p_sequence,
    add,
)

from seq_deposit.logger import get_logger
from seq_deposit.settings import LIB_PATH

log = get_logger(__name__)


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
    """
    Generate a DNA dataframe based on the input dataframe.

    Args:
        df (pd.DataFrame): The input dataframe.
        ntype (str): The type of nucleotide sequence ('DNA' or 'RNA').
        ignore_missing_t7 (bool, optional): Whether to ignore missing T7 promoter sequence. Defaults to False.
        ignore_missing_rt_seq (bool, optional): Whether to ignore missing RT sequence. Defaults to False.
        t7_seq (str, optional): The T7 promoter sequence. Defaults to "TTCTAATACGACTCACTATA".

    Returns:
        pd.DataFrame: The generated DNA dataframe.
    """
    df_dna = df.copy()
    if "name" not in df_dna.columns:
        df_dna = get_default_names(df_dna)
    if ntype == "DNA":
        if not ignore_missing_t7:
            if not has_5p_sequence(df_dna, t7_seq):
                log.error("Missing T7 promoter sequence")
                return None
    elif ntype == "RNA":
        if not ignore_missing_rt_seq:
            rt_seq = check_rt_seq(df_dna)
            if rt_seq is None:
                log.error("Missing RT sequence")
                return None
        df_dna = add(df_dna, t7_seq, "")
    df_dna = to_dna(df_dna)
    df_dna = df_dna[["name", "sequence"]]
    df_dna = get_length(df_dna)
    return df_dna


def generate_rna_dataframe(
    df: pd.DataFrame,
    ntype: str,
    ignore_missing_t7: bool = False,
    ignore_missing_rt_seq: bool = False,
    t7_seq: str = "TTCTAATACGACTCACTATA",
    threads: int = 1,
) -> pd.DataFrame:
    """
    Generate an RNA dataframe from the given DNA dataframe.

    Args:
        df (pd.DataFrame): The input DNA dataframe.
        ntype (str): The type of nucleotide sequence (e.g., "DNA").
        ignore_missing_t7 (bool, optional): Whether to ignore missing T7 promoter sequence. Defaults to False.
        ignore_missing_rt_seq (bool, optional): Whether to ignore missing reverse transcription sequence. Defaults to False.
        t7_seq (str, optional): The T7 promoter sequence. Defaults to "TTCTAATACGACTCACTATA".
        threads (int, optional): The number of threads to use for parallel processing. Defaults to 1.

    Returns:
        pd.DataFrame: The RNA dataframe with columns "name", "sequence", "structure", and "length".
    """
    df_rna = df.copy()
    if "name" not in df_rna.columns:
        df_rna = get_default_names(df_rna)
    if ntype == "DNA":
        if not ignore_missing_t7:
            if not has_5p_sequence(df_rna, t7_seq):
                log.error("Missing T7 promoter sequence")
                return None
            df_rna = trim(df_rna, len(t7_seq), 0)
    df_rna = to_rna(df_rna)
    if "structure" not in df_rna.columns or "ens_defect" not in df_rna.columns:
        if threads > 1:
            p = multiprocessing.Pool(threads)
            df_rna = run_func_in_parallel(fold, df_rna, p)
        else:
            df_rna = fold(df_rna)
    df_rna = df_rna[["name", "sequence", "structure", "ens_defect"]]
    df_rna = get_length(df_rna)
    return df_rna
