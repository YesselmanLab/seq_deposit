import io
from typing import Tuple

import pandas as pd
from tabulate import tabulate

from seq_tools import to_fasta, to_dna

from gsheets.sheet import get_next_code, get_sequence_sheet, get_oligo_sheet

from seq_deposit.logger import get_logger
from seq_deposit.construct_info import get_construct_entry


log = get_logger("process")


def process_primers(df: pd.DataFrame) -> pd.DataFrame:
    """
    Process primers in the given DataFrame and generate output files.

    Args:
        df (pd.DataFrame): Input DataFrame containing primer information.

    Returns:
        pd.DataFrame: Processed DataFrame containing primer information.
    """
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
    buffer = io.StringIO()
    df_assembly.to_csv(buffer, index=False)
    log.info("ASSEMBLIES----------------------------------:\n%s", buffer.getvalue())
    df_primers.drop(columns=["p_num", "name"], inplace=True)
    df_primers.rename(columns={"p_name": "name", "code": "c_code"}, inplace=True)
    df_primers.rename(columns={"p_code": "code"}, inplace=True)
    df_primers = df_primers[
        ["name", "code", "c_code", "useable", "a_temp", "sequence", "type"]
    ]
    df_primers.drop_duplicates(subset="name", inplace=True)
    df_primers.dropna(inplace=True)
    df_primers.to_csv("seq-deposit-output/primers.csv", index=False)
    buffer = io.StringIO()
    df_primers.to_csv(buffer, index=False)
    log.info("PRIMERS----------------------------------:\n%s", buffer.getvalue())
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
    return df_primers


def process_constructs(df: pd.DataFrame) -> pd.DataFrame:
    """
    Process constructs from a DataFrame.

    Args:
        df (pd.DataFrame): The input DataFrame containing construct information.

    Returns:
        pd.DataFrame: The processed DataFrame containing construct entries.

    Raises:
        None
    """
    required_cols = [
        "name",
        "sequence",
        "construct",
        "type",
        "dna_sequence",
        "structure",
        "mfe",
        "ens_defect",
    ]
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"missing required column: {col}")
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
    buffer = io.StringIO()
    df_constructs.to_csv(buffer, index=False)
    log.info("CONSTRUCTS----------------------------------:\n%s", buffer.getvalue())
    df_sub = df.query("type == 'ASSEMBLY'")
    if len(df_sub) > 0:
        process_primers(df_sub)
    df_constructs.to_csv("seq-deposit-output/constructs.csv", index=False)
    return df_constructs


def get_last_codes() -> Tuple[str, str]:
    """
    Get the last codes from the Google Sheet.

    Returns:
        A tuple containing the last construct code and the last primer code.
    """
    df_seqs = get_sequence_sheet()
    df_primers = get_oligo_sheet()
    last_code = df_seqs["code"].loc[df_seqs["code"].last_valid_index()]
    last_primer_code = df_primers["code"].loc[df_primers["code"].last_valid_index()]
    log.info("last construct code on sheet: %s", last_code)
    log.info("last primer code on sheet: %s", last_primer_code)
    return last_code, last_primer_code
