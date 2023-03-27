"""
handles processing of google sheets data
"""
import os
import pandas as pd
from dataclasses import dataclass
import wget

from seq_deposit.logger import get_logger
from seq_deposit.settings import LIB_PATH

import vienna
from seq_tools import to_dna, get_length, trim

log = get_logger("GSHEETS")

@dataclass(order=True)
class ConstructEntry:
    """
    construct entry
    """

    name: str
    code: str
    ctype: str
    arrived: str = "NO"
    usuable: str = "UNK"
    size: str = -1
    dna_len: str = -1
    dna_sequence: str = "LIBRARY"
    rna_len: str = -1
    rna_sequence: str = "LIBRARY"
    rna_structure: str = "LIBRARY"
    fwd_p: str = "NONE"
    rev_p: str = "NONE"
    rt_p: str = "RTB"
    seq_fwd_p: str = "P000R"
    seq_rev_p: str = "P000Y"
    project: str = ""
    comment: str = ""

    def to_csv_str(self):
        """
        converts to csv string
        """
        return ",".join(
            str(x)
            for x in [
                self.name,
                self.code,
                self.ctype,
                self.arrived,
                self.usuable,
                self.size,
                self.dna_len,
                self.dna_sequence,
                self.rna_len,
                self.rna_sequence,
                self.rna_structure,
                self.fwd_p,
                self.rev_p,
                self.rt_p,
                self.seq_fwd_p,
                self.seq_rev_p,
                self.project,
                self.comment,
            ]
        )

    def to_csv_header(self):
        """
        converts to csv string
        """
        return ",".join(
            str(x)
            for x in [
                "name",
                "code",
                "type",
                "arrived",
                "usable",
                "size",
                "dna_len",
                "dna_sequence",
                "rna_len",
                "rna_sequence",
                "rna_structure",
                "fwd_p",
                "rev_p",
                "rt_p",
                "seq_fwd_p",
                "seq_rev_p",
                "project",
                "comment",
            ]
        )


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

    log = get_logger("get_construct_entry")
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


def fetch_sequence_gsheet(params) -> pd.DataFrame:
    """
    fetches the sequence google sheet
    :param params: module params
    :return: dataframe of sequence sheet
    """
    gsheet_path = params["sequence_gsheet_url"]
    wget.download(gsheet_path, "temp.csv")
    df = pd.read_csv("temp.csv")
    os.remove("temp.csv")
    cols = (
        "name,code,type,arrived,usuable,size,dna_len,dna_sequence,rna_len,"
        "rna_sequence,rna_structure,fwd_p,rev_p,rt_p,seq_fwd_p,project,"
        "comment".split(",")
    )
    df = df[cols]
    path = os.path.join(LIB_PATH, "resources", "all.csv")
    df.to_csv(path, index=False)
    return df


def fetch_primers_gsheet(params) -> pd.DataFrame:
    """
    fetches the primer google sheet
    :param params: module params
    :return: dataframe of primer sheet
    """
    gsheet_path = params["primer_gsheet_url"]
    wget.download(gsheet_path, "temp.csv")
    cols = [
        "name",
        "code",
        "c_code",
        "useable",
        "a_temp",
        "sequence",
        "type",
        "",
    ]
    df = pd.read_csv("temp.csv", header=None, names=cols, skiprows=1)
    # remove last column that has nothing
    df = df.iloc[:, :-1]
    df.to_csv("temp.csv", index=False)
    os.remove("temp.csv")
    path = os.path.join(LIB_PATH, "resources", "primer.csv")
    df.to_csv(path, index=False)
    return df


def get_last_codes(params):
    """
    get the last codes from the google sheet
    :param params: module parameters
    """
    log = get_logger("get_last_codes")
    df_seqs = fetch_sequence_gsheet(params)
    df_primers = fetch_primers_gsheet(params)
    last_code = df_seqs.iloc[-1]["code"]
    last_primer_code = df_primers.iloc[-1]["code"]
    log.info("last construct code on sheet: %s", last_code)
    log.info("last primer code on sheet: %s", last_primer_code)
    return last_code, last_primer_code


# TODO finish implementing
def get_p5_valid_sequences(df: pd.DataFrame) -> pd.DataFrame:
    """
    gets the p5 sequence
    :param df: the dataframe with sequences
    """
    df = df[df["type"] == "SEQ_REV"]
    return df
