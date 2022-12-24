"""
test module for gsheets module
"""
import os

import pandas as pd

from seq_tools import to_dna, get_length

from seq_deposit.gsheets import (
    fetch_sequence_gsheet,
    fetch_primers_gsheet,
    get_p5_valid_sequences,
    get_construct_entry,
)

TEST_RESOURCES_PATH = os.path.join(os.path.dirname(__file__), "resources")


# generate test data ##########################################################


def get_test_data_rna() -> pd.DataFrame:
    """
    get test data for dna
    :return: pd.DataFrame
    """
    return pd.DataFrame(
        [
            ["seq_0", "GGGGUUUUCCCC", "((((....))))"],
        ],
        columns=["name", "sequence", "structure"],
    )


# tests  ######################################################################


def test_fetch_sequence_gsheet():
    """
    tests the fetch_sequence_gsheet function
    """
    params = {
        "sequence_gsheet_url": "https://docs.google.com/spreadsheets/d/1ihU99Dt"
        "K10J9yPDLiSQ9KdEC7KaddIn7Ta3tb_HFz8U/gviz/tq?tq"
        "x=out:csv&sheet=all"
    }
    df = fetch_sequence_gsheet(params)
    assert "name" in df.columns


def test_fetch_primers_gsheet():
    """
    tests the fetch_primers_gsheet function
    """
    params = {
        "primer_gsheet_url": "https://docs.google.com/spreadsheets/d/1_oGLJRZMG"
        "8zNoc6yF5lJV7k0MVP6XgVD9A2ya0XPAXQ/gviz/tq?"
        "tqx=out:csv&sheet=oligos"
    }
    df = fetch_primers_gsheet(params)
    assert "name" in df.columns


def test_get_p5_valid_sequences():
    """
    tests the get_p5_valid_sequences function
    """
    params = {
        "primer_gsheet_url": "https://docs.google.com/spreadsheets/d/1_oGLJRZMG"
        "8zNoc6yF5lJV7k0MVP6XgVD9A2ya0XPAXQ/gviz/tq?"
        "tqx=out:csv&sheet=oligos"
    }
    df = fetch_primers_gsheet(params)
    df = get_p5_valid_sequences(df)


def test_get_construct_entry_opool():
    """
    tests the get_construct_entry function
    """
    df = pd.read_csv(os.path.join(TEST_RESOURCES_PATH, "test.csv"))
    df = to_dna(df)
    df = get_length(df)
    centry = get_construct_entry(df, "test_lib", "C0001", no_t7=True)
    assert centry.code == "C0001"
    assert centry.ctype == "OPOOL"
    assert centry.arrived == "NO"
    assert centry.name == "test_lib"
    assert centry.usuable == "UNK"
    assert centry.size == len(df)
    assert centry.dna_len == df["length"].mean()
    assert centry.dna_sequence == "LIBRARY"
    assert centry.rna_len == df["length"].mean()
    assert centry.rna_sequence == "LIBRARY"
    assert centry.rna_structure == "LIBRARY"
    assert centry.rev_p == "P001F"
    assert centry.fwd_p == "P001E"


def test_get_construct_entry_assembly():
    """
    tests the get_construct_entry function
    """
    df = get_test_data_rna()
    df = get_length(df)
    df = to_dna(df)
    centry = get_construct_entry(df, "test_lib", "C0001", no_t7=True)
    assert centry.code == "C0001"
    assert centry.ctype == "ASSEMBLY"
    assert centry.dna_len == len(df.iloc[0]["sequence"])
    assert centry.dna_sequence == df.iloc[0]["sequence"]
    assert centry.rna_len == len(df.iloc[0]["sequence"])
    assert centry.rna_sequence == df.iloc[0]["sequence"].replace("T", "U")
    assert centry.rna_structure == "((((....))))"
    assert centry.fwd_p == "NONE"
    assert centry.rev_p == "NONE"
    assert centry.seq_fwd_p == "P000R"
    assert centry.seq_rev_p == ""


def test_get_construct_entry_assembly_t7():
    """
    tests the get_construct_entry function
    """
    df = get_test_data_rna()
    df = to_dna(df)
    df.iloc[0]["sequence"] = "TTCTAATACGACTCACTATA" + df.iloc[0]["sequence"]
    df = get_length(df)
    centry = get_construct_entry(df, "test_lib", "C0001", no_t7=False)
    assert centry.dna_len == len(df.iloc[0]["sequence"])
    assert centry.dna_sequence == df.iloc[0]["sequence"]
    assert centry.rna_len == len(df.iloc[0]["sequence"]) - 20
    assert centry.rna_sequence == df.iloc[0]["sequence"][20:].replace("T", "U")
    assert centry.seq_rev_p == ""


def test_get_construct_entry_assembly_p5():
    """
    tests the get_construct_entry function
    """
    df = get_test_data_rna()
    df = to_dna(df)
    df.loc[0]["sequence"] = (
        "TTCTAATACGACTCACTATA" + "GGCCAAAACAACGG" + df.iloc[0]["sequence"]
    )
    centry = get_construct_entry(df, "test_lib", "C0001")
    assert centry.seq_rev_p == "P0016"
    df.loc[1] = ["seq_1", "GGGGTTTTCCCC"]
    df.loc[1]["sequence"] = (
        "TTCTAATACGACTCACTATA" + "GGCCAATACAACGG" + df.iloc[0]["sequence"]
    )
    centry = get_construct_entry(df, "test_lib", "C0001")
    assert centry.seq_rev_p == ""
