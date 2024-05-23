"""
test module for gsheets module
"""

import os

import pandas as pd

from seq_tools import to_dna, get_length

from seq_deposit.gsheets import get_last_codes

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


def test_get_last_codes():
    last_construct_code, last_primer_code = get_last_codes()
    assert last_construct_code.startswith("C")
    assert last_primer_code.startswith("P")
