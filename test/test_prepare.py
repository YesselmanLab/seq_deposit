import pytest
import pandas as pd
import os
import multiprocessing

from seq_tools import fold
from seq_deposit.prepare import check_rt_seq, run_func_in_parallel

TEST_DIR = os.path.join(os.path.dirname(__file__))


def test_check_rt_seq():
    data = [
        {"sequence": "AAAAAAAAAAGAAACAACAACAACAAC", "name": "test"},
    ]
    df = pd.DataFrame(data)
    assert check_rt_seq(df) == "rt_tail"
    df.iloc[0]["sequence"] += "C"
    assert check_rt_seq(df) is None


def test_run_func_in_parallel():
    p = multiprocessing.Pool(2)
    path = os.path.join(TEST_DIR, "resources", "test.csv")
    df = pd.read_csv(path)
    df_result = run_func_in_parallel(fold, df, p)
    print(df_result)
