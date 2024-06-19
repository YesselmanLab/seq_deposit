"""
command line interface for seq_deposit script
"""

import os
import shutil
import cloup
import pandas as pd
import sys
import glob
from pathlib import Path
from typing import List
from tabulate import tabulate

from seq_tools import to_fasta, to_dna

from gsheets.sheet import get_next_code, get_sequence_sheet, get_oligo_sheet

from seq_deposit.logger import setup_logging, get_logger
from seq_deposit.deposit import deposit_files, deposit_rna_csv
from seq_deposit.prepare import generate_dna_dataframe, generate_rna_dataframe
from seq_deposit.construct_info import get_construct_entry


log = get_logger("CLI")


# helper functions ############################################################


def get_rna_dataframe_from_row(row, columns):
    row_df = pd.DataFrame([row], columns=columns, index=[0])
    row_df = row_df[["name", "rna_sequence", "rna_structure"]]
    row_df.rename(
        {"rna_sequence": "sequence", "rna_structure": "structure"},
        axis=1,
        inplace=True,
    )
    row_df = generate_rna_dataframe(row_df, "RNA")
    return row_df


def get_dna_dataframe_from_row(row, columns):
    row_df = pd.DataFrame([row], columns=columns, index=[0])
    row_df = row_df[["name", "dna_sequence"]]
    row_df.rename(
        {"dna_sequence": "sequence"},
        axis=1,
        inplace=True,
    )
    row_df = generate_dna_dataframe(row_df, "DNA")
    return row_df


def fasta_to_dataframe(fasta_file):
    """
    Reads a FASTA file and converts it into a DataFrame with columns 'name' and 'sequence'.

    Parameters:
    fasta_file (str): Path to the FASTA file.

    Returns:
    pd.DataFrame: DataFrame with columns 'name' and 'sequence'.
    """
    names = []
    sequences = []
    sequence = []

    with open(fasta_file, "r") as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if sequence:
                    sequences.append("".join(sequence))
                    sequence = []
                names.append(line[1:])
            else:
                sequence.append(line)
        if sequence:
            sequences.append("".join(sequence))

    data = {"name": names, "sequence": sequences}
    return pd.DataFrame(data)


def validate_rna_csv(row, seq_path):
    if not os.path.exists(f"{seq_path}/rna/{row['code']}.csv"):
        log.warning(f"rna csv does not exist for {row['code']}")
        return False
    df_rna = pd.read_csv(f"{seq_path}/rna/{row['code']}.csv")
    if row["type"] == "ASSEMBLY":
        csv_row = df_rna.iloc[0]
        if row["rna_sequence"] != csv_row["sequence"]:
            log.warning(f"rna sequence does not match for {row['code']}")
            return False
        return True


def validate_dna_csv(row, seq_path):
    if not os.path.exists(f"{seq_path}/dna/{row['code']}.csv"):
        log.warning(f"dna csv does not exist for {row['code']}")
        return False
    df_rna = pd.read_csv(f"{seq_path}/dna/{row['code']}.csv")
    if row["type"] == "ASSEMBLY":
        csv_row = df_rna.iloc[0]
        if row["dna_sequence"] != csv_row["sequence"]:
            log.warning(f"dna sequence does not match for {row['code']}")
            return False
        return True


def setup(
    ignore_missing_t7: bool,
    ignore_missing_rt_seq: bool,
    t7_seq: str,
    rt_seq: str,
    overwrite: bool,
) -> None:
    """
    Set up the necessary configurations and directories for seq-deposit.

    Args:
        ignore_missing_t7 (bool): Flag indicating whether to ignore missing T7 promoter sequence.
        ignore_missing_rt_seq (bool): Flag indicating whether to ignore missing RT sequence.
        t7_seq (str): The T7 promoter sequence.
        overwrite (bool): Flag indicating whether to overwrite existing output directory.

    Returns:
        None
    """
    # setup logger
    setup_logging()
    log.info("setting up seq-deposit")
    log.info("ignore missing t7: %s", ignore_missing_t7)
    log.info("ignore missing rt seq: %s", ignore_missing_rt_seq)
    log.info("t7 seq: %s", t7_seq)
    log.info("rt seq: %s", rt_seq)
    # create output directory
    if os.path.exists("seq-deposit-output") and not overwrite:
        log.error("seq-deposit-output already exists use --overwrite")
        exit(1)
    log.info("outputs will be written in seq-deposit-output")
    if not os.path.exists("seq-deposit-output"):
        os.makedirs("seq-deposit-output")
    log.info("ran at commandline as: ")
    log.info(" ".join(sys.argv))
    os.makedirs("seq-deposit-output/dna", exist_ok=True)
    os.makedirs("seq-deposit-output/rna", exist_ok=True)
    os.makedirs("seq-deposit-output/fastas", exist_ok=True)
    if ignore_missing_t7:
        log.warning("ignore missing t7 promoter. This is not recommended!!")
    else:
        log.info("t7 promoter sequence: %s", t7_seq)
    if ignore_missing_rt_seq:
        log.warning("ignore missing rt seq. This is not recommended!!")


def get_last_codes():
    """
    get the last codes from the google sheet
    :param params: module parameters
    """
    log = get_logger("get_last_codes")
    df_seqs = get_sequence_sheet()
    df_primers = get_oligo_sheet()
    last_code = df_seqs["code"].loc[df_seqs["code"].last_valid_index()]
    last_primer_code = df_primers["code"].loc[df_primers["code"].last_valid_index()]
    log.info("last construct code on sheet: %s", last_code)
    log.info("last primer code on sheet: %s", last_primer_code)
    return last_code, last_primer_code


def setup_output_dir(df):
    pass


# logging #####################################################################
# these functions are to create informationa about the new sequences


def log_new_libaries(params, new_libs, codes):
    """
    logs the new libraries
    :param params: module params
    :param new_libs: the new libraries
    :return: None
    """
    path = os.path.join(params["deposit_path"], "org_data_csvs", "log.csv")
    df_log = pd.read_csv(path)
    # create backup in case something goes wrong
    df_log.to_csv(
        os.path.join(params["deposit_path"], "org_data_csvs", "log_old.csv"),
        index=False,
    )
    pos = len(df_log) - 1
    for lib, code in zip(new_libs, codes):
        df_log.loc[pos] = [
            Path(lib).stem,
            code,
            os.path.abspath(lib),
        ]
        shutil.copy(
            lib,
            os.path.join(
                params["deposit_path"], "org_data_csvs", "csvs", f"{code}.csv"
            ),
        )
        pos += 1
    df_log.to_csv(path, index=False)


def log_constructs(constructs):
    """
    logs the new constructs
    :param constructs: the new constructs
    :return: None
    """
    log.info("Add seq-deposit-output/constructs.csv to sequences.gsheet on all sheet")
    with open("seq-deposit-output/constructs.csv", "w") as f:
        f.write(constructs[0].to_csv_header() + "\n")
        for c in constructs:
            f.write(c.to_csv_str() + "\n")


def log_primers(primers):
    """
    logs the new primers
    :param primers: the new primers
    :return: None
    """
    log.info("seq-deposit-output/primers.csv to the primers.gsheet on oligo sheet")
    with open("seq-deposit-output/primers.csv", "w") as f:
        f.write("name,code,c_code,useable,a_temp,sequence,type\n")
        for p in primers:
            f.write(",".join([str(x) for x in p]) + "\n")


# cli functions ###############################################################
def common_options():
    """
    Returns a Cloup option group with common options.

    Returns:
        cloup.OptionGroup: A Cloup option group containing the common options.

    Options:
        --ignore-missing-t7: Ignore the t7 promoter. (default: False)
        --ignore-missing-rt-seq: Ignore the rt sequence. (default: False)
        --overwrite: Overwrite existing files. (default: False)
        --t7-seq: The t7 promoter sequence. (default: "TTCTAATACGACTCACTATA")
        --rt-seq: The rt sequence. (default: "AAAGAAACAACAACAACAAC")
    """
    return cloup.option_group(
        "common options",
        cloup.option("--ignore-missing-t7", is_flag=True, help="ignore t7 promoter"),
        cloup.option("--ignore-missing-rt-seq", is_flag=True, help="ignore rt seq"),
        cloup.option("--overwrite", is_flag=True, help="overwrite existing files"),
        cloup.option(
            "--t7-seq", default="TTCTAATACGACTCACTATA", help="t7 promoter sequence"
        ),
        cloup.option("--rt-seq", default="AAAGAAACAACAACAACAAC", help="rt sequence"),
    )


@cloup.group()
def cli():
    """
    command line interface for seq_deposit script
    """
    pass


@cli.command(help="")
@common_options()
def from_final_dir(**kwargs):
    setup(**kwargs)
    if os.path.isdir("final"):
        log.info("final directory exists")
    else:
        log.error("final directory does not exist")
        exit()
    method = None
    try:
        with open("final/METHOD") as f:
            method = f.read().rstrip()
    except:
        pass
    last_code, last_primer_code = get_last_codes()
    df_seqs = get_sequence_sheet()
    df = pd.read_json("final/summary_jsons/without_codes.json")
    df["primer_codes"] = [[] for _ in range(len(df))]
    constructs = []
    primers = []
    for construct, g in df.groupby("construct"):
        df_dna = g[["name", "dna_sequence"]].copy()
        df_dna.rename(columns={"dna_sequence": "sequence"}, inplace=True)
        df_rna = g[["name", "sequence", "structure", "mfe", "ens_defect"]].copy()
        # standardize
        df_dna = generate_dna_dataframe(
            df_dna,
            "DNA",
            kwargs["ignore_missing_t7"],
            kwargs["ignore_missing_rt_seq"],
            kwargs["t7_seq"],
        )
        df_rna = generate_rna_dataframe(
            df_rna,
            "RNA",
            kwargs["ignore_missing_t7"],
            kwargs["ignore_missing_rt_seq"],
            kwargs["t7_seq"],
        )
        log.info("processing %s it has %d sequence(s)", construct, len(df_dna))
        code = get_next_code(last_code)
        construct_entry = get_construct_entry(df_dna, df_rna, construct, method, code)
        constructs.append(construct_entry)
        df_dna.to_csv(f"seq-deposit-output/dna/{code}.csv", index=False)
        df_rna.to_csv(f"seq-deposit-output/rna/{code}.csv", index=False)
        df_fasta = df_rna[["name", "sequence"]]
        df_fasta = to_dna(df_fasta)
        to_fasta(df_fasta, f"seq-deposit-output/fastas/{code}.fasta")
        last_code = code
        if not g.iloc[0]["has_primers"]:
            continue
        for i, row in g.iterrows():
            pcodes = []
            pos = 1
            for pname, pseq in zip(row["primer_names"], row["primers"]):
                pcode = get_next_code(last_primer_code)
                primers.append(
                    [pname, pcode, construct, code, "UNK", "N/A", pseq, "ASSEMBLY", pos]
                )
                pcodes.append(pcode)
                last_primer_code = pcode
                pos += 1
            df.at[i, "primer_codes"] = pcodes
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
    df_assembly = df_primers[["construct", "code", "p_num", "p_name", "p_code"]]
    df_assembly.to_csv("seq-deposit-output/assemblies.csv", index=False)
    df_primers.drop(columns=["p_num", "construct"], inplace=True)
    df_primers.rename(columns={"p_name": "name", "code": "c_code"}, inplace=True)
    df_primers.rename(
        columns={
            "p_code": "code",
        },
        inplace=True,
    )
    df_primers = df_primers[
        ["name", "code", "c_code", "useable", "a_temp", "sequence", "type"]
    ]
    df_primers.to_csv("seq-deposit-output/primers.csv", index=False)
    df_constructs = pd.DataFrame(constructs, columns=df_seqs.columns)
    log.info(
        "\n"
        + tabulate(
            df_constructs[["name", "code", "type", "size", "dna_len"]],
            headers="keys",
            tablefmt="psql",
        )
    )
    df_constructs.to_csv("seq-deposit-output/constructs.csv", index=False)


@cli.command(help="generate required files from final directory")
@cloup.argument("final_dir")
def finalize(final_dir):
    """
    finalize constructs that will be ordered as agilent libraries
    """


@cli.command(help="generate required files for opools")
@cloup.option("--ignore-missing-t7", is_flag=True, help="ignore t7 promoter")
@cloup.option("--ignore-missing-rt-seq", is_flag=True, help="ignore rt seq")
@cloup.option("--overwrite", is_flag=True, help="overwrite existing files")
@cloup.option("--ntype", default="DNA", help="type of nucleic acid")
@cloup.option("--t7-seq", default="TTCTAATACGACTCACTATA", help="t7 promoter sequence")
@cloup.option("--threads", default=1, help="number of threads to use")
@cloup.argument("csvs", nargs=-1)
def opools(
    csvs: List[str],
    ntype: str,
    ignore_missing_t7: bool,
    ignore_missing_rt_seq: bool,
    t7_seq: str,
    threads: int,
    overwrite: bool,
) -> None:
    """
    Generate required files for opools.

    Args:
        csvs (List[str]): List of CSV file paths.
        ntype (str): Type of nucleic acid. Default is "DNA".
        ignore_missing_t7 (bool): Flag to ignore t7 promoter. Default is False.
        ignore_missing_rt_seq (bool): Flag to ignore rt seq. Default is False.
        t7_seq (str): t7 promoter sequence. Default is "TTCTAATACGACTCACTATA".
        threads (int): Number of threads to use. Default is 1.
        overwrite (bool): Flag to overwrite existing files. Default is False.

    Returns:
        None
    """
    setup(ignore_missing_t7, ignore_missing_rt_seq, t7_seq, overwrite)
    log.info("processing %d csvs", len(csvs))
    log.info("ntype: %s", ntype)
    last_code, _ = get_last_codes()
    constructs = []
    df_seqs = get_sequence_sheet()
    for csv in csvs:
        if ntype == "DNA":
            df_dna = pd.read_csv(csv)
            df_rna = generate_rna_dataframe(
                df_dna, ntype, ignore_missing_t7, ignore_missing_rt_seq, t7_seq, threads
            )
        else:
            df_rna = pd.read_csv(csv)
            df_dna = generate_dna_dataframe(
                df_rna, ntype, ignore_missing_t7, ignore_missing_rt_seq, t7_seq
            )
        ctype = "OPOOL"
        if len(df_dna) > 100:
            ctype = "AGILENT"
        log.info("processing %s it has %d sequences", csv, len(df_dna))
        code = get_next_code(last_code)
        fname = Path(csv).stem
        construct_entry = get_construct_entry(df_dna, df_rna, fname, ctype, code)
        constructs.append(construct_entry)
        df_dna.to_csv(f"seq-deposit-output/dna/{code}.csv", index=False)
        df_rna.to_csv(f"seq-deposit-output/rna/{code}.csv", index=False)
        df_fasta = df_rna[["name", "sequence"]]
        df_fasta = to_dna(df_fasta)
        to_fasta(df_fasta, f"seq-deposit-output/fastas/{code}.fasta")
        last_code = code
    df_constructs = pd.DataFrame(constructs, columns=df_seqs.columns)
    df_constructs.to_csv("seq-deposit-output/constructs.csv", index=False)


@cli.command(help="generate required files for primer assemblies")
@cloup.argument("construct_csv")
@cloup.argument("primer_csv")
@cloup.option("--ignore-missing-t7", is_flag=True, help="ignore t7 promoter")
@cloup.option("--ignore-missing-rt-seq", is_flag=True, help="ignore rt seq")
@cloup.option("--t7-seq", default="TTCTAATACGACTCACTATA", help="t7 promoter sequence")
@cloup.option("--overwrite", is_flag=True, help="overwrite existing files")
def assembly(
    construct_csv,
    primer_csv,
    ignore_missing_t7,
    ignore_missing_rt_seq,
    t7_seq,
    overwrite,
):
    setup(ignore_missing_t7, ignore_missing_rt_seq, t7_seq, overwrite)
    last_code, last_primer_code = get_last_codes()
    df_constructs = pd.read_csv(construct_csv)
    df_primers = pd.read_csv(primer_csv)
    log.info("processing %s it has %d sequences", construct_csv, len(df_constructs))
    constructs = []
    p_codes = {}
    codes = {}
    all_primers = []
    for _, row in df_constructs.iterrows():
        df_row = pd.DataFrame([row], columns=df_constructs.columns, index=[0])
        df_primers_sub = df_primers[df_primers["name"] == row["name"]]
        code = get_next_code(last_code)
        codes[row["name"]] = code
        for _, p_row in df_primers_sub.iterrows():
            if p_row["primer_name"] not in p_codes:
                p_codes[p_row["primer_name"]] = get_next_code(last_primer_code)
                last_primer_code = p_codes[p_row["primer_name"]]
                all_primers.append(
                    [
                        p_row["primer_name"],
                        last_primer_code,
                        code,
                        "UNK",
                        "N/A",
                        p_row["primer_sequence"],
                        "ASSEMBLY",
                    ]
                )
        deposit_files(df_row, code, params["deposit_path"], dry_run, ignore_missing_t7)
        centry = get_construct_entry(df_row, row["name"], code, ignore_missing_t7)
        last_code = code
        constructs.append(centry)
    df_primers["primer_code"] = df_primers["primer_name"].apply(lambda x: p_codes[x])
    df_primers["code"] = df_primers["name"].apply(lambda x: codes[x])
    df_primers = df_primers[
        ["name", "code", "primer_num", "primer_name", "primer_code"]
    ]
    df_primers.to_csv("seq-deposit-output/assemblies.csv", index=False)


# TODO refactor this should supply either an csv and code
@cli.command(help="update libraries")
@cloup.argument("construct_csv")
@cloup.argument("csvs", nargs=-1)
def update_libraries(construct_csv, csvs):
    params = get_params()
    df_construct = pd.read_csv(construct_csv)
    for i, row in df_construct.iterrows():
        name = row["name"]
        current_csv = None
        for csv in csvs:
            if name in csv:
                current_csv = csv
                break
        df = pd.read_csv(current_csv)
        log.info("processing %s it has %d sequences", csv, len(df))
        if not has_t7_promoter(df):
            raise ValueError(f"no t7 promoter found in sequences in {csv}")
        code = row["code"]
        deposit_files(df, code, params["deposit_path"], False, False)


@cli.command(help="get last code")
def get_last_code():
    params = get_params()
    df_seqs = fetch_sequence_gsheet(params)
    df_primers = fetch_primers_gsheet(params)
    print(df_seqs.iloc[-1])
    print(df_primers.iloc[-1])


# TODO search by name size etc return csv with results
@cli.command()
@cloup.option("--code", default=None, help="code of the sequence")
def get_sequence_info(code):
    setup_logging()
    df = get_sequence_sheet()
    if code is not None:
        df = df[df["code"] == code]
        if len(df) == 0:
            log.error(f"code {code} not found")
            return
    log.info("num of constructs %d", len(df))
    seq_path = os.environ["SEQPATH"]
    for i, row in df.iterrows():
        if row["size"] == 1:
            print("name: ", row["name"])
            print("sequence: ", row["rna_sequence"])
            print("structure: ", row["rna_structure"])
        else:
            df_rna = pd.read_csv(f"{seq_path}/rna/{row['code']}.csv")
            rna_row = df_rna.iloc[0]
            print("example sequence ->")
            print("name: ", rna_row["name"])
            print("rna sequence: ", rna_row["sequence"])
            print("rna structure: ", rna_row["structure"])
        print("code: ", row["code"])
        print("type: ", row["type"])
        print("-" * 40)


# TODO takes a seq-desposit directory and actually transfers files!
@cli.command()
@cloup.option("--deposit-path", default=None, help="deposit path")
@cloup.option("--force", is_flag=True, help="force copy of files")
def deposit(deposit_path, force):
    # check if the deposit path exists if not use supplied option
    # deposit_files(df, code, params["deposit_path"], dry_run, ignore_missing_t7)
    setup_logging()
    output = input(
        "are you sure you want to deposit all files (i.e. did you check they are corect)? [y/n]"
    )
    if output == "y":
        print("depositing files")
    else:
        return
    if not os.path.exists("seq-deposit-output"):
        log.error("seq-deposit-output does not exist")
        return
    if deposit_path is None:
        pass


@cli.command()
@cloup.option("--debug", is_flag=True, help="debug mode")
def check_sequence_sheet(debug):
    setup_logging()
    if not debug:
        log.info("backing up sequence and oligo directory")
        os.system("zip -r sequences_and_oligos.zip $SEQPATH")
    df = get_sequence_sheet()
    seq_path = os.environ["SEQPATH"]
    for i, row in df.iterrows():
        if row["type"] == "ASSEMBLY" and row["usuable"] != "NO":
            if not validate_rna_csv(row, seq_path):
                row_df = get_rna_dataframe_from_row(row, df.columns)
                deposit_rna_csv(row_df, row["code"], seq_path)
            if not validate_dna_csv(row, seq_path):
                print("not a valid DNA sequence")
            # df_rna = pd.read_csv(f"{seq_path}/rna/{row['code']}.csv")
            # df_dna = pd.read_csv(f"{seq_path}/dna/{row['code']}.csv")
            # df_fasta = fasta_to_dataframe(f"{seq_path}/fastas/{row['code']}.fasta")


@cli.command()
def get_invalid_files():
    """
    remove files that are not in the sequence sheet
    """
    setup_logging()
    df = get_sequence_sheet()
    seq_path = os.environ["SEQPATH"]
    rna_csvs = glob.glob(f"{seq_path}/rna/*.csv")
    dna_csvs = glob.glob(f"{seq_path}/dna/*.csv")
    fasta_files = glob.glob(f"{seq_path}/fastas/*.fasta")
    files = rna_csvs + dna_csvs + fasta_files
    f = open("README_REMOVE", "w")
    for file in files:
        code = Path(file).stem
        if code not in df["code"].values:
            log.info(f"code {code} not in sequence sheet")
            f.write(f'rm "{file}"\n')
    f.close()


if __name__ == "__main__":
    cli()
