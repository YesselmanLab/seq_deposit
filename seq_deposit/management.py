import cloup

from seq_deposit.logger import get_logger, setup_logging

log = get_logger("management")


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


@cloup.group()
def cli():
    """
    command line interface for seq_deposit script
    """
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
