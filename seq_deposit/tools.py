import os
import pandas as pd

import vienna
from seq_deposit import settings
from seq_tools import sequence


def main():
    seq_path = os.environ["SEQPATH"]
    path = settings.LIB_PATH + "/seq_deposit/resources/all.csv"
    df = pd.read_csv(path)
    # df = df[pd.isnull(df["rna_sequence"])]
    df = df[df["usuable"] != "NO"]
    df = df[df["type"] == "ASSEMBLY"]
    cols = "rna_len,rna_sequence,rna_structure,fasta_link,comment".split(",")
    for i, row in df.iterrows():
        dna_seq = row["dna_sequence"]
        if not dna_seq.startswith("TTCTAATACGACTCACTATA"):
            print(row["name"] + " does not contain a t7 promoter")
            continue
        else:
            dna_seq = dna_seq[20:]
        rna_seq = sequence.convert_to_rna(dna_seq)
        rna_ss = vienna.folded_structure(rna_seq)
        fasta_path = f"{seq_path}/fastas/{row['code']}.fasta"
        csv_path = f"{seq_path}/rna/{row['code']}.csv"
        if not os.path.isfile(fasta_path):
            print(row["name"])
            f = open(fasta_path, "w")
            f.write(f">{row['name']}\n")
            f.write(dna_seq)
            f.close()
        if not os.path.isfile(csv_path):
            f = open(csv_path, "w")
            f.write("name,sequence,structure\n")
            f.write(f"{row['name']},{rna_seq},{rna_ss}\n")
            f.close()
        # print(f"{len(rna_seq)},{rna_seq},{rna_ss},$SEQPATH/fastas/{row['code']}.fasta")

    # s = f"{}"


if __name__ == "__main__":
    main()
