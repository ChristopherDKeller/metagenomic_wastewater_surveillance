import os
import numpy as np
import pandas as pd

INPUT_FOLDER = "kraken2_run"

TAXA = {
    "Caudoviricetes",
    "Vidaverviricetes",
    "Leviviricetes",
    "Ainoaviricetes",
    "Picobirnaviridae",
    "Partitiviridae",
    "Matsushitaviridae",
    "Obscuriviridae",
    "Plasmaviridae",
    "Vinavirales"
}

def parse_kraken2_report(path):
    """Liest Kraken2-Report ein und gibt DataFrame aller Taxa zurück"""
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["percent", "reads_clade", "reads_direct", "rank_code", "ncbi_taxid", "name"],
        dtype={"name": str}
    )
    df["name"] = df["name"].str.strip()
    return df

fractions = []

for report_file in os.listdir(INPUT_FOLDER):
    if not report_file.endswith("_report.txt") or report_file in ["ERR2356165_report.txt", "ERR12510732_report.txt"]:
        continue
    path = os.path.join(INPUT_FOLDER, report_file)
    df = parse_kraken2_report(path)

    viruses_row = df[df["name"] == "Viruses"]
    if viruses_row.empty:
        continue
    virus_reads = viruses_row["reads_clade"].iloc[0]

    taxa_df = df[df["name"].isin(TAXA)]
    reads = taxa_df["reads_clade"].sum()
    fractions.append(reads / virus_reads)

print(f"Mittelwert über Proben: {np.mean(fractions):.2%}")
