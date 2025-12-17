import os
import pandas as pd

INPUT_FOLDER = "kraken2_run"

# all possible phage taxa at various levels
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

phage_reads_total = 0
virus_reads_total = 0

for report_file in os.listdir(INPUT_FOLDER):
    if not report_file.endswith("_report.txt"):
        continue
    path = os.path.join(INPUT_FOLDER, report_file)
    df = parse_kraken2_report(path)

    viruses_row = df[df["name"] == "Viruses"]
    if viruses_row.empty:
        continue
    virus_reads = viruses_row["reads_clade"].iloc[0]
    virus_reads_total += virus_reads

    taxa_df = df[df["name"].isin(TAXA)]
    reads_total += taxa_df["reads_clade"].sum()

if virus_reads_total > 0:
    fraction = reads_total / virus_reads_total
    print(f"Anteil über alle Proben relativ zu allen Virus-Reads: {fraction:.2%}")
else:
    print("Keine viralen Reads gefunden.")
