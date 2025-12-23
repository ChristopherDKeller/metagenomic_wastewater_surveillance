import os
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import braycurtis

from util import PLANT_NAME_MAP

# This script calculaes the Bray-Curtis similarity between treatment plants on the basis of viral taxonomic profiles.
# The taxonomic profiles are aggregated per plant from Kraken2 reports.
# It can be configured by changing the constants below.

# ============================================================
# KONFIGURATION
# ============================================================
REPORT_DIR = "kraken2_run"
META_CSV = "samples.csv"
TAXON_LEVEL = None
# ============================================================


def parse_kraken2_report(path, taxon_level):
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["percent", "reads_clade", "reads_direct",
               "rank_code", "taxid", "name"],
        dtype={"name": str}
    )
    df["name"] = df["name"].str.strip()

    virus_row = df[df["name"] == "Viruses"]
    if virus_row.empty:
        return None

    virus_reads = virus_row["reads_clade"].iloc[0]

    if taxon_level is not None:
        df = df[df["rank_code"] == taxon_level].copy()

    df_direct = df[df["reads_direct"] > 0].copy()
    df_direct["rel"] = df_direct["reads_clade"] / virus_reads

    return df_direct.set_index("name")["rel"]


def load_metadata(csv_path):
    df = pd.read_csv(csv_path, sep=";")
    return df


def build_sample_matrix(metadata):
    rows = []
    meta_rows = []

    for _, row in metadata.iterrows():
        run = row["ENA_RUN_ACCESSION"]
        report = os.path.join(REPORT_DIR, f"{run}_report.txt")

        if not os.path.exists(report):
            continue

        rel = parse_kraken2_report(report, TAXON_LEVEL)
        if rel is None:
            continue

        rows.append(rel)
        meta_rows.append(row)

    df = pd.DataFrame(rows).fillna(0)
    meta_df = pd.DataFrame(meta_rows).reset_index(drop=True)

    return df, meta_df


def aggregate_by_plant(df, meta_df):
    df["PLANT"] = meta_df["PLANT"].values
    plant_profiles = df.groupby("PLANT").mean()
    return plant_profiles


def compute_similarity_matrix(plant_profiles):
    plants = plant_profiles.index
    sim = pd.DataFrame(index=plants, columns=plants, dtype=float)

    for p1 in plants:
        for p2 in plants:
            bc = braycurtis(plant_profiles.loc[p1], plant_profiles.loc[p2])
            sim.loc[p1, p2] = 1 - bc

    return sim

def plot_similarity_heatmap(similarity_df):
    fig, ax = plt.subplots(figsize=(6, 5))

    im = ax.imshow(
        similarity_df.values,
        vmin=0,
        vmax=1,
        aspect="auto"
    )

    # Achsenbeschriftungen
    ax.set_xticks(np.arange(len(similarity_df.columns)))
    ax.set_yticks(np.arange(len(similarity_df.index)))

    ax.set_xticklabels((PLANT_NAME_MAP.get(col, col) for col in similarity_df.columns), rotation=45, ha="right")
    ax.set_yticklabels((PLANT_NAME_MAP.get(row, row) for row in similarity_df.index))

    # Farbskala
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("Bray-Curtis Similarity", rotation=90)

    ax.set_title("Viral Similarity between Treatment Plants")

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    metadata = load_metadata(META_CSV)
    abundance, meta_df = build_sample_matrix(metadata)

    plant_profiles = aggregate_by_plant(abundance, meta_df)
    similarity = compute_similarity_matrix(plant_profiles)

    print("\nBray-Curtis Similarity zwischen Klärwerken:\n")
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):
        print(similarity.round(3))
    plot_similarity_heatmap(similarity)
