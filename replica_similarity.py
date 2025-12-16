import os
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist
import matplotlib.pyplot as plt
import seaborn as sns

# ============================================================
# EINSTELLUNGEN
# ============================================================
INPUT_FOLDER = "kraken2_run"
TAXON_LEVEL = "F"
META_CSV = "samples.csv"
MIN_TOTAL_READS = 1_000_000
# ============================================================

def parse_kraken2_report(path, taxon_level):
    """
    Liest einen Kraken2-Report ein, filtert auf ein Taxonomie-Level
    und berechnet relative Häufigkeiten relativ zu allen Virus-Reads.
    """
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["percent", "reads_clade", "reads_direct", "rank_code", "ncbi_taxid", "name"],
        dtype={"name": str}
    )
    df["name"] = df["name"].str.strip()

    # QC-FILTER
    total_reads = df["reads_clade"].sum()

    if total_reads < MIN_TOTAL_READS:
        return None

    viruses_row = df[df["name"] == "Viruses"]
    if viruses_row.empty:
        raise RuntimeError(f"Keine Virus-Reads in Datei {path} gefunden.")
    virus_reads_total = viruses_row["reads_clade"].iloc[0]

    df_level = df[df["rank_code"] == taxon_level].copy()
    df_level["rel"] = df_level["reads_clade"] / virus_reads_total

    # Return: Series mit Index = Taxon, Werte = relative Häufigkeit
    return df_level.set_index("name")["rel"]

def load_metadata(meta_csv):
    df = pd.read_csv(meta_csv, sep=";")
    return df

def compute_replicate_similarity(metadata, input_folder, taxon_level):
    """
    Berechnet die Ähnlichkeit zwischen Replicates
    """
    # Zuerst alle Reports laden und relative Häufigkeiten merken
    rel_abundances = {}
    for _, row in metadata.iterrows():
        run = row["ENA_RUN_ACCESSION"]
        report_file = os.path.join(input_folder, f"{run}_report.txt")
        if os.path.exists(report_file):
            series = parse_kraken2_report(report_file, taxon_level)
            if series is None:
                print(f"⚠ Run {run} wegen QC ausgeschlossen")
                continue
            rel_abundances[run] = series
        else:
            print(f"⚠ Report für {run} nicht gefunden")

    similarities = []

    # Nach Probe gruppieren, dann Replicate-Paare bilden
    grouped = metadata.groupby("ENA_ALIAS")  # gleiche Probe = gleiche ENA_ALIAS
    for alias, group in grouped:
        runs = group["ENA_RUN_ACCESSION"].tolist()
        if len(runs) < 2:
            continue  # keine Replicate
        # Alle Paare bilden
        for i in range(len(runs)):
            for j in range(i+1, len(runs)):
                r1, r2 = runs[i], runs[j]
                if r1 in rel_abundances and r2 in rel_abundances:
                    # Union aller Taxa für Alignment
                    taxa = rel_abundances[r1].index.union(rel_abundances[r2].index)
                    vec1 = rel_abundances[r1].reindex(taxa, fill_value=0)
                    vec2 = rel_abundances[r2].reindex(taxa, fill_value=0)
                    # Bray-Curtis Similarity: 1 - Bray-Curtis Distance
                    bc_dist = pdist([vec1.values, vec2.values], metric="braycurtis")[0]
                    similarity = 1 - bc_dist
                    similarities.append({
                        "ENA_ALIAS": alias,
                        "Replicate1": r1,
                        "Replicate2": r2,
                        "Similarity": similarity
                    })
    return pd.DataFrame(similarities)


if __name__ == "__main__":
    metadata = load_metadata(META_CSV)
    sim_df = compute_replicate_similarity(metadata, INPUT_FOLDER, TAXON_LEVEL)
    print(sim_df)
    print(f"  Mittelwert: {np.mean(sim_df['Similarity']):.4f}")
    print(f"  Median: {np.median(sim_df['Similarity']):.4f}")
    print(f"  Min: {np.min(sim_df['Similarity']):.4f}, Max: {np.max(sim_df['Similarity']):.4f}")
    
    ax = sns.histplot(sim_df["Similarity"])
    ax.set_xlabel("Bray-Curtis Similarity")
    ax.set_ylabel("Häufigkeit")
    ax.set_title("Similarities between Replicates")

    plt.tight_layout()
    plt.show()

#  Species Level:
#  Mittelwert: 0.8940
#  Median: 0.8999
#  Min: 0.7874, Max: 0.9429

#  Genus Level:
#  Mittelwert: 0.9305
#  Median: 0.9398
#  Min: 0.8469, Max: 0.9660

#  Family Level:
#  Mittelwert: 0.9739
#  Median: 0.9762
#  Min: 0.9378, Max: 0.9900