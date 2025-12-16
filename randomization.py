import pandas as pd
import numpy as np
from scipy.spatial.distance import braycurtis
import matplotlib.pyplot as plt
import seaborn as sns

# ==========================
# KONFIGURATION
# ==========================
REPORT_FILE = "kraken2_run/ERR12510713_report.txt"
TAXON_LEVEL = "F"
N_ITER = 100           # Anzahl Randomization Iterationen
PLOT_HIST = True
# ==========================


def parse_kraken2_report(path, taxon_level):
    """
    Liest einen Kraken2-Report ein, filtert auf Taxonomie-Level und
    gibt ein Dictionary: Taxon → reads_clade zurück.
    """
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["percent", "reads_clade", "reads_direct", "rank_code", "ncbi_taxid", "name"],
        dtype={"name": str}
    )
    df["name"] = df["name"].str.strip()
    df_level = df[df["rank_code"] == taxon_level].copy()
    taxon_reads = df_level.set_index("name")["reads_clade"].to_dict()
    return taxon_reads


def randomize_similarity(taxon_reads, n_iter=100):
    """
    Teilt die Reads zufällig in zwei Hälften und berechnet Bray-Curtis-Similarity
    """
    # Erstelle ein Array aller Reads, jedes Taxon entsprechend seiner Reads
    all_reads = []
    for taxon, count in taxon_reads.items():
        all_reads.extend([taxon] * int(count))
    
    all_reads = np.array(all_reads)
    similarities = []

    for _ in range(n_iter):
        np.random.shuffle(all_reads)
        half = len(all_reads) // 2
        part1, part2 = all_reads[:half], all_reads[half:]

        # Zähle Reads pro Taxon
        unique_taxa = set(all_reads)
        counts1 = np.array([np.sum(part1 == t) for t in unique_taxa])
        counts2 = np.array([np.sum(part2 == t) for t in unique_taxa])

        # Bray-Curtis Similarity = 1 - Bray-Curtis Distance
        sim = 1 - braycurtis(counts1, counts2)
        similarities.append(sim)
    
    return similarities


if __name__ == "__main__":
    taxon_reads = parse_kraken2_report(REPORT_FILE, TAXON_LEVEL)
    sims = randomize_similarity(taxon_reads, n_iter=N_ITER)

    print(f"Randomized Bray-Curtis Similarities ({N_ITER} Iterationen):")
    print(f"  Mittelwert: {np.mean(sims):.4f}")
    print(f"  Median: {np.median(sims):.4f}")
    print(f"  Min: {np.min(sims):.4f}, Max: {np.max(sims):.4f}")

    if PLOT_HIST:
        ax = sns.histplot(sims)
        ax.set_xlabel("Bray-Curtis Similarity")
        ax.set_ylabel("Häufigkeit")
        ax.set_title("Randomized Similarities")

        plt.tight_layout()
        plt.show()

# Randomized Bray-Curtis Similarities (100 Iterationen) on Species Level:
#   Mittelwert: 0.9466
#   Median: 0.9465
#   Min: 0.9434, Max: 0.9508

# Randomized Bray-Curtis Similarities (100 Iterationen) on Genus Level:
#   Mittelwert: 0.9663
#   Median: 0.9663
#   Min: 0.9610, Max: 0.9699

# Randomized Bray-Curtis Similarities (100 Iterationen) on Family Level:
#   Mittelwert: 0.9902
#   Median: 0.9901
#   Min: 0.9862, Max: 0.9934
