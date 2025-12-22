import os
import pandas as pd
from scipy.spatial.distance import pdist
import matplotlib.pyplot as plt
import seaborn as sns

# ============================================================
# EINSTELLUNGEN
# ============================================================
INPUT_FOLDER = "kraken2_run"
META_CSV = "samples.csv"

TAXON_LEVEL = "O"
MODE = "temporal"           # "replicate" oder "temporal"

MIN_TOTAL_READS = 1_000_000
DATE_COLUMN = "COLLECTION_DATE"
PLANT_COLUMN = "PLANT"
# ============================================================


def parse_kraken2_report(path, taxon_level):
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["percent", "reads_clade", "reads_direct",
               "rank_code", "ncbi_taxid", "name"],
        dtype={"name": str}
    )
    df["name"] = df["name"].str.strip()

    total_reads = df["reads_clade"].sum()
    if total_reads < MIN_TOTAL_READS:
        return None

    viruses_row = df[df["name"] == "Viruses"]
    if viruses_row.empty:
        return None

    virus_reads_total = viruses_row["reads_clade"].iloc[0]

    if taxon_level:
        df = df[df["rank_code"] == taxon_level]

    df["rel"] = df["reads_clade"] / virus_reads_total

    return df.set_index("ncbi_taxid")["rel"]


def load_all_profiles(metadata):
    rel_abundances = {}

    for _, row in metadata.iterrows():
        run = row["ENA_RUN_ACCESSION"]
        report = os.path.join(INPUT_FOLDER, f"{run}_report.txt")

        if not os.path.exists(report):
            continue

        series = parse_kraken2_report(report, TAXON_LEVEL)
        if series is not None:
            rel_abundances[run] = series

    return rel_abundances


def bray_curtis_similarity(s1, s2):
    taxa = s1.index.union(s2.index)
    v1 = s1.reindex(taxa, fill_value=0)
    v2 = s2.reindex(taxa, fill_value=0)

    return 1 - pdist([v1.values, v2.values], metric="braycurtis")[0]


# ============================================================
# VERGLEICHSMODI
# ============================================================

def compare_replicates(metadata, profiles):
    results = []

    for alias, group in metadata.groupby("ENA_ALIAS"):
        runs = group["ENA_RUN_ACCESSION"].tolist()
        if len(runs) < 2:
            continue

        for i in range(len(runs)):
            for j in range(i + 1, len(runs)):
                r1, r2 = runs[i], runs[j]
                if r1 in profiles and r2 in profiles:
                    sim = bray_curtis_similarity(profiles[r1], profiles[r2])
                    results.append({
                        "Type": "Replicate",
                        "Group": alias,
                        "Run1": r1,
                        "Run2": r2,
                        "Similarity": sim
                    })

    return pd.DataFrame(results)


def compare_temporal(metadata, profiles):
    results = []

    metadata[DATE_COLUMN] = pd.to_datetime(metadata[DATE_COLUMN])

    for plant, group in metadata.groupby(PLANT_COLUMN):
        group = group.sort_values(DATE_COLUMN)

        runs = group["ENA_RUN_ACCESSION"].tolist()
        dates = group[DATE_COLUMN].tolist()

        for i in range(len(runs) - 1):
            r1, r2 = runs[i], runs[i + 1]
            if r1 in profiles and r2 in profiles:
                sim = bray_curtis_similarity(profiles[r1], profiles[r2])
                delta_days = (dates[i + 1] - dates[i]).days

                results.append({
                    "Type": "Temporal",
                    "Plant": plant,
                    "Run1": r1,
                    "Run2": r2,
                    "Delta_days": delta_days,
                    "Similarity": sim
                })

    return pd.DataFrame(results)


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    metadata = pd.read_csv(META_CSV, sep=";")
    profiles = load_all_profiles(metadata)

    if MODE == "replicate":
        sim_df = compare_replicates(metadata, profiles)
        title = "Bray-Curtis Similarity (Technical Replicates)"

    elif MODE == "temporal":
        sim_df = compare_temporal(metadata, profiles)
        title = "Bray-Curtis Similarity (Adjacent Timepoints)"
        summary = (
            sim_df
            .groupby("Plant")
            .agg(
                mean_similarity=("Similarity", "mean"),
                median_similarity=("Similarity", "median"),
                min_similarity=("Similarity", "min"),
                max_similarity=("Similarity", "max"),
            )
        )
        print(summary.round(4))

    else:
        raise ValueError("MODE must be 'replicate' or 'temporal'")

    #print(sim_df.describe())
    print(f"Mittelwert: {sim_df['Similarity'].mean():.4f}")
    print(f"Median:     {sim_df['Similarity'].median():.4f}")
    print(f"Min:        {sim_df['Similarity'].min():.4f}")
    print(f"Max:        {sim_df['Similarity'].max():.4f}")

    ax = sns.histplot(sim_df["Similarity"])
    ax.set_xlabel("Bray-Curtis Similarity")
    ax.set_ylabel("Count")
    ax.set_title(title)

    plt.tight_layout()
    plt.show()


# Results replica similarity:
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

#  Overall:
#  Mittelwert: 0.9867
#  Median: 0.9866
#  Min: 0.9668, Max: 0.9986

# Results temporal similarity:
#  Species Level:
#  Mittelwert: 0.7824
#  Median:     0.8226
#  Min:        0.2438
#  Max:        0.9452

#  Genus Level:
#  Mittelwert: 0.8121
#  Median:     0.8482
#  Min:        0.2759
#  Max:        0.9671

#  Oder Level:
#  Mittelwert: 0.8603
#  Median:     0.9064
#  Min:        0.2842
#  Max:        0.9961

#  Overall:
#  Mittelwert: 0.8745
#  Median:     0.9207
#  Min:        0.3220
#  Max:        0.9980