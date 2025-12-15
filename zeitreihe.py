import os
import pandas as pd
import matplotlib.pyplot as plt

from colormaps import ORDER_COLOR_MAP

# ============================================================
# KONFIGURATION
# ============================================================
INPUT_FOLDER = "kraken2_run"
REPORTS_TO_USE = []
TAXON_LEVEL = "O"
MIN_REL_ABUNDANCE = 0.05
META_CSV = "samples.csv"
# ============================================================

def parse_kraken2_report(path, taxon_level, sample_mapping):
    df = pd.read_csv(
        path, sep="\t", header=None,
        names=["percent", "reads_clade", "reads_direct", "rank_code", "ncbi_taxid", "name"],
        dtype={"name": str}
    )

    df["name"] = df["name"].str.strip()

    # Virus-Reads
    viruses_row = df[df["name"] == "Viruses"]
    if viruses_row.empty:
        raise RuntimeError(f"Keine Virus-Reads in Datei {path} gefunden.")
    virus_reads_total = viruses_row["reads_clade"].iloc[0]

    # Filter Taxonomie-Level
    df_level = df[df["rank_code"] == taxon_level].copy()

    # relative Häufigkeiten
    df_level["rel"] = df_level["reads_clade"] / virus_reads_total

    # Sample Label bestimmen
    basename = os.path.basename(path)
    run_id = basename.replace("_report.txt", "")
    if run_id in sample_mapping:
        sample_label = sample_mapping[run_id]
    else:
        sample_label = run_id

    df_level["sample"] = sample_label
    df_level["run"] = run_id
    return df_level[["sample", "run", "name", "rel"]]

def load_sample_metadata(csv_path):
    """
    Lädt die Metadaten-CSV und baut ein Mapping:
    RUN_ACCESSION → {'DATE': ..., 'PLANT': ...}
    """
    df = pd.read_csv(csv_path, sep=";")
    mapping = {}
    for _, row in df.iterrows():
        run = str(row["ENA_RUN_ACCESSION"]).strip()
        date = pd.to_datetime(row["COLLECTION_DATE"])
        plant = str(row["PLANT"]).strip()
        mapping[run] = {"DATE": date, "PLANT": plant}
    return mapping

def load_reports(input_folder, reports_to_use, taxon_level, sample_mapping, reports_to_skip=None):
    all_files = sorted([f for f in os.listdir(input_folder) if f.endswith("_report.txt")])

    if reports_to_skip:
        all_files = [f for f in all_files if f not in reports_to_skip]

    if reports_to_use:
        selected = [f for f in all_files if f in reports_to_use]
    else:
        selected = all_files

    if not selected:
        raise RuntimeError("Keine passenden Reports gefunden!")

    dfs = []
    for report in selected:
        path = os.path.join(input_folder, report)
        print(f"Lade {path}")
        df_rel = parse_kraken2_report(path, taxon_level, sample_mapping)
        # Metadaten hinzufügen
        meta = sample_mapping[df_rel["run"].iloc[0]]
        df_rel["DATE"] = meta["DATE"]
        df_rel["PLANT"] = meta["PLANT"]
        dfs.append(df_rel)
    return pd.concat(dfs, ignore_index=True)

def prepare_time_series(df, plant, min_rel_abundance=MIN_REL_ABUNDANCE, top_n=None):
    """Stacked area plot der Virusfamilien über Zeit für eine Kläranlage"""
    plant_df = df[df["PLANT"] == plant].copy()
    if plant_df.empty:
        print(f"⚠ Keine Daten für Kläranlage {plant}")
        return

    # Pivot: DATE x Taxon
    pivot = plant_df.pivot_table(index="DATE", columns="name", values="rel", aggfunc="sum").fillna(0)

    # Optional: nur Top-N Taxa
    if top_n is not None:
        top_taxa = pivot.sum().sort_values(ascending=False).head(top_n).index
        pivot = pivot[top_taxa]

    # Alle NaNs auf 0
    pivot = pivot.fillna(0)

    # Taxa nach globaler Summe filtern
    taxa_to_keep = pivot.columns[pivot.sum(axis=0) >= min_rel_abundance]
    taxa_to_other = pivot.columns.difference(taxa_to_keep)

    # Pivot für Plot
    pivot_plot = pivot[taxa_to_keep].copy()
    pivot_plot["Other"] = pivot[taxa_to_other].sum(axis=1)

    pivot_plot = pivot_plot.div(pivot_plot.sum(axis=1), axis=0)

    #Replikate mitteln
    pivot_plot = pivot_plot.groupby(["DATE"]).mean()

    #sortieren nach Gesamtanteil
    cols_no_other = [c for c in pivot_plot.columns if c != "Other"]
    sorted_cols = (
        pivot_plot[cols_no_other]
        .sum(axis=0)
        .sort_values(ascending=False)
        .index
        .tolist()
    )
    sorted_cols.append("Other")

    return pivot_plot[sorted_cols]

def plot_all_plants(df):
    plants = sorted(df["PLANT"].unique())
    n = len(plants)

    fig, axes = plt.subplots(n, 1, figsize=(12, 3*n), sharex=True)

    if n == 1:
        axes = [axes]

    for ax, plant in zip(axes, plants):
        pivot_plot = prepare_time_series(df, plant)
        if pivot_plot is None:
            continue

        colors = [ORDER_COLOR_MAP.get(taxon, "#BBBBBB") for taxon in pivot_plot.columns]
        pivot_plot.plot.area(ax=ax, color=colors, legend=False)
        ax.set_title(plant)
        ax.set_ylabel("Rel. Abundance")

    axes[-1].set_xlabel("Date")
    plt.tight_layout()
    plt.show()

def plot_single_plant(pivot_plot, plant):
    colors = [ORDER_COLOR_MAP.get(taxon, "#BBBBBB") for taxon in pivot_plot.columns]

    pivot_plot.plot.area(color=colors)
    plt.title(f"Relative Virus-Häufigkeiten über Zeit ({plant})")
    plt.ylabel("Relative Häufigkeit")
    plt.xlabel("Datum")
    plt.legend(title="Order", bbox_to_anchor=(1.05,1), loc="upper left")
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    sample_mapping = load_sample_metadata(META_CSV)
    df = load_reports(INPUT_FOLDER, REPORTS_TO_USE, TAXON_LEVEL, sample_mapping, reports_to_skip=["ERR2356165_report.txt"])

    plot_all_plants(df)

    #for plant in df["PLANT"].unique():
        pivot_plot = prepare_time_series(df, plant)
        plot_single_plant(pivot_plot, plant)
