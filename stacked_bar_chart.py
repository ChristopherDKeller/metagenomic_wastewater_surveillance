import os
import pandas as pd
import matplotlib.pyplot as plt

from util import FAMILY_COLOR_MAP, PLANT_NAME_MAP

# This script creates stacked bar charts of viral taxonomic compositions across samples.
# It can be configured by changing the constants below.

# ============================================================
# KONFIGURATION
# ============================================================
INPUT_FOLDER = "kraken2_run"
REPORTS_TO_USE = ["ERR12510709_report.txt", "ERR12510710_report.txt", "ERR12510724_report.txt", "ERR12510725_report.txt"] # [] = alle Reports im Ordner
TAXON_LEVEL = "F"
MIN_REL_ABUNDANCE = 0.03             # Taxa <x% werden zu "Other" zusammengefasst
META_CSV = "samples.csv"
# ============================================================

def parse_kraken2_report(path, taxon_level, sample_mapping):
    """
    Liest einen Kraken2-Report ein, filtert auf ein Taxonomie-Level
    und berechnet relative Häufigkeiten relativ zu allen Virus-Reads.
    """

    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=[
            "percent",
            "reads_clade",
            "reads_direct",
            "rank_code",
            "ncbi_taxid",
            "name"
        ],
        dtype={"name": str}
    )

    df["name"] = df["name"].str.strip()

    viruses_row = df[df["name"] == "Viruses"]
    if viruses_row.empty:
        raise RuntimeError(f"Keine Virus-Reads in Datei {path} gefunden.")
    else:
        virus_reads_total = viruses_row["reads_clade"].iloc[0]

    # taxonomisches Level filtern
    df_level = df[df["rank_code"] == taxon_level].copy()

    reads_at_level = df_level["reads_clade"].sum()
    unassigned_reads = virus_reads_total - reads_at_level
    
    # Erstelle eine Zeile für die nicht zugeordneten Viren
    unassigned_row = pd.DataFrame({ 
        "name": ["Unassigned"], 
        "rel": [unassigned_reads / virus_reads_total]
    })

    df_level["rel"] = df_level["reads_clade"] / virus_reads_total
    df_level = pd.concat([df_level[["name", "rel"]], unassigned_row], ignore_index=True)

    # Sample Label bestimmen
    basename = os.path.basename(path)
    run_id = basename.replace("_report.txt", "")

    if run_id in sample_mapping:
        sample_label = sample_mapping[run_id]
    else:
        sample_label = run_id

    df_level["sample"] = sample_label

    return df_level[["sample", "name", "rel"]]


def load_reports(input_folder, reports_to_use, taxon_level, sample_mapping):
    """
    Lädt alle Reports einzeln, normalisiert sie, und kombiniert erst danach.
    """
    all_files = sorted([f for f in os.listdir(input_folder) if f.endswith("_report.txt")])

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
        dfs.append(df_rel)

    return pd.concat(dfs, ignore_index=True)


def plot_stacked(df):
    """
    Erstellt einen gestapelten Balkenplot mit 'Other' und 'Unassigned'.
    """

    # Pivot: Zeilen = Samples, Spalten = Taxa, Werte = rel
    pivot = df.pivot_table(
        index="sample",
        columns="name",
        values="rel",
        fill_value=0
    )

    unassigned = pivot["Unassigned"]
    taxa_cols = pivot.columns.drop("Unassigned", errors="ignore")

    taxa_sums = pivot[taxa_cols].sum(axis=0)
    taxa_to_keep = taxa_sums[taxa_sums >= MIN_REL_ABUNDANCE].index
    taxa_to_other = taxa_cols.difference(taxa_to_keep)

    pivot_plot = pivot[taxa_to_keep].copy()
    pivot_plot["Other"] = pivot[taxa_to_other].sum(axis=1)

    if isinstance(unassigned, pd.Series):
        pivot_plot["Unassigned"] = unassigned

    main_cols = [c for c in pivot_plot.columns if c not in ("Other", "Unassigned")]

    sorted_main = (
        pivot_plot[main_cols]
        .sum(axis=0)
        .sort_values(ascending=False)
        .index
        .tolist()
    )

    final_cols = sorted_main
    if "Other" in pivot_plot.columns:
        final_cols.append("Other")
    if "Unassigned" in pivot_plot.columns:
        final_cols.append("Unassigned")

    pivot_plot = pivot_plot[final_cols]

    pivot_plot.plot(
        color=FAMILY_COLOR_MAP,
        kind="bar",
        stacked=True,
    )

    plt.ylabel("Relative Abundance")
    plt.title(f"Viral composition at taxonomic level '{TAXON_LEVEL}'")
    plt.legend(title="Taxon", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    plt.show()

def load_sample_metadata(csv_path):
    """
    Lädt die Metadaten-CSV und baut ein Mapping:
    RUN_ACCESSION → 'PLANT DATE (R1)'
    """
    df = pd.read_csv(csv_path, sep=";")

    mapping = {}
    for _, row in df.iterrows():
        run = str(row["ENA_RUN_ACCESSION"]).strip()
        plant = PLANT_NAME_MAP[str(row["PLANT"]).strip()]
        date = str(row["COLLECTION_DATE"]).strip()
        replica = row["REPLICA"]
        label = f"{plant} {date} (R{replica})"
        mapping[run] = label

    return mapping


if __name__ == "__main__":
    sample_mapping = load_sample_metadata(META_CSV)
    df = load_reports(INPUT_FOLDER, REPORTS_TO_USE, TAXON_LEVEL, sample_mapping)
    plot_stacked(df)
