import os
from matplotlib.patches import Patch
import matplotlib.dates as mdates
import pandas as pd
import matplotlib.pyplot as plt

from util import FAMILY_COLOR_MAP, GENUS_COLOR_MAP, ORDER_COLOR_MAP, PLANT_NAME_MAP

# This script creates stacked area plots of virus taxonomic levels over time for wastewater treatment plants.
# It can be configured by changing the constants below. 

# ============================================================
# KONFIGURATION
# ============================================================
INPUT_FOLDER = "kraken2_run"
REPORTS_TO_USE = []
TAXON_LEVEL = "G"               # "O"=Order, "F"=Family, "G"=Genus, None=all levels
MIN_REL_ABUNDANCE = 0.06
META_CSV = "samples.csv"
# ============================================================

def parse_kraken2_report(path, taxon_level):
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
    df_level["run"] = run_id
    return df_level[["run", "name", "rel"]]

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
        df_rel = parse_kraken2_report(path, taxon_level)
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
    pivot = plant_df.pivot_table(index=["DATE", "run"], columns="name", values="rel", aggfunc="sum").fillna(0)

    # Optional: nur Top-N Taxa
    if top_n is not None:
        top_taxa = pivot.sum().sort_values(ascending=False).head(top_n).index
        pivot = pivot[top_taxa]

    # Alle NaNs auf 0
    pivot = pivot.fillna(0)

    unassigned = pivot.get("Unassigned", pd.Series(0, index=pivot.index))
    remaining_cols = pivot.columns.drop("Unassigned", errors='ignore')

    taxa_to_keep = remaining_cols[pivot[remaining_cols].max(axis=0) >= min_rel_abundance]
    taxa_to_other = remaining_cols.difference(taxa_to_keep)

    pivot_plot = pivot[taxa_to_keep].copy()
    pivot_plot["Other"] = pivot[taxa_to_other].sum(axis=1)
    pivot_plot["Unassigned"] = unassigned

    #Replikate mitteln
    pivot_plot = pivot_plot.groupby(["DATE"]).mean()

    #sortieren nach Gesamtanteil
    cols_no_other = [c for c in pivot_plot.columns if c != "Other" and c != "Unassigned"]
    sorted_cols = (
        pivot_plot[cols_no_other]
        .sum(axis=0)
        .sort_values(ascending=False)
        .index
        .tolist()
    )
    sorted_cols.append("Other")
    sorted_cols.append("Unassigned")

    return pivot_plot[sorted_cols]

def plot_all_plants(df):
    plants = sorted(df["PLANT"].unique())
    n = len(plants)

    fig, axes = plt.subplots(n, 1, figsize=(12, 3*n), sharex=True)
    plt.subplots_adjust(hspace=0.4)
    if n == 1:
        axes = [axes]

    all_taxa = set()

    for ax, plant in zip(axes, plants):
        pivot_plot = prepare_time_series(df, plant)
        if pivot_plot is None or pivot_plot.empty:
            continue

        all_taxa.update(pivot_plot.columns)

        if TAXON_LEVEL == "O":
            colors = [ORDER_COLOR_MAP.get(t, "#BBBBBB") for t in pivot_plot.columns]
        elif TAXON_LEVEL == "F":
            colors = [FAMILY_COLOR_MAP.get(t, "#BBBBBB") for t in pivot_plot.columns]
        elif TAXON_LEVEL == "G":
            colors = [GENUS_COLOR_MAP.get(t, "#BBBBBB") for t in pivot_plot.columns]
        else:
            colors = None

        pivot_plot.plot.area(
            ax=ax,
            color= colors if colors else None,
            legend=False
        )

        ax.set_title(PLANT_NAME_MAP.get(plant, plant))
        ax.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[3,6,9,12], interval=1))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))

        ax.xaxis.set_minor_locator(mdates.MonthLocator())


    axes[int(n/2)].set_ylabel("Relative Abundance")
    axes[-1].set_xlabel("Date")

    ordered_taxa_O = ["Crassvirales", "Timlovirales", "Chitovirales", "Imitervirales","Herpesvirales", "Tubulavirales", "Lefavirales", "Bunyavirales","Pimascovirales", "Algavirales", "Halopanivirales","Mononegavirales", "Rowavirales", "Picornavirales","Other"]

    ordered_taxa_F = ["Intestiviridae","Suoliviridae","Peduoviridae","Crevaviridae","Herelleviridae","Kyanoviridae","Straboviridae","Schitoviridae","Steitzviridae","Steigviridae","Autographiviridae","Demerecviridae","Mimiviridae","Poxviridae","Arenbergviridae","Other"]

    ordered_taxa_G = [ "Carjivirus", "Burzaovirus", "Punavirus", "Agtrevirus", "Betabaculovirus", "Gihfavirus", "Casadabanvirus", "Purivirus", "Baikalvirus", "Pamexvirus", "Immutovirus", "Other", "Unassigned" ]

    if TAXON_LEVEL == "O":
        order = ordered_taxa_O
    elif TAXON_LEVEL == "F":
        order = ordered_taxa_F
    elif TAXON_LEVEL == "G":
        order = ordered_taxa_G

    legend_taxa = [t for t in order if t in all_taxa]
    
    add_global_legend(axes[-1], legend_taxa)
    plt.show()


def add_global_legend(ax, taxa):
    if TAXON_LEVEL == "O":
        color_map = ORDER_COLOR_MAP
    elif TAXON_LEVEL == "F":
        color_map = FAMILY_COLOR_MAP
    elif TAXON_LEVEL == "G":
        color_map = GENUS_COLOR_MAP

    handles = [
        Patch(facecolor=color_map.get(t, "#BBBBBB"), label=t)
        for t in taxa
    ]

    ax.legend(
        handles=handles,
        title="Viral order",
        loc="lower left",
        bbox_to_anchor=(1.01, 3)
    )

def plot_single_plant(pivot_plot, plant):
    if TAXON_LEVEL == "O":
        colors = [ORDER_COLOR_MAP.get(taxon, "#BBBBBB") for taxon in pivot_plot.columns]
    elif TAXON_LEVEL == "F":
        colors = [FAMILY_COLOR_MAP.get(taxon, "#BBBBBB") for taxon in pivot_plot.columns]
    elif TAXON_LEVEL == "G":
        colors = [GENUS_COLOR_MAP.get(taxon, "#BBBBBB") for taxon in pivot_plot.columns]
    else:
        colors = None

    pivot_plot.plot.area(color=colors if colors else None)
    plt.title(f"Relative Virus-Häufigkeiten über Zeit ({plant})")
    plt.ylabel("Relative Häufigkeit")
    plt.xlabel("Datum")
    plt.legend(title="Order", bbox_to_anchor=(1.05,1), loc="upper left")
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    sample_mapping = load_sample_metadata(META_CSV)
    df = load_reports(INPUT_FOLDER, REPORTS_TO_USE, TAXON_LEVEL, sample_mapping, reports_to_skip=["ERR2356165_report.txt", "ERR12510732_report.txt"])

    plot_all_plants(df)

    #for plant in df["PLANT"].unique():
    #    pivot_plot = prepare_time_series(df, plant)
    #    plot_single_plant(pivot_plot, plant)
