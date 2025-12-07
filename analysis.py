import os
import pandas as pd
import matplotlib.pyplot as plt

# ============================================================
# KONFIGURATION
# ============================================================
INPUT_FOLDER = "kraken2_run"         # Ordner mit Kraken2-Reports
REPORTS_TO_USE = ["ERR12510713_report.txt", "ERR12510871_report.txt"] # [] = alle Reports im Ordner
TAXON_LEVEL = "F"                    # "F"=Familie, "G"=Genus, "S"=Species ...
MIN_REL_ABUNDANCE = 0.03             # Taxa <x% werden zu "Other" zusammengefasst
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

    # Virus-Gesamtreads
    viruses_row = df[df["name"] == "Viruses"]

    if viruses_row.empty:
        raise RuntimeError("Keine Virus-Reads in Datei {path} gefunden.")
    else:
        virus_reads_total = viruses_row["reads_clade"].iloc[0]

    # taxonomisches Level filtern
    df_level = df[df["rank_code"] == taxon_level].copy()

    # relative Häufigkeiten im Verhältnis zu allen Viren
    df_level["rel"] = df_level["reads_clade"] / virus_reads_total
    
    # Probenname
    df_level["sample"] = os.path.basename(path)

    return df_level[["sample", "name", "rel"]]


def load_reports(input_folder, reports_to_use, taxon_level):
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
        df_rel = parse_kraken2_report(path, taxon_level)
        dfs.append(df_rel)

    return pd.concat(dfs, ignore_index=True)


def plot_stacked(df):
    """
    Erstellt einen gestapelten Balkenplot mit einer 'Other'-Kategorie.
    """

    # Pivot: Zeilen = Samples, Spalten = Taxa, Werte = rel
    pivot = df.pivot_table(
        index="sample",
        columns="name",
        values="rel",
        fill_value=0
    )

    # Taxa auswählen, die global >= MIN_REL_ABUNDANCE sind
    taxa_to_keep = pivot.columns[pivot.sum(axis=0) >= MIN_REL_ABUNDANCE]
    taxa_to_other = pivot.columns.difference(taxa_to_keep)

    pivot_plot = pivot[taxa_to_keep].copy()

    # "Other" spalte hinzufügen
    pivot_plot["Other"] = pivot[taxa_to_other].sum(axis=1)

    # Jede Zeile sauber auf 1 normieren
    pivot_plot = pivot_plot.div(pivot_plot.sum(axis=1), axis=0)

    pivot_plot.plot(
        kind="bar",
        stacked=True,
    )

    plt.ylabel("Relative Abundance (within all viruses)")
    plt.title(f"Viral composition at taxonomic level '{TAXON_LEVEL}'")
    plt.legend(title="Taxon", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    df = load_reports(INPUT_FOLDER, REPORTS_TO_USE, TAXON_LEVEL)
    plot_stacked(df)
