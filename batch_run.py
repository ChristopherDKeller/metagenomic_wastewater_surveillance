import pandas as pd
import subprocess
import sys

# This script automates the download of FASTQ files from ENA, runs Kraken2 in a Docker container,
# and manages the output files for a batch of samples specified in a CSV file.

def main():
    if len(sys.argv) != 2:
        print("Verwendung:")
        print("  python batch_run.py <CSV_DATEI>")
        sys.exit(1)

    csv_path = sys.argv[1]

    print(f"→ Lese CSV: {csv_path}")
    df = pd.read_csv(csv_path, sep=";")

    if "ENA_RUN_ACCESSION" not in df.columns:
        raise ValueError("CSV hat keine Spalte 'ENA_RUN_ACCESSION'!")

    run_ids = df["ENA_RUN_ACCESSION"].dropna().unique()

    print(f"→ Gefundene {len(run_ids)} Runs\n")

    for run_accession in run_ids:
        print(f"====================================================")
        print(f" Starte Analyse für {run_accession}")
        print(f"====================================================")

        subprocess.run(
            ["python", "ena_kraken_automate.py", run_accession],
            check=True
        )

    print("\n✔ Alle Samples verarbeitet!\n")


if __name__ == "__main__":
    main()
