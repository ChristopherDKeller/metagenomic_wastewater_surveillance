import os
import subprocess
import requests
import sys
import shutil
import time

KRAKEN2_IMAGE = "staphb/kraken2:2.1.6-viral-20250402"
OUTPUT_DIR = "kraken2_run"
THREADS = 4


def ena_fastq_links(run_accession):
    url = (
        "https://www.ebi.ac.uk/ena/portal/api/filereport"
        f"?accession={run_accession}"
        "&result=read_run"
        "&fields=fastq_ftp"
        "&format=tsv"
        "&download=true"
    )

    print(f"→ Anfrage an ENA: {url}")

    response = requests.get(url)
    if response.status_code != 200:
        raise RuntimeError(f"ENA API Fehler: {response.status_code}")

    lines = response.text.strip().split("\n")
    if len(lines) < 2:
        raise RuntimeError("Keine FASTQ-Dateien gefunden!")

    fastq_ftp = lines[1].split("\t")[1]
    fastq_files = fastq_ftp.split(";")
    fastq_https = ["https://" + f for f in fastq_files]

    return fastq_https

def download_fastqs(urls, output_dir, max_retries=1000, chunk_size=1024*1024):
    os.makedirs(output_dir, exist_ok=True)
    local_paths = []

    for url in urls:
        filename = url.split("/")[-1]
        local_path = os.path.join(output_dir, filename)
        local_paths.append(local_path)

        if os.path.exists(local_path):
            print(f"→ Datei existiert bereits, überspringe: {filename}")
            continue

        print(f"→ Lade herunter: {filename}")

        for attempt in range(1, max_retries+1):
            try:
                headers = {}
                # Resume if file partially exists
                if os.path.exists(local_path):
                    downloaded = os.path.getsize(local_path)
                    headers["Range"] = f"bytes={downloaded}-"
                else:
                    downloaded = 0

                with requests.get(url, headers=headers, stream=True, timeout=60) as r:
                    r.raise_for_status()

                    mode = "ab" if downloaded > 0 else "wb"

                    with open(local_path, mode) as f:
                        for chunk in r.iter_content(chunk_size=chunk_size):
                            if chunk:
                                f.write(chunk)

                print(f"✔ Download abgeschlossen: {filename}")
                break

            except Exception as e:
                print(f"⚠ Download-Fehler (Versuch {attempt}/{max_retries}): {e}")
                time.sleep(3)

                if attempt == max_retries:
                    raise RuntimeError(f"❌ Abbruch nach {max_retries} Fehlversuchen bei {filename}")

    return local_paths

def run_kraken2(docker_image, run_accession, fastq_files, output_dir, threads):
    docker_mount = os.path.abspath(output_dir)
    fastq_inside = [os.path.basename(f) for f in fastq_files]

    output_report = os.path.join(output_dir, f"{run_accession}_report.txt")
    output_output = os.path.join(output_dir, f"{run_accession}_output.txt")

    kraken_cmd = [
        "docker", "run", "--rm",
        "-v", f"{docker_mount}:/data",
        docker_image,
        "kraken2",
        "--db", "/kraken2-db",
        "--threads", str(threads),
        "--report", f"/data/{run_accession}_report.txt",
        "--output", f"/data/{run_accession}_output.txt",
    ]

    if len(fastq_inside) == 2:
        kraken_cmd += ["--paired"] + fastq_inside
    else:
        kraken_cmd += fastq_inside

    print("\n→ Starte Kraken2 im Docker-Container:")
    print(" ".join(kraken_cmd))
    print()

    subprocess.run(kraken_cmd, check=True)

    print("\n✔ Kraken2 abgeschlossen")
    print(f"  Report: {output_report}")
    print(f"  Output: {output_output}")

    # -------------------------------------------------------------
    # CLEANUP: FASTQ-Dateien + Kraken-Output löschen
    # -------------------------------------------------------------
    print("→ Cleanup...")
    for fq in fastq_files:
        try:
            os.remove(fq)
            print(f"  gelöscht: {fq}")
        except Exception as e:
            print(f"  konnte {fq} nicht löschen: {e}")

    try:
        os.remove(output_output)
        print(f"  gelöscht: {output_output}")
    except Exception as e:
        print(f"  konnte {output_output} nicht löschen: {e}")
    print("✔ Cleanup abgeschlossen")


def run_pipeline(run_accession):
    print(f"\n=== Kraken2 Pipeline für ENA Run {run_accession} ===\n")

    report_path = os.path.join(OUTPUT_DIR, f"{run_accession}_report.txt")

    if os.path.exists(report_path):
        print(f"✔ Report existiert bereits: {report_path}")
    else:
        fastq_urls = ena_fastq_links(run_accession)
        print(f"→ Gefundene FASTQ-Dateien: {fastq_urls}")

        fastq_paths = download_fastqs(fastq_urls, OUTPUT_DIR)
        print(f"→ Downloads gespeichert in {OUTPUT_DIR}")

        run_kraken2(KRAKEN2_IMAGE, run_accession, fastq_paths, OUTPUT_DIR, THREADS)

    print("\n=== Fertig! ===\n")


def main():
    if len(sys.argv) != 2:
        print("Verwendung:")
        print("  python kraken_single.py <ENA_RUN_ACCESSION>")
        sys.exit(1)

    run_accession = sys.argv[1]
    run_pipeline(run_accession)


if __name__ == "__main__":
    main()
