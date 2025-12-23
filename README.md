# Metagenomic Wastewater Surveillance

## Overview

This project automates the analysis of viral sequences from wastewater samples using Kraken2. The pipeline downloads FASTQ data from ENA, performs taxonomic classification, and analyzes viral composition over time and across different wastewater treatment plants.

## Components

### Pipeline Automation

**`batch_run.py`**

- Automates execution of multiple samples
- Reads ENA Run Accessions from CSV file
- Triggers `ena_kraken_automate.py` for each sample

**`ena_kraken_automate.py`**

- Downloads FASTQ files from ENA
- Runs Kraken2 in Docker container
- Automatic cleanup of raw data and output files

### Analysis Modules

**`plant_similarity.py`**

- Computes Bray-Curtis similarity between treatment plants
- Aggregates viral profiles per plant
- Visualizes similarity matrices as heatmap

**`similarity.py`**

- Compares similarity between technical replicates or temporally adjacent samples
- Modes: "replicate" or "temporal"
- Calculates statistics (mean, median, min/max)

**`randomization.py`**

- Randomly partitions reads and computes Bray-Curtis similarity

**`proportion.py`**

- Calculates proportions of specific viral taxa

### Visualization

**`stacked_bar_chart.py`**

- Creates stacked bar charts of viral composition

**`zeitreihe.py`**

- Creates stacked area charts over time
- Separate plots per treatment plant
- Supports various taxonomic levels (Order, Family, Genus)
- Averages replicates by date

**`util.py`**

- Central configuration for color and name mappings
- `PLANT_NAME_MAP`: Treatment plant aliases
- `FAMILY_COLOR_MAP`, `GENUS_COLOR_MAP`, `ORDER_COLOR_MAP`: Visualization colors

## Workflow

1. Prepare CSV with ENA Run Accessions
2. Execute `batch_run.py`
3. Pipeline downloads data and runs Kraken2
4. Reports are saved in `kraken2_run/`
5. Run analysis scripts on the reports

## Requirements

- Python 3.x
- see `requirements.txt`
