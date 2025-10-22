# Boltz-2 Analysis for Hayward, et al 2025

This repository contains code used for the analysis of Boltz-2 structures generated for the article "Tryptoline Stereoprobe Elaboration Identifies Inhibitors of the GRPEL1-HSPA9 Chaperone Complex". The [preprint](https://www.biorxiv.org/content/10.1101/2025.10.20.683548v1) is available on bioRxiv.

## Installation

First, clone this repository:

```bash
git clone https://github.com/cravattlab/Hayward_2025_Boltz_analysis.git
cd Hayward_2025_Boltz_analysis
```

We recommend installing Python dependencies for the scripts in this repository using [uv](https://github.com/astral-sh/uv) or [conda](https://github.com/conda/conda):

With `uv`:
```bash
# Create and activate a virtual environment
uv venv .venv
source .venv/bin/activate

# Install dependencies using requirements.txt
uv pip install -r requirements.txt
```

With `conda`:
```bash
# Create and activate a virtual env using environment.yml
conda env create -f environment.yml
conda activate hayward-boltz2
```

[Phenix](https://www.phenix-online.org) (v1.21) is the sole external dependency. The Phenix and PoseBusters validation steps, which account for the majority of the runtime of the analysis pipeline, can optionally be skipped by using the `--skip-validation` argument.

Installation of all required dependencies should take less than 5 minutes; analysis takes on the order of 10 sec per structure with Phenix validation as tested on consumer-grade hardware (M1 Max, macOS Tahoe 26.0 or Intel Xeon E5-2660, CentOS 7).

## Recommendations for running the analysis pipeline

The modules in this repository are called from the main `process_boltz_from_index.py` script, which performs initial parsing of Boltz-2 cif files as well as all distance calculations. The script additionally orchestrates the application of orthosteric site analysis, orthostery classification for ligands and cysteines, stereochemical evaluation, and structure validation with PoseBusters and Phenix.

The `process_boltz_from_index.py` script expects a directory of cif files that are named according to `{accession}_{arbitrary_field}_{ligand_name}_model_{num}.cif`, along with an index file that specifies corresponding UniProt accession numbers and ligand names. For the stereochemistry analysis, an additional sheet in the index file containing ligand SMILES is required. The index file used in the manuscript associated with this repository is provided as Supplementary Dataset 1.

After activating the appropriate virtual environment, an example processing command might be:

```bash
(hayward_boltz2) ➜ ~ python -u process_boltz_from_index.py --predictions-dir=./example_data/ --index-file=./Supplementary_Dataset_1.xlsx --index-sheet='liganded sites' --ligand-sheet='compounds' --output=example_output.xlsx --index-header=0 --ligands-header=0 | tee example_output.log
```

All additional arguments available for adjusting the behavior of this analysis pipeline can be found in the `--help` option:

```bash
(hayward_boltz2) ➜ ~ python process_boltz_from_index.py --help
usage: process_boltz_from_index.py [-h] --predictions-dir PREDICTIONS_DIR [--output OUTPUT]
                                   --index-file INDEX_FILE --index-sheet INDEX_SHEET
                                   [--ligand-sheet LIGAND_SHEET] [--site-cutoff SITE_CUTOFF]
                                   [--filter-cys-domains] [--include-extended-fields]
                                   [--skip-uniprot] [--skip-orthosteric] [--skip-validation]
                                   [--include-metals] [--include-mutations]
                                   [--index-header INDEX_HEADER]
                                   [--ligands-header LIGANDS_HEADER] [--save-raw-data]

Process Boltz-2 prediction data

options:
  -h, --help            show this help message and exit
  --predictions-dir PREDICTIONS_DIR
                        Path to directory of cif files. This directory must contain cif files
                        with {accession}_{arbitrary_field}_{ligand_name}_model_{num}.cif
  --output OUTPUT       Output Excel file path
  --index-file INDEX_FILE
                        Path to index file with ligand information
  --index-sheet INDEX_SHEET
                        sites sheet name in index file
  --ligand-sheet LIGAND_SHEET
                        compounds sheet name in index file
  --site-cutoff SITE_CUTOFF
                        Distance cutoff for binding site
  --filter-cys-domains  Only include domains that contain the cysteine site in output
  --include-extended-fields
                        Include all verbose fields in output
  --skip-uniprot        Skip fetching data from UniProt API
  --skip-orthosteric    Skip fetching orthosteric site data from UniProt
  --skip-validation     Skip Phenix/PoseBusters validation
  --include-metals      Include metal ions in orthosteric sites analysis
  --include-mutations   Include mutagenesis sites in orthosteric sites analysis
  --index-header INDEX_HEADER
                        Header in the index sheet
  --ligands-header LIGANDS_HEADER
                        Header in the ligand sheet
  --save-raw-data       Skip Excel formatting and write full output to an Excel file
```
