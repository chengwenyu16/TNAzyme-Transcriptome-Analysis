TNAzyme Transcriptome Analysis

This repository contains the Python scripts and computational pipeline used for transcriptome scanning and DOE mapping analysis in the development of artificial ribozymes (TNAzymes).

*Repository Structure*

01_transcriptome_scan.py: Script for global transcriptome scanning. It analyzes human mRNA transcripts (GENCODE) to identify targetable motifs and calculates overall coverage and cleavage site frequencies.

02_doe_inhibitor_analysis.py: Script for mapping DOE inhibitor targets against the MANE Select transcript database. It outputs a comprehensive Excel summary report detailing cleavable genes, motif counts, and mapping rates.

doe.csv: The input dataset containing target genes and their predicted inhibitor scores (threshold used in analysis: >= 0.61).

*Data Availability*

Due to GitHub's file size limitations (100 MB maximum), the large customized biological databases required to run these scripts are hosted on Zenodo.

You must download the following large data files and place them in a data/ directory relative to the scripts:

gencode.v50.transcripts.fa.gz (~175 MB) - Required for 01_transcriptome_scan.py

MANE.GRCh38.v1.5.select_ensembl_rna.txt (~74 MB) - Required for 02_doe_inhibitor_analysis.py

MANE.GRCh38.v1.5.summary.txt - Required for 02_doe_inhibitor_analysis.py

Download the complete dataset here: []

*How to Run*

1. Prerequisites

Ensure you have Python 3 installed along with the required libraries. You can install the dependencies using pip:

pip install biopython openpyxl


2. Directory Setup

The scripts expect a specific directory structure. Make sure your working directory looks like this:

├── data/                                 # Create this folder and put downloaded Zenodo files here
│   ├── gencode.v50.transcripts.fa.gz
│   ├── doe.csv                           
│   ├── MANE.GRCh38.v1.5.select_ensembl_rna.txt
│   └── MANE.GRCh38.v1.5.summary.txt
├── scripts/                              # Put the python scripts here
│   ├── 01_transcriptome_scan.py
│   └── 02_doe_inhibitor_analysis.py
└── results/                              # This folder will be generated automatically


3. Execution

Run the scripts from your terminal:

To run the transcriptome motif scanner:

python scripts/01_transcriptome_scan.py


Output: A text summary (transcriptome_scan_summary.txt) will be printed to the console and saved in the results/ folder.

To run the DOE inhibitor analysis:

python scripts/02_doe_inhibitor_analysis.py


Output: A comprehensive Excel workbook (Tz2_DOE_inhibitor_MANE_Select_8motifs.xlsx) will be saved in the results/ folder.
