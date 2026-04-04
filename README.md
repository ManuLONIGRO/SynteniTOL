# SynteniToL

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Docker Image](https://img.shields.io/badge/Docker-manulonigro%2Fsyntenitol:latest-blue)](https://hub.docker.com/r/manulonigro/syntenitol)

A Nextflow DSL 2 pipeline for synteny analysis of protein sequences against genomic annotations, using HMMER for domain detection and generating visualizations for [iTOL](https://itol.embl.de/) (Interactive Tree of Life).

---

## Table of Contents

- [About](#about)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage](#usage)
- [Pipeline Overview](#pipeline-overview)
- [License](#license)
- [Contact](#contact)

---

## About

SynteniToL enables researchers to perform comparative micro-synteny analysis across multiple prokaryotic genomes by detecting homologous proteins (using COG profiles) and generating publication-ready visualizations for iTOL.

**Key features:**

- Automated retrieval of genome annotations from NCBI
- HMMER-based domain detection using COG (Clusters of Orthologous Groups) profiles
- Synteny context extraction and visualization
- iTOL-compatible output files for interactive phylogenetic trees and synteny plots
- Taxonomy-aware visualizations
- Fully containerized with Docker for reproducible execution

This pipeline is particularly useful for microbiologists studying gene cluster evolution, horizontal gene transfer, and genomic rearrangements across related species.

---

## Installation

### Prerequisites

- [Nextflow](https://www.nextflow.io/) (version 21.0 or higher)
- [Docker](https://www.docker.com/)

---

## Quick Start

**Input FASTA file:**

The input FASTA file should contain protein sequences. Each sequence header should include an accession number that can be resolved by NCBI (e.g., from BLASTp results):

```fasta file
>WP_145111153.1:6-544 formylmethanofuran...
MVLSPADKTNVKAAVKGEG...
>WP_200486823.1:5-533 formylmethanofuran...
MLSGIV...
```

**Run the pipeline** with minimal required parameters:

```bash
nextflow run main.nf \
    --inputFASTA input_proteins.faa \
    --cogs COG1152,COG1795 \
```

**Output:**

```
results/run_timestamp/
├── class.itol.txt              # iTOL taxonomic info file: Classes
├── phylum.itol.txt             # iTOL taxonomic info file: Phylums
├── domain.itol.txt             # iTOL taxonomic info file: Domains
├── tree_timestamp.nwk          # Phylogenetic tree
├── itol_binary.txt             # iTOL dataset: Mapped genes
├── itol_synteny_oriented.txt   # iTOL dataset: Synteny
├── formatted_headers.fasta     # header format: protein_id|specie_name|strain|isolate
└── run_command.txt             # Reproducibility record
```

---

## Usage

```
Do a blastp in NCBI (https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome)
Get the COG id of the query protein

Go to your terminal, move to SynteniToL dir:
nextflow run main.nf [OPTIONS]
```

### Required Parameters

| Parameter          | Description                                                               |
|--------------------|---------------------------------------------------------------------------|
| `--inputFASTA`     | Input FASTA file containing protein sequences (e.g., from BLASTp results) |
| `--cogs`           | Comma-separated list of COG identifiers (e.g., `COG1152,COG1795`)         |
**IMPORTANT: the first cog in --cogs is the COG of the protein homologs in the FASTA file**

### Optional Parameters

| Parameter          | Description                                                               |
|--------------------|---------------------------------------------------------------------------|
| `--ncbi_api_key`   | NCBI API key to increase rate limits (recommended for large queries)      |
| `--color_by_group` | Group COGs for coloring in iTOL (format: `COG1229-COG1029,COG2218`)       |
| `--outdir`         | Output directory for results                                              |
| `--help`           | Display help message                                                      |

### Examples

**Basic usage with NCBI API key:**

```bash
nextflow run main.nf \
    -profile docker \
    --inputFASTA proteins.faa \
    --cogs COG1152,COG1795,COG2218 \
    --ncbi_api_key your_ncbi_api_key
```

**With group coloring (COG1229 and COG1029 share the same color):**

```bash
nextflow run main.nf \
    -profile docker \
    --inputFASTA proteins.faa \
    --cogs COG1152,COG1229,COG1029,COG2218 \
    --color_by_group "COG1229-COG1029,COG2218,COG1152"
```

### Getting an NCBI API Key

1. Sign in to your [NCBI account](https://www.ncbi.nlm.nih.gov/account/)
2. Go to [API Keys Settings](https://www.ncbi.nlm.nih.gov/account/settings/)
3. Generate a new API key

> **Note:** Without an API key, NCBI limits requests to 3 per second. With an API key, the limit increases to 10 per second.

---

## Pipeline Overview

```
┌─────────────┐
│  Input      │
│  FASTA      │
└──────┬──────┘
       │
           ▼
┌─────────────────────────────────────────────────────────────────┐
│  1. Fetch genome metadata from NCBI (efetch_to_tsv)             │
└─────────────────────────────────────────────────────────────────┘
       │
           ▼
┌─────────────────────────────────────────────────────────────────┐
│  2. Download genome annotations from NCBI Assembly Database     │
│     - Split into batches of 50 for efficient downloading        │
│     - Handle missing assemblies gracefully                      │
└─────────────────────────────────────────────────────────────────┘
       │
           ▼
┌─────────────────────────────────────────────────────────────────┐
│  3. Download COG profiles from NCBI                             │
└─────────────────────────────────────────────────────────────────┘
       │
           ▼
┌─────────────────────────────────────────────────────────────────┐
│  4. Build HMM profiles from COG multiple sequence alignments    │
└─────────────────────────────────────────────────────────────────┘
       │
           ▼
┌─────────────────────────────────────────────────────────────────┐
│  5. HMMER search against proteomes                              │
│     - Scan all genomes for COG domain presence                  │
│     - Extract locus tags and protein mappings                   │
└─────────────────────────────────────────────────────────────────┘
       │
           ▼
┌─────────────────────────────────────────────────────────────────┐
│  6. Build results dataframe and extract synteny contexts        │
│     - First COG in list = reference for synteny plots           │
│     - Calculate gene neighborhoods and conservation             │
└─────────────────────────────────────────────────────────────────┘
       │
           ▼
┌─────────────────────────────────────────────────────────────────┐
│  7. Generate iTOL visualization files                           │
│     - Synteny plots, presence/absence matrices, taxonomy trees  │
└─────────────────────────────────────────────────────────────────┘
```

**Pipeline steps in detail:**

1. **NCBI Data Fetching**: Parses input FASTA and retrieves genome metadata from NCBI using E-utilities
2. **Genome Download**: Downloads genome annotations (GBFF) from NCBI Assembly Database
3. **COG Profile Download**: Fetches COG multiple sequence alignments from NCBI
4. **HMM Profile Building**: Builds HMM profiles from COG alignments using `hmmbuild`
5. **Domain Detection**: Scans all proteomes with HMMER using `hmmsearch`
6. **Gene Mapping**: Maps protein IDs to locus tags and gene positions
7. **Synteny Analysis**: Extracts syntenic contexts around the first COG of --cogs
8. **Visualization**: Generates iTOL-compatible files for interactive visualization

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Contact

**Author:** ManuLONIGRO
**Email:** lonigromanuel@gmail.com

For questions, issues, or feature requests, please open an issue on the project repository.

