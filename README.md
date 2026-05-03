# CRISPR-Cas9 sgRNA Designer

## Overview
This is a computational biology project that acts as an automated tool for designing CRISPR-Cas9 Single Guide RNAs (sgRNAs). Utilizing **Python** and **Biopython**, it accesses the global DNA database to download target human genes and scans them for viable genome-editing targets.

## What it does:
1. **API Integration**: Connects to the NCBI database via `Bio.Entrez` to programmatically download target genes (Default is the human `HBB` gene responsible for Sickle Cell Anemia).
2. **PAM Recognition**: Uses regular expressions to scan the genomic sequence for the SpCas9 PAM motif (`NGG`).
3. **Biological Filtration**: Extracts 20-bp upstream guide RNA candidates and evaluates their viability by calculating GC-content via `Bio.SeqUtils`. Filtering outputs to only display optimal targets (40-60% GC).

## Requirements
* Python 3.x
* Biopython (`pip install biopython`)

## How to use:
Simply execute the python script in your terminal:
`python crispr_designer.py`
