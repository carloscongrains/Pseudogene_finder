# Pseudogene_finder

This script is designed to identify potential pseudogenes within a query sequence dataset by comparing them against a reference sequence. It focuses on detecting common pseudogene characteristics, specifically premature stop codons and frameshift mutations, which disrupt the open reading frame (ORF) and lead to non-functional protein products.

How it Works
The script takes a query FASTA file (containing sequences you want to analyze) and a reference FASTA file (containing known functional genes or proteins). It is recommended to use a complete reference sequence (including start codon and stop codon). The program aligns the query sequence with the reference using MAFFT algorithm (Katoh et al. 2002) and scans for:

Premature Stop Codons: Stop codons appearing before the expected end of the coding sequence, indicating a truncated protein.

Frameshift Mutations: Insertions or deletions not multiple of 3 that leads to a completely different downstream amino acid sequence and often premature stop codons.

Getting Started
Prerequisites

Python 3.x
Biopython 
Mafft 7.x (It should be in the path)

Installation using conda

This is the recommended way to install the program. It requires conda to be already installed. Then follow the steps:

#Create a environment
conda create --name pseudogene_finder
#Open the environment
conda activate pseudogene_finder
#Install biopython 
conda install conda-forge::biopython
#Install MAFFT
conda install bioconda::mafft

#Activate the environment any time the user runs the script

---
Usage
To run the script, use the following command structure:

```
python get_stop_codon_frameshifts.py --input_fasta query.fasta --input_reference reference.fasta --genetic_code_table 5 --output_dir output_dir
```

Command-Line Arguments

--input_fasta <path/to/query.fasta> (Required)

Specifies the path to the input multi FASTA file containing the query nucleotide sequences to be analyzed for pseudogene features. This file can include more than one sequence

--input_reference <path/to/reference.fasta> (Required)

Specifies the path to the input FASTA file containing the reference in nucleotide sequence (only one reference is allowed). The reference must start in the 1st codon position and its length must be a multiple of 3. These are used for comparison and determining expected coding regions.

--genetic_code_table <table_number> (Required)

Specifies the genetic code table to use for translation. Genetic code tables available in this script are: 1 to 6, and 9 to 16 (for more information check this website: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=tgencodes#SG). This is crucial for correctly identifying start and stop codons.

Example: 5 refers to the Invertebrate Mitochondrial genetic code.

--output_dir <path/to/output_directory> (Required)

Specifies the path to the directory where all output files will be saved. The script will create this directory if it does not already exist.

---
Example
Let's say you have a file my_genes.fasta with sequences to check and known_proteins.fasta as your reference. You want to use the Standard Genetic Code (table 1) and save results to a folder named pseudogene_results.

```
python get_stop_codon_frameshifts.py --input_fasta my_genes.fasta --input_reference known_proteins.fasta --genetic_code_table 1 --output_dir pseudogene_results
```
---
Output
The script will generate various output files within the specified --output_dir, including

A summary report detailing identified pseudogenes in a tab-delimited file named: output_results.txt and located into output_dir location.

Fasta files of every combination of query and reference and their corresponding alignments.

Reference
