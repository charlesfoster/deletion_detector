# deletion_detector
A simple program that uses parallel processing to detect deletions in fasta files relative to a reference genome. 

## Installation
Firstly, clone this repository:

```
git clone https://github.com/charlesfoster/deletion_detector.git

cd deletion_detector
```

Create a new conda environment:

```
# install mamba if not already installed
conda install -c conda-forge mamba
mamba env create -f environment.yml
```

Install `deletion_detector`:

```
conda activate deletion_detector
pip install .
```

## Usage
Each time you wish to run the program, make sure to activate the conda environment first:

```
conda activate deletion_detector
```

Check the usage options:

```
deletion_detector -h
```

Output:

```
usage: deletion_detector [-h] [-c COORDINATES] [-d DELETION_OF_INTEREST] [-f] [-p] [-r REFERENCE] [-t THREADS] [-o OUTFILE] [--quick]
                         [--version]
                         [fasta ...]

positional arguments:
  fasta                 Fasta file containing sequences to check (default: None)

options:
  -h, --help            show this help message and exit
  -c COORDINATES, --coordinates COORDINATES
                        Survey subset of genome, e.g. an open reading frame. Must be specified in the form of start:end (1-based
                        coordinates; inclusive). If specified, stats will be given for protein truncation etc. Currently assumes standard
                        translation table. (default: False)
  -d DELETION_OF_INTEREST, --deletion_of_interest DELETION_OF_INTEREST
                        Specify a deletion of interest in the form of 'start:end' (1-based coordinates; inclusive). If specified, outfile
                        can be easily filtered to find 'target_deletion'. (default: None)
  -f, --force           Force overwrite outfile (default: False)
  -p, --parse_gisaid    If analysing fasta files from GISAID, parse out the date of collection and country (default: False)
  -r REFERENCE, --reference REFERENCE
                        Reference genome in fasta format (default: /home/cfos/Programs/deletion_detector/MN908947.3.fasta)
  -t THREADS, --threads THREADS
                        Number of threads to use in parallel processing. Defaults to all available threads. (default: 20)
  -o OUTFILE, --outfile OUTFILE
                        Name of the outfile to store results (default: deletion_results.tsv)
  --quick               Run in 'quick' mode: faster than default mode, but less rich output (default: False)
  --version             show program's version number and exit

Analysis quits without running if outfile already exists and --force not specified. Limitations: (a) gap detection is always subject to how
easy/difficult a region is to align to the reference genome, (b) if there are large deletions close to the beginning/end of a sequence, soft
clipping of the alignment might introduce false missing data (padded as 'N's prior to deletion analysis).
```

## What does the program do?
Let's say you have a fasta-format file with one or more sequences that you wish to compare to a reference sequence to find any deletions. That's when this program comes in handy. All input query sequences are compared to the reference sequence, and the output file (TSV format) will list all deletions in query sequences relative to the reference and metrics about the number of Ns and gaps. 

Note: the input sequences do NOT need to be aligned: `deletion_detector` assumes they aren't aligned, and aligns them as part of the analysis. At a future stage I might implement an option to skip alignment (for pre-aligned sequences).

## Benchmark
All sequences within the provided file are compared to the specified reference sequence using parallel computation (tested on Ubuntu 22.04 and Mac OS), running in chunks so as to (hopefully) not exhaust all RAM. This method allows very large input files to be processed, but will vary from computer to computer based on your setup. If you receive out of memory errors, you probably need to split your input fasta into smaller files - judge based on your available RAM. 

There are two 'modes' to the program:
- "Rich" mode: outfile has the deletion results and various QC metrics (see ["What do you get?"](#What-do-you-get) below)
- "Quick" mode: outfile has the deletion results, but no extra QC metrics

As a guide, with my 10-core (20-thread) processor (Intel(R) Core(TM) i9-10900X CPU @ 3.70GHz) and 64 GB RAM, I can process a 10.7 GB file with 351,995 sequences in:
- "Rich" mode: ~582 seconds (just over nine and a half minutes), equating to an average of 867.53 iterations (sequences) per second
- "Quick" mode: ~399 seconds (just over six and a half minutes), equating to an average of 1289.65 iterations (sequences) per second

## What do you need to run the program?
At minimum:
* <sequences>.fasta: a fasta-format file with at least one sequence (mandatory)
* <reference>.fasta: a fasta-format file with one reference sequence for comparison (assumes SARS-CoV-2 reference genome by default)

## What do you get?
Assuming you are running in "rich" mode (default), the output file (in .tsv format) contains statistics about the deletions in query sequences, and basic QC metrics. The main output columns, in order, contain:

* The full query sequence ID
* The stretches of gap(s) identified in a query sequence, in the form of '(start coordinate, end coordinate, number of consecutive gap sites)'
* Number of inferred deletions (+1 for each "gap stretch" identified)
* A summary of the 'deletion status' of the query sequence: 'deletion', 'deletion_plus' (multiple deletions), or 'no_deletion'.
* Percentage of gaps (-) in the sequence
* Percentage of Ns in the sequence
* Overall QC based on the percentage of Ns: 'NO_N', 'LOW_N', 'MILD_N', 'MEDIUM_N', 'HIGH_N', 'FAIL'

If running in "quick" mode, you will have no QC criteria in the output. So: you will have your results in ~2/3 of the time, but downstream filtering of the outfile (e.g., by coverage of an open reading frame) won't be possible.

Additionally, coordinates for a deletion of interest can be specified, in which case the 'deletion status' has two additional possibilities: 'target_deletion' or 'target_deletion_plus'. The outfile can then be easily filtered to find the 'target_deletion'.

If coordinates of an open reading frame are specified, the QC metrics above are only for that open reading frame (note the overall sequence). Additional columns will also be provided:

* The inferred amino acid stop position
* The percentage of the full-length peptide sequence (<100% if early stop codon identified)
* Whether the peptide sequence is truncated (Boolean: true/false)

The latter stats only make sense if the coordinates are for an open reading frame, and assume the 'Standard' translation table.

## SARS-CoV-2 specific option
I wrote this program with SARS-CoV-2 in mind since I was interested in querying consensus genomes to find the prevalence of deletions of interest in global data sets. Because of the SARS-CoV-2 focus, I've provided an optional argument to parse GISAID-format sequence headers in input sequences. If this option is used, the output TSV file will have additional columns with useful information extracted for each sample: the date of collection (DOC), country, and a shortened sequence ID. The `--parse_gisaid` option assumes that sequences have been downloaded from GISAID, either from the Search section or a whole-clade download, and that the sequence headers have not been modified since the download. 

## Citations
If you use this program and find it useful, I'd appreciate some kind of attribution. For example, I have uploaded the initial release to Zenodo, which can be cited like so:

Foster, C.S.P. (2022). deletion_detector: a simple and quick way to detect deletions relative to a reference sequence (v0.1.0). Zenodo. https://doi.org/10.5281/zenodo.7049120

I will also upload future releases to Zenodo.

Deletion Detector depends on:

* `gofasta`: B. Jackson, gofasta: command-line utilities for genomic epidemiology research. Bioinformatics. 38, 4033–4035 (2022).
* `minimap2`:H. Li, Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics. 34, 3094–3100 (2018).
* `tqdm`: C. da Costa-Luis, et al., tqdm: A fast, Extensible Progress Bar for Python and CLI (2022), doi:10.5281/zenodo.7046742.
* `biopython`: P. J. A. Cock, et al., Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics. 25, 1422–1423 (2009).


