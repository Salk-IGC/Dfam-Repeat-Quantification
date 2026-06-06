# Dfam Repeat Quantification in Subtelomeric Reads

A pipeline for quantifying transposable element and repeat motifs from the [Dfam database](https://dfam.org/) in subtelomeric long reads, developed at the Salk Institute IGC (Karlseder Lab).

## Overview

This workflow characterizes the repeat content of subtelomeric long reads by scanning them against the Dfam curated HMM library using `nhmmscan`. It then counts the number of distinct reads carrying each repeat motif and calculates per-motif coverage as a percentage of total reads. Overlap between sample sets is visualized as Venn diagrams in R.

The pipeline was developed using IMR90 fibroblast telomere-containing long reads from two cell lines:
- **IMR90SV** – SV40-immortalized IMR90 (ALT-negative)
- **IMR90ALT** – ALT-positive IMR90 cells

## Repository Contents

| File | Description |
|------|-------------|
| `dfam_analysis.Rmd` | End-to-end analysis notebook — database download, `nhmmscan` commands, AWK-based counting, and R Venn diagram generation |
| `unique_readsv2.awk` | AWK script to count distinct reads per motif and compute percentages from `nhmmscan` tabular output |

## Dependencies

**Bioinformatics tools (via Singularity):**
- [Dfam TE Tools container](https://hub.docker.com/r/dfam/tetools) — includes HMMER 3.4 (`nhmmscan`, `hmmpress`) and TRF
- [Dfam curated HMM library](https://dfam.org/releases/current/families/) (release 3.9, March 2025)

**R packages:**
- `VennDiagram`

## Workflow

### 1. Download and Prepare the Dfam HMM Database

```bash
# Download Dfam curated motifs (release 3.9)
wget https://dfam.org/releases/current/families/Dfam-curated_only-1.hmm.gz
gunzip Dfam-curated_only-1.hmm.gz

# Press the database for fast searching
hmmpress Dfam-curated_only-1.hmm
mkdir Dfam-curated-only-1-hmm
mv Dfam-curated* Dfam-curated-only-1-hmm/
```

### 2. Convert FASTQ to FASTA

`nhmmscan` requires FASTA input. Convert your telomere-containing reads with:

```bash
sed -n '1~4s/^@/>/p;2~4p' your_reads.fq > your_reads.fasta
```

### 3. Run nhmmscan (inside the Singularity container)

```bash
# Pull and enter the TE Tools container
singularity pull dfam-tetools-latest.sif docker://dfam/tetools:latest
singularity run dfam-tetools-latest.sif

# Set PATH inside the container
export PATH=/opt/hmmer/bin:/opt/trf:$PATH

# Scan reads against the Dfam HMM library
nhmmscan --noali --tblout output_tbl.out --cpu=56 \
    Dfam-curated-only-1-hmm/Dfam-curated_only-1.hmm \
    your_reads.fasta
```

### 4. Count Unique Reads per Motif

Use the provided AWK script to count how many **distinct reads** contain at least one hit for each motif, along with the percentage of total reads:

```bash
awk -f unique_readsv2.awk output_tbl.out > motif_unique_read_counts.txt
```

The script:
- Skips the 2-line header and the 10-line footer of `nhmmscan` tabular output
- Deduplicates `(motif, read)` pairs so each read is counted once per motif
- Reports distinct read counts and percentages across all motifs

**Example output:**
```
TANDEM MOTIF              DISTINCT READS PERCENT (%)
------------------------- --------------- ----------
HSATII                              1823      40.52
TAR1                                 945      21.01
...
------------------------- --------------- ----------
PROCESSED LINE RANGE: 3 through 9572158
TOTAL UNIQUE READS (in range): 4498
```

You can also get simple per-motif hit counts (without read-level deduplication) using:

```bash
# Total unique motifs detected
awk 'NR > 2 {seen[$1]++} END {print length(seen)}' output_tbl.out

# Count table
awk 'NR > 2 { count[$1]++ } END { for (v in count) print v, count[v] }' \
    output_tbl.out > motif_counts.txt
```

### 5. Visualize Overlap Between Samples (R)

```r
library(VennDiagram)

sampleA <- read.delim("sampleA_counts.txt", sep=" ", header=FALSE)
sampleB <- read.delim("sampleB_counts.txt", sep=" ", header=FALSE)

overlap <- list(sampleA$V1, sampleB$V1)
names(overlap) <- c("SampleA", "SampleB")

venn.diagram(
    x = overlap,
    category.names = c("SampleA", "SampleB"),
    filename = "sampleA_sampleB_venn.png",
    main = "Overlap of Dfam Motifs",
    imagetype = "png"
)
```

A helper function `simpleVennCSV()` (defined in `dfam_analysis.Rmd`) exports the Venn partition membership to a CSV for downstream analysis.

## Input / Output Summary

| Step | Input | Output |
|------|-------|--------|
| `nhmmscan` | `.fasta` reads + Dfam HMM db | `*_tbl.out` tabular hits file |
| `unique_readsv2.awk` | `*_tbl.out` | `*_motif_unique_read_counts.txt` |
| AWK one-liners | `*_tbl.out` | `*_counts.txt` (motif hit counts) |
| R / VennDiagram | `*_counts.txt` files | Venn PNG + overlap CSV |

## Notes

- The `nhmmscan` `--tblout` format includes a 2-line comment header and a trailing summary block (~10 lines); `unique_readsv2.awk` handles this automatically and exits with an error if the file is too small to parse safely.
- Adjust `--cpu` to match your available cores.
- The Singularity container approach ensures reproducibility across HPC environments without requiring local HMMER installation.

## Citation / Attribution

Dfam database: Storer J, Hubley R, Rosen J, Wheeler TJ, Smit AF. *The Dfam community resource of transposable element families, sequence models, and genome annotations.* Mobile DNA. 2021.
