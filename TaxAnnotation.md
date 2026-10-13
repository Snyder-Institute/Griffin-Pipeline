# Taxonomic Classification
> Read-level taxonomic assignment and abundance profiling using a custom Kraken2 database

## Rationale
Accurate taxonomic profiling of the long-read metagenomic data requires both high-quality reads and context-aware reference databases. Long-read sequencing introduces higher base error rates, and naïve classification against large, generic databases can lead to ambiguous assignments and inflated false positives.

This module ensures robust and reproducible estimation of microbial composition and engraftment dynamics:
  * Reduces misclassification caused by sequencing errors
  * Eliminates spurious taxa absent from the biological system
  * Improves abundance estimation accuracy
  * Enhances interpretability in engraftment and colonization studies

## Workflows

1. **Input**: Host-free reads ([Custom Database Generation](./HMAcustomDB.md) Step 01)
2. **Taxonomic classification**: Assign reads using Kraken2 with a custom database
3. **_k_-mer-based read filtering** (optional): Remove ambiguous read assignments
4. **Quantification**: Estimate species-level relative abundance with **BASEN v1.0.0**

### Environment setup
Set the following environment variables to define input data, outputs, and database locations.

```bash
# Define base folders
export HOSTFREE_FOLDER=/path/to/hostfree        # Host-free reads (Custom Database Generation Step 01)
export CLASSIFICATION_FOLDER=/path/to/kraken2
export CUSTOMDB_FOLDER=/path/to/custom_db
export PROFILE_FOLDER=/path/to/basen
```

### Dependencies

Install or load the following tools prior to execution:

| Tool | Version | Purpose | Reference |
|------|---------|----------|------------|
| [Kraken2](https://github.com/DerrickWood/kraken2) | 2.1.6 | High-throughput taxonomic classification of sequencing reads | Wood et al. |
| [BASEN](https://github.com/Snyder-Institute/BASEN) | **v1.0.0** | Base-level Abundance estimation with Species-assigned Evidence using Nanopore | _Not published_ |

> **BASEN version.** All commands and results here use **BASEN v1.0.0** (GitHub commit `fad69ed`), in which the effective genome length is Gₛ = L − k + 1 (L = genome length in bp, k = 35). Check the installed version before running: `Rscript -e 'packageVersion("basen")'` should return `1.0.0`.

Kraken2 and BASEN have no random-seed option; Kraken2 classification is deterministic.

### Example directory structure
```bash
├── hostfree/               # Host-free reads (input)
├── kraken2/                # Kraken2 outputs: <sample>.kraken, <sample>.report
├── custom_db/              # Custom Kraken2 database and genome_stats.txt
└── basen/                  # Relative abundance profiles
```

### Schematic overview

![RelativeAbundance](./.github/profiles.png)

### Pipeline steps
#### Step 01. Input reads
Reads are optionally filtered by length in the [Preprocessing](./Preprocessing.md) module (Step 09), and host reads are removed in the [Custom Database Generation](./HMAcustomDB.md) module (Step 01). Classifying reads that still contain host sequences inflates the unclassified fraction and lets host reads compete for _k_-mer assignments.

#### Step 02. Taxonomic classification with Kraken2

Host-free reads are classified against the custom Kraken2 database generated from donor-derived MAGs (see [Custom Database Generation](./HMAcustomDB.md)). Kraken2's own `--confidence` option is left at its default (0); ambiguous assignments can be removed in Step 03 (optional) instead.

```bash
k2 classify \
  --db $CUSTOMDB_FOLDER/myDB \
  --threads 8 \
  --report $CLASSIFICATION_FOLDER/F0.report \
  --output $CLASSIFICATION_FOLDER/F0.kraken \
  $HOSTFREE_FOLDER/F0.fastq.gz
```

#### Step 03. _k_-mer-based read filtering (optional)

We evaluated taxonomic classification accuracy using a defined positive control dataset, OMM12 (oligo-mouse microbiota), through comparative benchmarking of Kraken2 and Minimap2 coupled with CoverM. Each method was tested against two custom reference databases: (i) OMM12 alone and (ii) OMM12 augmented with decoy genomes (e.g., the human reference genome).

This benchmarking revealed two key patterns. First, Kraken2 produced more reliable taxonomic assignments than Minimap2 when validated against BLAST results. This discrepancy is likely attributable to Minimap2’s reliance on minimizers as alignment anchors; loss of these anchors can disrupt chaining and result in missed alignments, particularly in repetitive genomic regions. Second, we conducted _k_-mer–level profiling using Kraken2’s default _k_ value (_k_ = 35) and observed that the number of _k_-mers assigned to decoy genomes increased when the proportion of unclassified reads was high. This suggests that unclassified reads may compete with dominant species for _k_-mer assignments.

To quantify this effect, we applied the Shannon diversity index at the read level, treating each read as a discrete distribution of _k_-mers across assignment categories (i.e., assigned to a specific reference genome versus unassigned). This formulation enables direct comparison of the relative abundance of assigned and unassigned _k_-mers within individual reads.

A classified read is excluded if **either** of the following holds:
  1. more than 75% of its _k_-mers are unassigned (unassigned proportion > 0.75); **or**
  2. its Shannon diversity index exceeds the cut-off for its sample type: 1.0 for the high-depth donor samples (one per donor, `<donor>_NA_0_1`), and 0.75 for all other samples (mouse samples, donor replicates `<donor>_NA_0_2`, and controls).

Unclassified reads are removed as well. BASEN computes the read-level _k_-mer profile (`profile_kmer_composition()`) and returns the IDs of the reads that pass both rules (`filter_low_info_reads()`); these IDs are passed to Step 04.

```R
library(basen) # BASEN v1.0.0, https://github.com/Snyder-Institute/BASEN
stopifnot(packageVersion("basen") == "1.0.0")

kmer_profile <- profile_kmer_composition("/path/to/kraken2/F0.kraken", k_mer_length = 35)
max_shannon  <- 1.0 # high-depth donor sample (<donor>_NA_0_1); 0.75 for all other samples
hq_reads <- filter_low_info_reads(kmer_profile,
                                  max_unassigned_perc   = 0.75, # the BASEN default is 0.5
                                  max_shannon_diversity = max_shannon)
```

#### Step 04. Quantification and abundance estimation with BASEN v1.0.0
Bracken is based on a Bayesian framework and requires a fixed read length as a prior, which makes it suboptimal for long-read sequencing data. To address this limitation, we developed a _k_-mer–based approach for estimating relative abundance, termed BASEN (Base-level Abundance estimation with Species-assigned Evidence using Nanopore).

BASEN estimates species-level relative abundance from long-read metagenomic data by restricting classification to species-level assignments and counting only species-assigned _k_-mers as base-level evidence. For each species s, the evidence Bₛ (species-assigned _k_-mers summed over the classified reads, or over the reads that passed Step 03 if it is applied) is divided by the effective genome length Gₛ = L − k + 1, giving a coverage proxy Cₛ = Bₛ / Gₛ. Within each sample, the coverage proxies are scaled to sum to 1: RAₛ = Cₛ / ΣᵢCᵢ.

The code below computes relative abundance with BASEN v1.0.0 from the `.report` and `.kraken` files of Step 02 and `genome_stats.txt` from the database build; `read_filter` keeps only the reads that passed Step 03; without it, all classified reads are used. The files must be named `<sample>.report` and `<sample>.kraken`.

```R
resultsL <- kraken_relative_abundance(
        report_dir = "/path/to/kraken2",                        # <sample>.report (Step 02)
        kraken_dir = "/path/to/kraken2",                        # <sample>.kraken (Step 02)
        k_mer_length = 35,
        genome_stats_file = "/path/to/custom_db/genome_stats.txt",
        sample_ids = "F0",
        read_filter = hq_reads,                                 # optional: reads that passed Step 03; omit to use all classified reads
        output_folder = "/path/to/basen"
)
```

---
