# Griffin-Pipeline
> Nanopore sequencing data analysis pipeline

This repository provides a reproducible and modular pipeline for **Oxford Nanopore Technologies (ONT) sequencing data preprocessing and analysis**, from basecalling, quality assessment, generating high-quality metagenome-assembled genomes (MAGs), constructing a tailored Kraken2 reference database, and taxonomic classification.

The goal is to promote **research reproducibility** by documenting each step clearly and using widely available open-source tools.

---

## Repository structure
 * Each module contains a detailed README file. Click on the README for more information.

### 1. Preprocessing

This module prepares raw Oxford Nanopore sequencing data for downstream analyses by converting signal-level data into high-quality, sample-resolved reads. It emphasizes rigorous quality control and standardized reporting to ensure data integrity and reproducibility.

  * [README](./Preprocessing.md)
  * The preprocessing module includes the following core components:
    1. **Basecalling and demultiplexing**: Converts raw electrical signals into nucleotide sequences and assigns reads to their corresponding samples.
    2. **Quality control**: Assesses read quality and summarizes key metrics.
    3. **Quality reporting**: Generates an integrated QC report to facilitate data assessment and visualization.
    4. **Read filtering and trimming** (optional, data-dependent): Removes ultra-long reads (>100,000 bp) and trims long reads to 14,999 bp. The filtered reads are the input to the module below.

### 2. Custom Database Generation

This module provides a comprehensive workflow for constructing a custom reference database based on metagenome-assembled genomes (MAGs) generated from long-read sequencing data. The resulting database is optimized for taxonomic classification of Oxford Nanopore reads.

  * [README](./HMAcustomDB.md)
  * The workflow includes:
    1. **Host read removal**: Removes reads that align to the host (mouse) or human genome. The host-free reads are the input to the assembly and to the Taxonomic Classification module.
    2. **Assembly**: Generates contigs from host-free long reads (five independent runs per sample).
    3. **Polishing**: Corrects INDELs and other sequencing errors in the draft assembly.
    4. **Binning**: Groups contigs into MAGs based on sequence composition and coverage patterns.
    5. **Quality assessment**: Evaluates MAG completeness and contamination to retain high-quality genomes.
    6. **Pooling**: Keeps species recovered reproducibly across the independent runs.
    7. **Taxonomic classification**: Assigns GTDB taxonomy to the pooled MAGs.
    8. **Dereplication**: Collapses MAGs shared across samples into one MAG per species.
    9. **Database construction**: Builds a custom Kraken2 database from the MAGs themselves, with a GTDB-based taxonomy.

### 3. Taxonomic Classification

This module performs read-level taxonomic assignment and abundance profiling using a custom reference database derived from curated MAGs. Ambiguous read assignments can optionally be removed with a _k_-mer-based read filter before abundance estimation.

  * [README](./TaxAnnotation.md)
  * The module includes the following core components:
    1. **Input**: Host-free reads from the Custom Database Generation module
    2. **Taxonomic classification**: Assigns reads to taxonomic lineages using Kraken2 with a custom reference database
    3. **_k_-mer-based read filtering** (optional): Removes ambiguous classifications based on unassigned _k_-mer proportions and read-level _k_-mer diversity
    4. **Quantification and abundance estimation**: Calculates species-level relative abundance with **BASEN v1.0.0** ([Snyder-Institute/BASEN](https://github.com/Snyder-Institute/BASEN))

---

## Key Features
  * Modular design for flexible integration into diverse projects
  * Fully documented steps for transparent, reproducible workflows
  * Support for open-source, widely adopted bioinformatics tools
  * Scalable to small and large sequencing datasets

## Getting Started
 * Instructions for installation, dependencies, and example commands are provided in individual module READMEs.
 * A short tutorial is available on the [wiki](https://github.com/Snyder-Institute/Griffin-Pipeline/wiki).

