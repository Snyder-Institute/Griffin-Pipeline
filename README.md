# Griffin-Pipeline
> Nanopore sequencing data analysis pipeline

This repository provides a reproducible and modular pipeline for **Oxford Nanopore Technologies (ONT) sequencing data preprocessing and analysis**, from basecalling, quality assessment, generating high-quality metagenome-assembled genomes (MAGs), constructing a tailored Kraken2 reference database, and taxonomy classification. 

The goal is to promote **research reproducibility and reproductivity** by documenting each step clearly and using widely available open-source tools.

---

## Repository structure
 * Each module contains a detailed README file. Click on the README for more information.

### 1. Preprocessing

This module prepares raw Oxford Nanopore sequencing data for downstream analyses by converting signal-level data into high-quality, sample-resolved reads. It emphasizes rigorous quality control and standardized reporting to ensure data integrity and reproducibility.

  * [README](./Preprocessing.md)
  * The preprocessing module includes the following core components:
    1. **Basecalling and demultiplexing**: Converts raw electrical signals into nucleotide sequences and assigns reads to their corresponding samples.
    2. **Quality control**: Assesses read quality, summarizes key metrics, and performs length- and quality-based filtering.
    3. **Read filtering and trimming**: Removes ultra-long or low-quality reads and trims suboptimal regions to improve downstream analyses.
    4. **Quality reporting**: Generates an integrated QC report to facilitate data assessment and visualization.
        
### 2. Custom Database Generation

This module provides a comprehensive workflow for constructing a custom reference database based on metagenome-assembled genomes (MAGs) generated from long-read sequencing data. The resulting database is optimized for taxonomic classification of Oxford Nanopore reads.

  * [README](./HMAcustomDB.md)
  * This module provides a complete workflow for constructing a custom reference database using metagenome-assembled genomes (MAGs). The workflow includes:
    1. **Assembly**: Generates contigs from long-read sequencing data.
    2. **Polishing**: Corrects INDELs and other sequencing errors in the draft assembly.
    3. **Binning**: Groups contigs into MAGs based on sequence composition and coverage patterns.
    4. **Quality assessment**: Evaluates MAG completeness and contamination to retain high-quality genomes.
    5. **Taxonomic classification**: Assigns taxonomic identities to the curated MAGs.
    6. **Genome retrieval**: Downloads or extracts genome sequences corresponding to the classified MAGs.
    7. **Database construction**: Builds a custom reference database using the validated and taxonomically annotated MAGs.

### 3. Taxonomic Classification

This module performs read-level taxonomic assignment and abundance profiling using a custom reference database derived from curated MAGs. The workflow emphasizes stringent quality filtering and confidence-based classification to ensure accurate downstream interpretation.

  * [README](./TaxAnnotation.md)
  * The module includes the following core components:
    1. **Preprocessing**: Removes ultra-long reads (>100,000 bp) and trims long reads (>15,000 bp)
    2. **Taxonomic classification**: Assigns reads to taxonomic lineages using Kraken2 with a custom reference database
    3. **Confidence-based read filtering**: Removes low-confidence classifications based on assigned _k_-mer proportions
    4. **Quantification and abundance estimation**: Calculates taxon-level coverage and relative abundance

---

## Key Features
  * Modular design for flexible integration into diverse projects
  * Fully documented steps for transparent, reproducible workflows
  * Support for open-source, widely adopted bioinformatics tools
  * Scalable to small and large sequencing datasets

## Getting Started
 * Instructions for installation, dependencies, and example commands are provided in individual module READMEs. 
 * A short tutorial is available on the [wiki](https://github.com/Snyder-Institute/Griffin-Pipeline/wiki).

