# Griffin-Pipeline
> Nanopore sequencing data analysis pipeline

This repository provides a reproducible and modular pipeline for **Oxford Nanopore Technologies (ONT) sequencing data preprocessing and analysis**, from basecalling, quality assessment, and taxonomy classification. The goal is to promote **research reproducibility and reproductivity** by documenting each step clearly and using widely available open-source tools.

---
## Rationale

MinKNOW from Oxford Nanopore frequently crashes and does not provide the most up-to-date _Dorado_ basecaller. Therefore, establishing an external computational environment is recommended for scalability and reproducibility. The Dorado version reported in the Web summary is the _Dorado Basecall Server_ version, which is misleading. In addition, the latest _Dorado_ release (v1.X.X), following the major update, provides an expanded quality score range based on improved signal-to-sequence predictors.

## Preprocessing

The pipeline performs the following tasks:

1. **Basecalling** of raw POD5 files using Dorado.  
2. **Merging** individual BAM files into a single dataset.  
3. **Demultiplexing** by barcode.  
4. **Generating sequencing summaries** for merged and demultiplexed BAM files.  
5. **Converting BAM to FASTQ** format.  
6. **Computing quality statistics** using NanoPlot and FastQC.  
7. **Aggregating reports** with MultiQC.

### Environment setup

Before running the pipeline, set the following environment variables to define input locations, output directories, and software paths.

```bash
# Define tool locations
export DORADO=/path/to/dorado           # Path to Dorado executable directory

# Define base folders
export POD5_FOLDER=/path/to/pod5_files  # Directory containing input POD5 files
export BAM_FOLDER=/path/to/bam_output   # Directory for basecalled BAM files
export DEMUX_FOLDER=/path/to/demuxed    # Directory for demultiplexed BAMs
export SUMMARY_FOLDER=/path/to/summary  # Directory for Dorado summary outputs
export FASTQ_FOLDER=/path/to/fastq      # Directory for FASTQ output files
export STAT_FOLDER=/path/to/nanoplot    # Directory for NanoPlot statistics
export FASTQC_FOLDER=/path/to/fastqc    # Directory for FastQC results
export MULTIQC_FOLDER=/path/to/multiqc  # Directory for aggregated MultiQC report

# Define sequencing kit information (update accordingly)
export KIT_INFO="SQK-NBD-24"
export PLATE="Plate01"
```

Make sure all tools are installed and available in your system `PATH`, or update the variables to their absolute paths.


### Dependencies

Install or load the following tools before execution:

| Tool | Purpose | Reference |
|------|----------|------------|
| [Dorado](https://github.com/nanoporetech/dorado) | Basecalling, demultiplexing, and summarization | Oxford Nanopore Technologies |
| [SAMTools](http://www.htslib.org/) | BAM file manipulation and FASTQ conversion | Heng Li et al. |
| [NanoPlot](https://github.com/wdecoster/NanoPlot) | Read length and quality visualization | Wouter De Coster |
| [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | Per-base quality assessment | Babraham Bioinformatics |
| [MultiQC](https://multiqc.info/) | Aggregation of QC reports | Phil Ewels et al. |

### Example directory structure

```bash
├── pod5/             # Raw input POD5 files
├── bam/              # Basecalled BAM outputs
├── demux/            # Demultiplexed BAM files
├── fastq/            # FASTQ outputs
├── summary/          # Dorado summary files
├── stats/            # NanoPlot results
├── fastqc/           # FastQC reports
└── multiqc/          # MultiQC aggregated reports
```

### Workflow
![Preprocessing](/.github/preprocessing.png)  

### Pipeline steps
#### Step 01. Basecalling with Dorado

```bash
mkdir -p $BAM_FOLDER
$DORADO/dorado basecaller sup $POD5_FOLDER/ONT_5e22_a7b0_1.pod5 --kit-name $KIT_INFO --barcode-both-ends > $BAM_FOLDER/Unaligned_1.bam
$DORADO/dorado basecaller sup $POD5_FOLDER/ONT_5e22_a7b0_2.pod5 --kit-name $KIT_INFO --barcode-both-ends > $BAM_FOLDER/Unaligned_2.bam
```

#### Step 02. Merge BAM files

```bash
samtools merge $BAM_FOLDER/Merged.bam $BAM_FOLDER/*.bam
```

#### Step 03. Demultiplex reads by barcode

```bash
mkdir -p $DEMUX_FOLDER
$DORADO/dorado demux --output-dir $DEMUX_FOLDER --no-classify $BAM_FOLDER/Merged.bam
```

#### Step 04. Generate sequencing summary

```bash
mkdir -p $SUMMARY_FOLDER
$DORADO/dorado summary $BAM_FOLDER/Merged.bam > $SUMMARY_FOLDER/$PLATE.tsv
```

#### Step 05. Convert BAM to FASTQ

```bash
mkdir -p $FASTQ_FOLDER
samtools fastq -T"*" $DEMUX_FOLDER/ONT_5e22_a7b0_1_barcode96.bam | gzip > $FASTQ_FOLDER/barcode96.fastq.gz
```

#### Step 06. Summarize read quality using NanoPlot

```bash
mkdir -p $STAT_FOLDER
NanoPlot --fastq $FASTQ_FOLDER/barcode96.fastq.gz --no_static --tsv_stats --raw --outdir $STAT_FOLDER/barcode96
```

#### Step 07. Evaluate per-base quality using FastQC

```bash
mkdir -p $FASTQC_FOLDER/barcode96
fastqc -o $FASTQC_FOLDER/barcode96 --extract -f fastq $FASTQ_FOLDER/barcode96.fastq.gz
```

#### Step 08. Aggregate reports with MultiQC

```bash
mkdir -p $MULTIQC_FOLDER
multiqc -z -o $MULTIQC_FOLDER $FASTQC_FOLDER
```

### Notes

- Each step can be **parallelized** for batch processing across multiple samples or barcodes.  
- Environment variables such as `$BAM_FOLDER`, `$DEMUX_FOLDER`, `$FASTQ_FOLDER`, etc., should be defined before execution.  
- The pipeline assumes **Dorado output in BAM format** and subsequent **conversion to FASTQ** for compatibility with standard QC tools.  

---
