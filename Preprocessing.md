# Preprocessing
> Advancing the preprocessing of Nanopore sequencing data

## Rationale

_MinKNOW_ from Oxford Nanopore frequently crashes and does not provide the most up-to-date _Dorado_ basecaller. Therefore, establishing an external computational environment is recommended for scalability and reproducibility. The _Dorado_ version reported in the Web summary from MinKNOW is the _Dorado Basecall Server_ version, which is misleading. In addition, the latest _Dorado_ release (v1.X.X), following the major update, provides an expanded quality score range based on improved signal-to-sequence predictors.

## Workflows

The pipeline performs the following tasks:

1. **Basecalling** of raw POD5 files using Dorado.  
2. **Merging** individual BAM files into a single dataset.  
3. **Demultiplexing** by barcode.  
4. **Generating sequencing summaries** for merged and demultiplexed BAM files.  
5. **Converting BAM to FASTQ** format.  
6. **Computing quality statistics** using NanoPlot and FastQC.  
7. **Aggregating reports** with MultiQC.
8. **Filtering and trimming reads** by length.
9. **Removing host reads** (mouse and human) with minimap2.

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
export FILTERED_FOLDER=/path/to/filtered # Directory for length-filtered FASTQ files
export HOSTFREE_FOLDER=/path/to/hostfree # Directory for host-free FASTQ files
export HOST_REF=/path/to/host_reference  # Directory for the combined host reference and minimap2 index

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
| [cutadapt](https://cutadapt.readthedocs.io/) | Trimming reads to a specified length | Martin |
| [minimap2](https://github.com/lh3/minimap2) (v2.30) | Alignment of reads to the host genomes | Heng Li |

### Example directory structure

```bash
├── pod5/             # Raw input POD5 files
├── bam/              # Basecalled BAM outputs
├── demux/            # Demultiplexed BAM files
├── fastq/            # FASTQ outputs
├── summary/          # Dorado summary files
├── stats/            # NanoPlot results
├── fastqc/           # FastQC reports
├── multiqc/          # MultiQC aggregated reports
├── filtered/         # Length-filtered FASTQ files
├── host_reference/   # Combined host genome and minimap2 index
└── hostfree/         # Filtered, host-free FASTQ files (input to the MAG and classification modules)
```

### Schematic overview
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

#### Step 09. Filter and trim reads by length

The selection of filtering thresholds is inherently arbitrary and can substantially influence downstream analyses. Therefore, we recommend that users justify cutoff values based on empirical evidence derived from their own data.

To inform our filtering strategy, we evaluated per-base read quality across read-length bins defined as follows: 1–999 bp, 1,000–9,999 bp, 10,000–99,999 bp (long reads), and ≥100,000 bp (ultra-long reads). This analysis revealed two notable patterns: (1) ultra-long reads exhibited a pronounced decline in base quality after the first 5–10 bases, and (2) long reads showed a marked decrease in base quality beyond approximately 15,000 bp. Based on these observations, we chose to discard ultra-long reads entirely and to trim long reads to a maximum length of 15,000 bp.

With respect to read quality filtering, Dorado discards reads with a mean Phred quality score below 20. Given that our length-based filtering and trimming strategy substantially improved overall read quality, we did not apply additional quality-based filtering steps.

The code below removes ultra-long reads and trims long reads according to the criteria described above, starting from FASTQ files split into the length bins above.

```bash
mkdir -p $FILTERED_FOLDER
cutadapt -l 14999 -o $FILTERED_FOLDER/barcode96_10000-15000.fastq.gz $FASTQ_FOLDER/barcode96_10000-99999.fastq.gz
# barcode96_gtet100000.fastq.gz (gtet, greater than or equal to) is not carried forward
cat $FASTQ_FOLDER/barcode96_1-999.fastq.gz $FASTQ_FOLDER/barcode96_1000-9999.fastq.gz $FILTERED_FOLDER/barcode96_10000-15000.fastq.gz > $FILTERED_FOLDER/barcode96.fastq.gz
```

#### Step 10. Remove host reads

In human microbiota-associated mouse models, reads from the mouse host (and from human contamination) should be removed before assembly and taxonomic classification. Otherwise, host reads enter the assemblies and compete for _k_-mer assignments during classification. Here, the GENCODE primary assemblies of mouse (GRCm39) and human (GRCh38) are combined, with sequence names prefixed by species so that host reads can be counted per species, and any read with a primary alignment to either genome is removed.

**Keep the read headers.** Dorado writes basecalling information as tags in each FASTQ header, and Medaka 2.1.0 reads these tags to check that its model matches the basecaller. Use `minimap2 -y` to copy the tags into the alignment and `samtools fastq -T '*'` to write them back. Without both options the host-free reads have no tags, and Medaka stops with `TypeError: argument of type 'NoneType' is not iterable`.

```bash
# Build the combined host reference once
mkdir -p $HOST_REF
sed 's/^>/>mm_/' GRCm39.primary_assembly.genome.fa >  $HOST_REF/GRCm39_GRCh38.fa
sed 's/^>/>hs_/' GRCh38.primary_assembly.genome.fa >> $HOST_REF/GRCm39_GRCh38.fa
minimap2 -x map-ont -d $HOST_REF/GRCm39_GRCh38.mmi $HOST_REF/GRCm39_GRCh38.fa

# Align each sample and keep the unmapped (host-free) reads, with header tags
mkdir -p $HOSTFREE_FOLDER
minimap2 -y -ax map-ont --secondary=no $HOST_REF/GRCm39_GRCh38.mmi $FILTERED_FOLDER/barcode96.fastq.gz | samtools view -b -o barcode96.host.bam -
samtools fastq -f 4 -T '*' barcode96.host.bam | gzip > $HOSTFREE_FOLDER/barcode96.fastq.gz

# Count host reads per species (primary alignments only)
samtools view -F 0x904 barcode96.host.bam | cut -f3 | grep -c '^mm_'   # mouse
samtools view -F 0x904 barcode96.host.bam | cut -f3 | grep -c '^hs_'   # human
```

### Notes

- Each step can be **parallelized** for batch processing across multiple samples or barcodes.  
- Environment variables such as `$BAM_FOLDER`, `$DEMUX_FOLDER`, `$FASTQ_FOLDER`, etc., should be defined before execution.  
- The pipeline assumes **Dorado output in BAM format** and subsequent **conversion to FASTQ** for compatibility with standard QC tools.  

---
