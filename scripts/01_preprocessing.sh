#### #### #### #### #### #### #### #### #### #### #### #### 
# ONT Data Preprocesing Pipeline
# by Heewon Seo <Heewon.Seo@ucalgary.ca>
# on October 28, 2025
#### #### #### #### #### #### #### #### #### #### #### #### 
## The following script is provided for demonstration purposes only. 
## Dorado is optimized for GPU usage, whereas SAMTools, FastQC, and MultiQC do not require a GPU. 
## Additionally, a few steps below require parallelization, which is not reflected in this example.
#### #### #### #### #### #### #### #### #### #### #### #### 
# Settings
PLATE="Plate10"
KIT_INFO="SQK-NBD-96"

DORADO="~/tools/dorado-1.1.1-linux-x64/bin/"
POD5_FOLDER="~/RAWDATA/Benchmarking/POD5/$PLATE"
BAM_FOLDER="~/projects/Benchmarking/BAM/$PLATE"
DEMUX_FOLDER="~/projects/Benchmarking/DEMUX/$PLATE"
FASTQ_FOLDER="~/projects/Benchmarking/FASTQ/$PLATE"
SUMMARY_FOLDER="~/projects/Benchmarking/Summary"
STAT_FOLDER="~/projects/Benchmarking/Stats/$PLATE"
FASTQC_FOLDER="~/projects/Benchmarking/FASTQC/$PLATE"
MULTIQC_FOLDER="~/projects/Benchmarking/MULTIQC/$PLATE"

#### #### #### #### #### #### #### #### #### #### #### #### 
# [Step 01] Basecalling with Dorado requires parallelization to process multiple POD5 files efficiently
## The code below demonstrates an example of basecalling for two POD5 files
mkdir -p $BAM_FOLDER
$DORADO/dorado basecaller sup $POD5_FOLDER/ONT_5e22_a7b0_1.pod5 --kit-name $KIT_INFO --barcode-both-ends > $BAM_FOLDER/Unaligned_1.bam
$DORADO/dorado basecaller sup $POD5_FOLDER/ONT_5e22_a7b0_2.pod5 --kit-name $KIT_INFO --barcode-both-ends > $BAM_FOLDER/Unaligned_2.bam

# [Step 02] Merge individual BAM files into a single consolidated BAM file
samtools merge $BAM_FOLDER/Merged.bam $BAM_FOLDER/*.bam

# [Step 03] Demultiplex the merged BAM file into separate barcode-specific BAM files
mkdir -p $DEMUX_FOLDER
$DORADO/dorado demux --output-dir $DEMUX_FOLDER --no-classify $BAM_FOLDER/Merged.bam

# [Step 04] Generate a summary of the sequencing output from the merged BAM file
mkdir -p $SUMMARY_FOLDER
$DORADO/dorado summary $BAM_FOLDER/Merged.bam > $SUMMARY_FOLDER/$PLATE.tsv

# [Step 05] Convert the BAM files to FASTQ format - parallelization
## The code below demonstrates an example of converting a single BAM file
mkdir -p $FASTQ_FOLDER
samtools fastq -T"*" $DEMUX_FOLDER/ONT_5e22_a7b0_1_barcode96.bam | gzip > $FASTQ_FOLDER/barcode96.fastq.gz

# [Step 06] Generate a summary of the sequencing output for each barcode-specific BAM file - parallelization
mkdir -p $STAT_FOLDER
NanoPlot --fastq $FASTQ_FOLDER/barcode96.fastq.gz --no_static --tsv_stats --raw --outdir $STAT_FOLDER/barcode96

# [Step 07] Assess per-base quality and other sequencing metrics using FastQC - parallelization
mkdir -p $FASTQC_FOLDER/barcode96
fastqc -o $FASTQC_FOLDER/barcode96 --extract -f fastq $FASTQ_FOLDER/barcode96.fastq.gz

# [Step 08] Aggregate all FastQC reports using MultiQC
mkdir -p $MULTIQC_FOLDER
multiqc -z -o $MULTIQC_FOLDER $FASTQC_FOLDER

#### #### #### #### #### #### #### #### #### #### #### #### 
# FIN
