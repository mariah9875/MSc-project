# This readme file contains step-by-step procedures that were done in this analysis

# Raw fastq files were run with FASTQC program to check if adapters needed to be removed along with the quality of the reads.


# Raw fastq files were trimmed with cutadapt program
# 4 bases of each end:
 ls *.fastq | sed s/.fastq// | while read file; do cutadapt -u 4 -u -4 -o $file.trimmed.fastq $file.fastq; echo $file processed; done
 #trim adapter
 ls *.trimmed.fastq | sed s/.trimmed.fastq// | while read file; do cutadapt -a what to trim -o $file.trimmed_1.fastq --minimum-length length $file.trimmed.fastq; echo $file processed; done

Trimmed reads ran with Snakefile:
The snakemake pipeline takes the trimmed files and performs the follwing steps:
  QC: FastQC and MultiQC
  Align reads to a reference genome (hg38) with STAR aligner Mapping 1 (M1) & Mapping 2 (M2)
  Sort and index bam files for both M1 & M2
  Quantification of reads with featureCounts for M1 & M2
  Transcription assembly and quantification for M1
  Convert bam files to BigWig files for M1 & M2

Directory structure to run snakemake pipeline; the Snakefile is suppose to be in the same directory as the data and bin directories.
Example:
  projectdir/

        Snakefile

        config.yaml

        GTF annotation file

        data/

        bin/


# A file named config.yaml is created by the user, and the file has to include a list called
# "samples" with the names of the fastq sample files used, as shown below.
#config.yaml
  samples: ["sample_1.fastq", "sample_2.fastq", sample_3.fastq"]
