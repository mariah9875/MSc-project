""" This readme file contains step-by-step procedures that were done in this analysis """

BEDTools/2.26.0

""" Raw fastq files were run with FASTQC program to check if adapters needed to be removed along with the quality of the reads. """


""" Raw fastq files were trimmed with cutadapt program
4 bases of each end:"""
 ls *.fastq | sed s/.fastq// | while read file; do cutadapt -u 4 -u -4 -o $file.trimmed.fastq $file.fastq; echo $file processed; done
 """trim adapter"""
 ls *.trimmed.fastq | sed s/.trimmed.fastq// | while read file; do cutadapt -a what to trim -o $file.trimmed_1.fastq --minimum-length length $file.trimmed.fastq; echo $file processed; done

"""
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


A file named config.yaml is created by the user, and the file has to include a list called
 "samples" with the names of the fastq sample files used, as shown below.
#config.yaml
samples: ["sample_1.fastq", "sample_2.fastq", "sample_3.fastq"]

The output files from the Snakemake pipeline are countmatrixes of .csv files, .bw and .bam files.
The .csv file can be used for DEA, the .bw files are uploaded to the UCSC genome browser and the bam files are used for the COMPSRA program to get countmatrix.

"""
""" COMPSRA program run with bam files, for each sample file 6 .txt annotation files will be produced for each small RNA type counted for.
repeat for each sample being counted """
java -jar /sw/pkg/compsra/1.0.2/COMPSRA.jar -ref hg38 -ann -ac 1,2,3,4,5,6 -in results/sorted_bam/version_2/*files.sorted.bam -out results/COMPSRA/version_2/*files_sorted

""" Merge COMPSRA annotation files from each sample into 1 .txt file for that small RNA type. A sample.list file is needed, which contains the pwd of the files being merged
together. Repeat 5 times for each type """
java -jar /sw/pkg/compsra/1.0.2/COMPSRA.jar -ref hg38 -fun -fm -fms 1-11 -fdclass 1 -pro COMPSRA_MERGE -inf results/COMPSRA/version_2/sample.list -out results/COMPSRA/version_2

""" The resulting 6 .txt countmatrix files are used as an input for the R script Data.prep.R. The resulting countmatrix generated from Data.prep.R is
used as an input for DEA.R.

The resulting up_regulated.txt and all_regulated_location.txt files generated from DEA.R were used for bedTools intersect, to overlap with TEs and SNPs locations """
bedtools intersect -a TEs.bed -b up_regulated_location.txt -wa -wb > TEs_overlapp.txt
bedtools intersect -a gwas.bed -b all_regulated_location.txt -wa -wb

""" The resulting .txt file from bedTools intersect was used as an input for R script TEintersect.R"""

SEtja inn líka fyrir multiqc
