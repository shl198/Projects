DE pipeline.

1. Required software
====================

* fastqc 
* python 
* trimmomatic
* STAR
* samtools
* picard tools (version should be lower than 1.121)


2. Steps before pipeline
========================

* run fastqc for one fastq file, if the Encoding version is smaller than 1.8 then trimming is necessary, and phred should be set as 64.

* If it detects any adapter, put the adapter into a fasta file format and provide to trimmomatic adapter file in the parameter file.


3. steps to run the quantification pipeline
===========================================

1) put all fastq files into one folder. (for paired end fastq file, make sure file names end with _1.fq.gz, _2.fq.gz, or _1.fastq.gz, _2.fastq.gz).

2) Define all parameters in the 01_DE_Parameters.txt.


4. steps to run the DE pipeline
===============================

1) In the 04_Diff_pairs.csv, define which samples you want to compare. The column is flexible, you can add/delete control and test columns. The input should be the index of the filename sorted by nature.

2) in terminal run:   Rscript 04_Diff_DESeq2.R /path/to/your/htseqcount/files

