Detect Virus pipeline

1. Required software
====================

* python
* fastqc
* trimmomatic
* gsnap
* samtools
* picard tools (version should be lower than 1.123)
* blast


2. Steps before pipeline
========================

* run fastqc for one fastq file, if the Encoding version is smaller than 1.8 then trimming is necessary, and phred should be set as 64.
* download virus genome seuqnce at ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/all.fna.tar.gz
* download virus genome annotation at ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/all.gff.tar.gz

* build host genome index for gsnap by running:
	gmap_index -D /path/to/store/the/index -d organism_name fastaFile.fa


3. Steps to run the pipeline
============================

1) do
