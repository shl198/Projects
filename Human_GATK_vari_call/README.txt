GATK pipeline for Human.

1. Required software
====================

* fastqc
* python
* trimmomatic
* bwa
* samtools
* picard tools
* gatk
* RStudio and gplots, ggplot2, reshape, gsalib


2. Steps before pipeline
========================

* run fastqc for one fastq file, if the Encoding version is smaller than 1.8 then trimming is necessary, and phred should be set as 64.

* build index for aligner bwa or STAR

* download necessary data on website ftp://ftp.broadinstitute.org/bundle/2.8/b37/ and unzip them.
  Files include:
	* human_g1k_v37.fasta.gz and human_g1k_v37.fasta.fai.gz
	* human_g1k_v37.dict.gz
	* human_g1k_v37.dict.gz and human_g1K_v37.dict.idx.gz
	* 1000G_phase1.indels.b37.vcf.gz and 1000G_phase1.indels.b37.vcf.idx.gz
	* dbsnp_138.b37.vcf and dbsnp_138.b37.vcf.idx.gz
	* hapmap_3.3.b37.vcf and hapmap_3.3.b37.vcf.idx.gz
	* Mills_and_1000G_gold_standard.indels.b37.vcf.gz and Mills_and_1000G_gold_standard.indels.b37.vcf.idx.gz
		

3. Steps to run the pipeline
============================

1) put all fastq files into one folder. (for paired end fastq file, make sure file names end with 1.fq.gz, 2.fq.gz, or 1.fastq.gz, 2.fastq.gz)

2) Define all parameters in the GATK_parameters.txt file.

3) In terminal run the code:
	python /pathway/to/Human_GATK_DNA_vari_call.py /pathway/to/DNA_parameters.txt > log.txt
or	nohup python /pathway/to/Human_GATK_DNA_vari_call.py /pathway/to/DNA_parameters.txt > log.txt & 
	(nohup helps keep running process when connection to server is closed)

4) In the log.txt file, the bottom will show the filename that stores the final results.
 
