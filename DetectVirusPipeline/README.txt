Detect Virus pipeline

1. Required software
====================

* python (should also install pandas, numpy packages)
* fastqc
* trimmomatic
* gsnap or STAR
* samtools
* picard tools (version should be lower than 1.123)
* blast+
* trinity


2. Steps before pipeline
========================

* run fastqc for one fastq file, if the Encoding version is smaller than 1.8 then trimming is necessary, and phred should be set as 64.
* download virus genome sequence at ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/all.fna.tar.gz
* download virus genome annotation at ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/all.gff.tar.gz
* download host genome sequence either at ncbi or ensembl
* if use gsnap: build host genome index for gsnap by running:
	gmap_index -D /path/to/store/the/index -d organism_name fastaFile.fa
* build virus genome index using blast

3. Steps to run the pipeline
============================

1) put all fastq files into one folder. (for paired end fastq file, make sure file names end with _1.fq.gz, _2.fq.gz, or _1.fastq.gz, _2.fastq.gz)

2) Define all parameters in the 01_DetectVirus_Parameters.txt file.

3) In terminal run the code:
	python /pathway/to/01_DetectVirus.py /pathway/to/01_DetectVirus_Parameters.txt > log.txt
or	python /pathway/to/01_DetectVirus.py /pathway/to/01_DetectVirus_Parameters.txt > log.txt & 
	(nohup helps keep running process when connection to server is closed)
4) If choose not to run blast after assembly, the final result would be the file Trinity.fasta in folder trinity_out_dir. If chosse to run blast after assembly, teh final result would be the file Trinity.blast.txt in folder trinity_out_dir
