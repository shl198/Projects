Ribosome Profiling Project
==========================

01_RibosomePipeline.py
----------------------
Pipeline to preprocess fastq files to get the final bam file.

Method to run the command:
 * define all parameters in 01_RiboPro_Parameters.txt
 * In the python file, change the folder in sys.path.append('/home/shangzhong/Codes/Pipeline') to your local Pipeline
 * In terminal, run the following code:
 		nohup python 01_RibosomePipeline.py 01_RiboPro_Parameters.txt > log.txt &

 02_RibosomeRecalibration.py
 ---------------------------
 Calculate distance between 5' end and A site in ribosome for each alignment length.

03_GeneralAnalysis.py
---------------------
Some general analysis: coverage at each position of protein; total count for each gene.

04_RNA_process.py
-----------------
Process RNAseq reads.

05_RiboSignalpCoverPipeline.py
------------------------------
Analysis regarding the proteins with signal peptide.

06_Proteomaps.py
----------------
Analysis for preparing data for proteomaps.

f02_RiboDataModule.py
---------------------
Basic modules and functions the other files would use.
