Differential Expression Analysis Pipeline
=========================================

Pipeline Description
--------------------

h3. Diff_Expression.py

Have the pipeline: Trim reads --> Alignment --> sort by samtools --> htseq count --> DESeq

Python file Description
-----------------------

* 01_DE_pipeline.py: the pipeline for getting htseq count
* 02_txedo_pipeline.py: txedo pipeline for only getting the fpkm step
* 03_sailfish_pipeline.py: saifish pipeline for geetigng the expression data

R file Description
------------------

* 04_Diff_DESeq2.R: main file to do DE analysis
* 04_Diff_DESeq2_QC.R: quality control before doing the analysis.
* 04_Diff_pairs.csv: files define control samples and test samples. Each line is one comparision DESeq2 will do.
