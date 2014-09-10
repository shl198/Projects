#!/bin/bash -l
ref_path=/home/shangzhong/Database
chok1_fa=/opt/genome/cho/cgr_ref_CriGri_1.0_chrUn.fa
chok1_annotation=/opt/genome/cho/chok1.gff3
run_path=/data/shangzhong/VariantCall
fastq1=/data/shangzhong/VariantCall/L1.fq.gz
fastq2=/data/shangzhong/VariantCall/L2.fq.gz
sjdb=$run_path/CHOS_A4SJ.out.tab
picard=/home/shangzhong/Installation/picard-tools-1.113/picard-tools-1.113
gatk=/home/shangzhong/Installation/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar
declare t=24
cd $run_path
#============  0.generate reference ================================
#java -jar $picard/CreateSequenceDictionary.jar R=$chok1_fa O=chok1.dict
#samtools faidx $chok1_fa 
#============  1.mapping to reference   ============================
#----- (1) index genome reference -----------
#STAR --runMode genomeGenerate --genomeDir $ref_path/STAR_chok1Db --genomeFastaFiles $chok1_fa --runThreadN $t --limitGenomeGenerateRAM 310000000000 --sjdbGTFfile $chok1_annotation --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 100
#----- (2) align reads ----------------------
#STAR --genomeDir $ref_path/STAR_chok1Db --readFilesCommand zcat --readFilesIn $fastq1 $fastq2 --runThreadN $t --outFileNamePrefix CHOS_A4 
#----- (3) generate genome using splice junction -------
#STAR --runMode genomeGenerate --genomeDir $ref_path/STAR_chok1Db --genomeFastaFiles $chok1_fa --runThreadN $t --limitGenomeGenerateRAM 310000000000 --sjdbFileChrStartEnd $sjdb --sjdbOverhang 100
#----- (4) 2nd-pass mapping -------------------
#STAR --genomeDir $ref_path/STAR_chok1Db --readFilesCommand zcat --readFilesIn $fastq1 $fastq2 --runThreadN $t --outFileNamePrefix CHOS_A4
#===== 2. Add read groups,sort,mark duplicates,and create index ====
#----- (1) group reads ------------------------
#java -jar $picard/AddOrReplaceReadGroups.jar I=CHOS_A4Aligned.out.sam O=rg_added_sorted.bam SO=coordinate RGID=L2CHOSA4 RGLB=CHOSA4 RGPL=illumina RGPU=machine RGSM=CHOS
#----- (2) mark duplicates --------------------
#java -jar $picard/MarkDuplicates.jar I=rg_added_sorted.bam O=dedupped.bam M=output.metrics CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT
#===== 3. Split 'N' Trim and reassign mapping qualiteies ===========
#java -jar $gatk -T SplitNCigarReads -R $chok1_fa -I dedupped.bam -o split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS -fixMisencodedQuals
#===== 4. Indel Realignment ========================================
#--------- (1) Realigner Target Creator -------
#java -Xmx2g -jar $gatk -T RealignerTargetCreator -R $chok1_fa -I split.bam -o forIndelRealigner.intervals
#--------- (2) Indel Realigner ----------------
#java -Xmx4g -jar $gatk -T IndelRealigner -R $chok1_fa -I split.bam -targetIntervals forIndelRealigner.intervals -o realignedBam.bam 
#===== 5. Base recalibration =======================================
#   I cannot do this step since I don't have dbSNP for CHO, this step
#   only has a marginal effect on the result
#--------- (1) Base recalibrator --------------
# cannot generate because we don't have the dbSNP for CHO.
#java -Xmx4g -jar $gatk -T BaseRecalibrator -I split.bam -R $chok1_fa -o recal_data.table
#java -jar $gatk -T PrintReads -R $chok1_fa -I realignedBam.bam -BQSR recalibration_report.grp -o recali.bam
#===== 6. Variant calling ===========================================
#java -jar $gatk -T HaplotypeCaller -R $chok1_fa -I realignedBam.bam -recoverDanglingHeads -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 -o variCall.vcf
#===== 7. Variant filtering =========================================
java -jar $gatk -T VariantFiltration -R $chok1_fa -V variCall.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o fianlVariCall.vcf
