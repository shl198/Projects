import sys,os
import subprocess
from Bio import SeqIO
import pandas as pd

def annotateVCF(vcf,snpEff,genome):
    """
    This function annotation vcf using snpEff,
    return annotated eff file
    
    * vcf: vcf file that need to be annotated
    * snpEff: snpEff software
    * genome: genome database built using snpEff 
    """
    configure = snpEff[:-3] + 'config'
    output = vcf[:-3] + 'eff.vcf'
    cmd = ('java -jar {snpEff} -c {configure} -v {genome} {vcf} > {output}').format(
                    snpEff=snpEff,configure=configure, genome=genome,vcf=vcf,output=output)
    subprocess.call(cmd,shell=True)
    return output

def filterVCF(annotatedVCF,snpSift,filters):
    """
    This function filter the vcf inputwith snpSift
    
    *vcf: annotated vcf file
    *snpSift: pathway to snpSift
    *filters': a list of arguments used to filter the file
    """
    outputFile = annotatedVCF[:-3] + 'target.vcf'    # variant only for interested genes
    filterCmd = ' '.join(filters)
    filterCmd = '\"' + filterCmd + '\"'
    # extract variants of target genes.
    cmd = ('java -jar {snpSift} filter {filterCmd} {inputVcf} | java -jar {snpSift} extractFields -s {sep} -e {empty} - '
           'CHROM \"ANN[*].EFFECT\" \"ANN[*].GENE\" \"ANN[*].IMPACT\" '
           '\"ANN[*].FEATUREID\" \"ANN[*].HGVS_P\" > {output}').format(snpSift=snpSift,
            filterCmd=filterCmd,inputVcf=annotatedVCF,sep='\',\'',empty='\'.\'',output = outputFile)
    subprocess.call(cmd,shell=True)
    return outputFile

def complement_reverse(sequence):
    """
    This fucntion gets the complement reverse sequence of a sequence
    
    * sequence: a string of sequence in ATCG format
    """
    dic = {'A':'T','T':'A','C':'G','G':'C'}
    result = ''
    for letter in sequence:
        result = dic[letter] + result
    return result

def RNA2Pr(RNA,CodonFile):
    """
    This function translate RNA to protein
    
    * RNA: rna sequence
    * CodonFile: filename. It has codon table
    """
    codon = {}
    res = open(CodonFile,'r')
    # build codon library
    for line in res:
        if line[-1] == '\n':
            item = line[:-1].split(' ')
        else:
            item = line.split(' ')
        try:
            codon[item[0]] = item[1]
        except:
            codon[item[0]] = ''
    # decode
    i = 0
    code = RNA[i:i+3]
    amino = codon[code]
    while codon[code] != '':
        i = i + 3
        try:
            code = RNA[i:i+3]
            amino = amino + codon[code]            
        except:
            break
    return amino


def get_AA_from_gff(refDNA_dic,gffFile,chrom,transcript,CodonFile):
    """
    This function extracts amino acid sequence from gff annotation file, given 
    reference sequence dictionary, chromosome name, gene name, transcript id
    and CodonFile.
    
    * refDNA_dic: dictionary get from command record_dict = SeqIO.index(fastaFile,'fasta')
    * gffFile: gff annotation file
    * chrom: chromosome name
    * gene: gene symbol
    * transcript: transcript id
    * CodonFile: file has two columns. 1st is codon, 2nd is amino acids.
    
    """
    # get chromosome sequence
    chr_seq = refDNA_dic[chrom].seq
    # get CDS, the coding sequence
    cmd = ('grep Parent={trans_id} {gff} | awk -F \'\\t\' \'{ifcmd} {printrow}\'').format(trans_id=transcript,gff=gffFile,
                                                                                          ifcmd='$3==\"CDS\"',printrow='{print $4FS$5}')
    exon_loci = subprocess.check_output(cmd,shell=True)
    exon_loci_list = exon_loci[:-1].split('\n')
    # combine CDS together
    start_loci = []; stop_loci = []
    for item in exon_loci_list:
        locis = item.split('\t')
        start_loci.append(locis[0])
        stop_loci.append(locis[1])
    transcript_seq = ''
    for i,j in zip(start_loci,stop_loci):
        i = int(i); j = int(j)
        if int(start_loci[0]) < int(start_loci[1]): # this means the transcript is in positive strand
            transcript_seq = transcript_seq + chr_seq[i-1:j]
        else:  # this transcript is in negtive strand
            transcript_seq = transcript_seq + complement_reverse(chr_seq[i-1:j])
            
    # get the amino acid sequence
    AA = RNA2Pr(transcript_seq,CodonFile)
    return AA

def vcf2input4provean(record_dict,gff3File,vcfFile,gene):
    """
    This function will get two input files for provean: [protein sequence fasta file, protein variant txt file]
    
    * record_dict: dictionary for reference DNAseq
    * gff3File: annotation file
    * vcfFile: variant result file
    * gene: gene name
    """
    outputFile = vcfFile[:-3] + 'target.vcf'    # variant only for interested genes
    # extract variants of target genes.
    cmd = ('java -jar {snpSift} filter \"ANN[*].GENE=\'{gene}\'\" {inputVcf} | java -jar {snpSift} extractFields -s {sep} -e {empty} - '
           'CHROM \"ANN[*].EFFECT\" \"ANN[*].IMPACT\" \"ANN[*].FEATUREID\" \"ANN[*].HGVS_P\" > {output}').format(snpSift=snpSift,
                gene=gene,inputVcf=vcfFile,sep=',',empty='.',output = outputFile)
    subprocess.call(cmd,shell=True)
    
    p_vari_dic = {}   # a dictionary, each item is   transcript id:[protein variant]
    result = open(outputFile,'r')      # filename.eff.target.vcf
    protein_fasta = vcfFile[:-3] + 'protein.fa'
    protein_variant = vcfFile[:-3] + 'variant.txt'
    output_protein = open(protein_fasta,'w')
    output_variant = open(protein_variant,'w')
    next(result)
    for line in result:
        item = line[:-1].split('\t')
        chrom = item[0]
        trids = item[2].split(',')
        HGVS_P = item[3].split(',')
        for trid,pr_vari in zip(trids,HGVS_P): # iterate protein variant
            if pr_vari == '.':
                continue
            else:
                try:
                    p_vari_dic[trid].append(pr_vari)
                except:
                    p_vari_dic[trid] = [pr_vari]
            # get AA sequence
            AA = get_AA_from_gff(record_dict,gffFile,chrom,gene,trid,CodonFile) 
            output_variant.write(gene + '_' + trid + '\t' + ','.join(p_vari_dic[trid]) + '\n')
            output_protein.write('>' + gene + '_' + trid + '\n' + AA + '\n')
    output_variant.close()
    output_protein.close()

    return [protein_fasta,protein_variant]


fastaFile = '/opt/genome/hamster/hamster.fa'
record_dict = SeqIO.index(fastaFile,'fasta')

gffFile = '/opt/genome/hamster/hamster.gff3'
gene = 'Rad51b'
CodonFile = '/data/shangzhong/VariantCall/RNA_codon_table.txt'

vcfFile = '/data/shangzhong/VariantCall/CHOS_ToT_RNA_seq/vcfResult/CHOS_Tot_4A/A4.merged.filter.vcf'
snpSift = '/data/shangzhong/snpEff/SnpSift.jar'
snpEff = '/data/shangzhong/snpEff/snpEff.jar'
genome = 'hamster201405'
#===============================================================================
#        Variant analysis pipeline
#===============================================================================
#============= 1. Annotate vcf results using snpEff ================
#annotated = annotateVCF(vcfFile,snpEff,genome)  # annotated: filename.eff.vcf
#============= 2. Filter the annotated file ========================
annotated = '/data/shangzhong/VariantCall/CHOS_ToT_RNA_seq/vcfResult/CHOS_Tot_4A/A4.merged.filter.eff.vcf'
filtered = filterVCF(annotated,snpSift,['(ANN[*].GENE=\'Rad51b\')','&','((ANN[*].IMPACT=\'HIGH\') | '
                     '(ANN[*].IMPACT=\'MODERATE\'))'])
#============= 3. Get input files for provean ======================
# [proteinFile,variantFile] = vcf2input4provean(record_dict,gffFile,annotated,gene)
# 3. Run provean

os.chdir('/data/shangzhong/VariantCall/CHOS_ToT_RNA_seq/vcfResult/CHOS_Tot_4A')
vcf_df = pd.read_csv(filtered,sep='\t')
# iterate lines of filtered vcf results
pro_vari_dic = {}   # dictionary. Format: transcriptID:[protein variants]
trid_gene_dic = {}  # dictionary. Format: transcriptID:[gene symbol]
for i in range(len(vcf_df['#CHROM'])):
    chrom = vcf_df['#CHROM'][i]
    impacts = vcf_df['ANN[*].IMPACT'][i].split(',')
    trids = vcf_df['ANN[*].FEATUREID'][i].split(',')
    HGVS_P = vcf_df['ANN[*].HGVS_P'][i].split(',')
    genes = vcf_df['ANN[*].GENE'][i].split(',')
    for impact,trid,vari,gene in zip(impacts,trids,HGVS_P,genes):
        if vari == '.':
            continue
        else:
            if impact == 'HIGH' or impact == 'MODERATE':
                try:
                    pro_vari_dic[trid].append(vari)
                except:
                    pro_vari_dic[trid] = [vari]
                # add trid:gene to the dictionary
                if trid not in trid_gene_dic:
                    trid_gene_dic[trid] = gene
            else:
                continue
for trid in pro_vari_dic:
    gene = trid_gene_dic[trid]
    protein_file = gene + '_' + trid + '.protein.fa'
    vari_file = gene + '_' + trid + 'variant.txt'
    protein_output = open(protein_file,'w')
    AA = get_AA_from_gff(record_dict,gffFile,chrom,trid,CodonFile)
    protein_output.write(('>{gene_trid} \\n {AA}').format(gene_trid=gene+'_'+trid,AA=AA))
    protein_output.close()
    
    variant_output = open(vari_file,'w')
    for vari in pro_vari_dic[trid]:
        variant_output.write(vari + '\n')
    variant_output.close()    
            
print 'done'