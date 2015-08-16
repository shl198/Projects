import pandas as pd
import subprocess
from Modules.f11_snpEff_provean import snpEff_annotateVCF

def snpSift_filterVCF(annotatedVCF,snpSift,filters):
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
           '\"ANN[*].FEATUREID\" \"ANN[*].HGVS_C\" \"ANN[*].HGVS_P\" > {output}').format(snpSift=snpSift,
            filterCmd=filterCmd,inputVcf=annotatedVCF,sep='\',\'',empty='\'.\'',output = outputFile)
    subprocess.call(cmd,shell=True)
    
    return outputFile

def IndelAnalysis(vcf,snpEff,snpSift,genome,geneFile):
    """
    This function annotates and extracts the target genes' variants
    
    * vcf: str. VCF file name
    * snpEff: str. Pathway to snpEff
    * snpSift: str. Pathway to snpSift
    * genome: str. Genome name for snpEff. eg: 'chok1201405'
    * geneFile: str. A file has one column, list all target genes
    """
    annotatedVCF = snpEff_annotateVCF(vcf,snpEff,genome)
    # extract gene indels 
    df = pd.read_csv(geneFile,header=0,sep='\t',names=['GeneName'])
    genes = df['GeneName'].tolist()
    res_df = pd.DataFrame()
    for gene in genes:
        gene_if = ('(ANN[*].GENE=\'{gene}\')').format(gene=gene)
        filteredVCF = snpSift_filterVCF(annotatedVCF,snpSift,[gene_if])
        df = pd.read_csv(filteredVCF,sep='\t',header=0)
        if df.empty:
            continue
        else:
            res_df = res_df.append(df,ignore_index=True)
    # separate the items, keep only one item in each field
    res = vcf[:-31] + '.final.csv'
    outHandle = open(res,'w')
    vcf_df = res_df
    for i in range(len(vcf_df['#CHROM'])):
        chrom = vcf_df['#CHROM'][i]
        effects = vcf_df['ANN[*].EFFECT'][i].split(',')
        genes = vcf_df['ANN[*].GENE'][i].split(',')
        impacts = vcf_df['ANN[*].IMPACT'][i].split(',')
        trids = vcf_df['ANN[*].FEATUREID'][i].split(',')
        HGVS_C = vcf_df['ANN[*].HGVS_C'][i].split(',')
        HGVS_P = vcf_df['ANN[*].HGVS_P'][i].split(',')
        for effect,gene,impact,trid,hgc,hgp in zip(effects,genes,impacts,trids,HGVS_C,HGVS_P):
            if (hgp == '.') or (gene not in genes):
                continue
            else:
                outHandle.write('\t'.join([chrom,effect,gene,impact,trid,hgc,hgp])+'\n')
    outHandle.close()
    df = pd.read_csv(res,sep='\t',header=0)
    df = df.drop_duplicates()
    df.to_csv(res,index=False,sep='\t')

# vcf = '/data/shangzhong/VariantCall/GlycoKO/127/filterResults/127.merged.filter.indel.recode.vcf'
# snpEff = '/data/shangzhong/VariantCall/snpEff/snpEff.jar'
# snpSift = '/data/shangzhong/VariantCall/snpEff/SnpSift.jar'
# genome = 'chok1201405'
# geneFile = '/data/shangzhong/VariantCall/GlycoKO/gs_genes.txt'
# IndelAnalysis(vcf,snpEff,snpSift,genome,geneFile)


