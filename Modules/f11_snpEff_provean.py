import sys
from Bio import SeqIO
from Modules.p05_ParseGff import *

#===============================================================================
#                 snpEff and snpSift
#===============================================================================
def snpEff_annotateVCF(vcf,snpEff,genome):
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
           '\"ANN[*].FEATUREID\" \"ANN[*].HGVS_P\" > {output}').format(snpSift=snpSift,
            filterCmd=filterCmd,inputVcf=annotatedVCF,sep='\',\'',empty='\'.\'',output = outputFile)
    subprocess.call(cmd,shell=True)
    
    return outputFile

#===============================================================================
#                 Provean analysis related functions
#===============================================================================
def vcf2input4provean(filteredVCF,record_dict,gffFile,CodonFile):
    """
    This function prepares input files for provean. protein.fa and variant.txt
    
    * filteredVCF: vcf file filtered using snpSift
    * record_dict: 
    """
    vcf_df = pd.read_csv(filteredVCF,sep='\t')
    if vcf_df.empty:
        return [[],[]]
    ##------- 1. define the dictionaries and file lists ------------
#     pro_vari_dic = {}   # dictionary. Format: transcriptID:[protein variants]
#     trid_gene_dic = {}  # dictionary. Format: transcriptID:[gene symbol]
#     trid_chrom = {}     # dictionary. Format: transcriptID:[chromosome]
    gene_trid2chrom = {}  # dictionary. Format: gene_transcriptID:chrome
    gene_trid2vari = {}   # dictionary. Format: gene_transcriptID:[protein variants]
    protein_files = [] # stores input files of protein sequences
    variant_files = [] # stores input files of variants
    ##------ 2. build the dictionaries. The reason for building the libraryies is that one ----------
    #    transcript may have many variants in different lines.
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
                # if trid is in the format of Transcript_genenumber, remove the transcript_
                if 'Transcript_' in trid:
                    trid = trid[11:] 
                
                gene_trid = gene + '_' + trid
                one_letter_vari = vari3letter2vari1letter(vari[2:])
                if impact == 'HIGH' or impact == 'MODERATE':
                    try:
                        gene_trid2vari[gene_trid].append(one_letter_vari)
                    except:
                        gene_trid2vari[gene_trid] = [one_letter_vari]
                    # add gene_trid:chrom to the dictionary
                    if gene_trid not in gene_trid2chrom:
                        gene_trid2chrom[gene_trid] = chrom
                else:
                    continue
    
    gene_tridList = sorted(gene_trid2vari.keys())
    #------- 3. generate a list of protein fa files and variant txt files ----
    for gene_trid in gene_tridList:
        item = gene_trid.split('_')
        gene = item[0]; trid = item[1]
        chrom = gene_trid2chrom[gene_trid]
        
        protein_file = gene_trid + '.protein.fa'
        vari_file = gene_trid + '.variant.txt'
        
        protein_files.append(protein_file)
        variant_files.append(vari_file)
        ### generate protein fa file
        protein_output = open(protein_file,'w')
        try:
            AA = get_AA_from_gff(record_dict,gffFile,chrom,gene,trid,CodonFile)
            protein_output.write(('>{gene_trid} \n{AA}').format(gene_trid=gene+'_'+trid,AA=AA))
        except:
            print 'fail to get',gene,'amino acid sequence.'
            raise
        protein_output.close()
        ### generate variant txt file
        variant_output = open(vari_file,'w')
        for vari in gene_trid2vari[gene_trid]:
            variant_output.write(vari + '\n')
        variant_output.close()    
    
    return [protein_files,variant_files]

def run_provean(provean_soft,protein,variant,support_set_path,support_set,thread = 1):
    """
    This function run provean program
    
    * provean_soft: pathway to provean.sh
    * protein: input fasta file which has query protein sequence
    * variant: input txt file which has variant represented by one letter amino acids
    * support_set_path: folder that has all the file.sss
    * support_set: a list of files in format of 'file.sss', if provided, provean would skip the blast step to save time.
    * thread: number of thread to run blast
    """
    output = protein[:-2] + 'sss'
    if output in support_set:
        cmd = ('{provean} -q {fasta} -v {vari} --num_threads {thread} --supporting_set {set}').format(
                    provean=provean_soft,fasta=protein,vari=variant,thread=str(thread),set=support_set_path+'/'+output)
    else:
        cmd = ('{provean} -q {fasta} -v {vari} --num_threads {thread} --save_supporting_set {output}').format(
                    provean=provean_soft,fasta=protein,vari=variant,thread=str(thread),output=output)
    result = subprocess.check_output(cmd.split(' '))
    return result

def capture_provean_scores(outputFile,provean,proteinFiles,variantFiles,support_set_path,support_set,thread = 1):
    """
    This function captures the provean scores originally output to terminal. Stores all the socres
    into a file.
    
    * outputFile: the file that stores all scores
    * provean: pathway to provean
    * protein_files: a list of protein fasta files
    * variant_files: a list of variant txt files
    * support_set_path: a folder that has all provided .sss file for provean
    * support_set: the .sss file name
    * thread: number of thread to run blast
    """
    output = open(outputFile,'w')
    output.write('# Gene_trID\tVARIATION\tSCORE\n')
    for proteinFile,variFile in zip(proteinFiles,variantFiles):
        result = run_provean(provean,proteinFile,variFile,support_set_path,support_set,thread)
        #==== 1. put all variants in vari into a list
        allVari = pd.read_csv(variFile,header=None,names=['Variant'])
        proveanVari = [] # proveanVari cannot process frameshift in allVari
        lines = result[:-1].split('\n')
        #==== 2. capture the output of provean
        for line in lines:
            if line.startswith('#') or line.startswith('[') or line.startswith('No'):
                continue
            else:
                item = line.split('\t')
                output.write(variFile[:-12] + '\t' + '\t'.join(item) + '\n')
                proveanVari.append(item[0])
        #==== 3. output the varis that provean cannot process
        for variant in allVari['Variant']:
            if variant not in proveanVari:
                output.write(variFile[:-12] + '\t' + variant + '\tNA' + '\n')
        print proteinFile,'provean analysis finishes\n'
    output.close()

# os.chdir('/data/shangzhong/VariantCall/CHOS_ToT_VS_chok1/filteredVcfResults/CHOS_TOT_A4')
# provean = '/home/shangzhong/Installation/provean-1.1.5/bin/provean.sh'
# proteinFiles = ['Nbn_rna17749.protein.fa']
# variantFiles = ['Nbn_rna17749.variant.txt']
# support_set_path = '/data/shangzhong/VariantCall/chok1_Provean_Supporting_Set'
# support_set = ['Nbn_rna17749.sss']
# outputFile = 'pro.txt'
# capture_provean_scores(outputFile,provean,proteinFiles,variantFiles,support_set_path,support_set,thread = 19)

def merge_all_provean_scores(previousFile,ref_name,*pathways):
    """
    This function combine all the provean results together in one table.
    if previousFile is provided, then the new result is appended to it.
    
    * previousFile: file that has provean results.
    * pathways: pathways that has folders for samples, with file named proveanScore.txt
    """
    columns = ['Gene_TRID','Variant','Score','Sample','Reference']
    if previousFile=='':
        df = pd.DataFrame(columns=columns)  # df stores all 5 columns information
    else:
        df = pd.read_csv(previousFile,sep='\t',header=None,
                         names=columns)
    for pathway in pathways:
        os.chdir(pathway)
        folders = [f for f in os.listdir(pathway) if os.path.isdir(os.path.join(pathway,f))]
        for folder in folders:
            f = folder + '/' + 'proveanScore.txt'
            data = pd.read_csv(f,sep='\t',header=0,names=['Gene_TRID','Variant','Score'])
            sample = [folder]*len(data['Score'])
            ref = [ref_name]*len(data['Score'])
            data['Sample'] = pd.Series(sample,index=data.index)
            data['Reference'] = pd.Series(ref,index=data.index)
            data = data.fillna(0)
            df = df.append(data,ignore_index=True)
    # get pivot results of column: gene_trID_vari, rows: samples, values: score.
    pivot = df[['Sample','Score']]
    pivot.insert(0,'Gene_TRID_Vari',df['Gene_TRID'].map(str) + '_' + df['Variant'])
    df_pivot = pd.pivot_table(pivot,values='Score',index='Gene_TRID_Vari',cols='Sample')
    # output whole results to file
    df.to_csv(previousFile+'merge.csv',sep='\t')
    df_pivot.to_csv(previousFile+'pivot.csv',sep='\t')
    

# merge_all_provean_scores('','hamster','/data/shangzhong/VariantCall/CHOSDNA_VS_hamster/finalVcfResult',
#                          '/data/shangzhong/VariantCall/CHOS_ToT_VS_hamster/filteredVcfResult',
#                          '/data/shangzhong/VariantCall/GlycoDeletion_VS_hamster/filteredVcfResults',
#                          '/data/shangzhong/VariantCall/GlycoMutant_VS_hamster/filteredVcfResults',
#                          '/data/shangzhong/VariantCall/hamster_VS_hamster/final_result',
#                          '/data/shangzhong/VariantCall/pgsa_VS_hamster',
#                          '/data/shangzhong/VariantCall/pgsaDNA_VS_hamster/filteredVcfResults'
#                          )
# print 'done'

def merge_all_genes_with_diff_pr_from_refseq(previousFile,outputFile,commonFileName,*pathways):
    """
    This function merge all transcripts in gff annotation file that has different protein sequence
    from refseq protein sequence
    
    * previousFile: fasta file that will be added sequence
    * ref_name: the organism name
    * *pathways: pathways that have folders of samples, which have name LocalOnlineProteinDiff.fa
    """
    columns = ['Gene','TrID']
    if previousFile=='':
        df = pd.DataFrame(columns=columns)  # df stores 2 columns
    else:
        df = pd.read_csv(previousFile,sep='\t',header=None,
                         names=columns)
    # 1. store all the description into one list
    for pathway in pathways:
        os.chdir(pathway)
        folders = [f for f in os.listdir(pathway) if os.path.isdir(os.path.join(pathway,f))]
        description = []
        for folder in folders:
            f = folder + '/' + commonFileName
            sequence = SeqIO.parse(open(f,'rU'),'fasta')
            for seq in sequence:
                des = seq.description
                if des not in description:
                    description.append(seq.description)
    # 2. extract gene name and transcript id
    for des in description:
        # get transcript id
        index = des.index("Parent=")
        inter = des[index:]
        index = inter.index(';')
        trid = inter[7:index]  # transcript id
        # get gene name
        index = des.index("gene=")
        inter = des[index:]
        index = inter.index(';')
        gene = inter[5:index]
        # build a df
        row = pd.DataFrame({'Gene':[gene],'TrID':[trid]})
        df = df.append(row)
    # 3. output to file
    df.to_csv(outputFile,sep='\t',index=False)    
# merge_all_genes_with_diff_pr_from_refseq('','/data/shangzhong/VariantCall/Results/hamster_DNA_repair_pr_different_from_refseq.txt',
#                                          'hamster_DNA_repair_pr_different_from_refseq.fa',
#                                          '/data/shangzhong/VariantCall/Results')
