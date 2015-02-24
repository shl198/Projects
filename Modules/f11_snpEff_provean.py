import subprocess,os
from Bio import SeqIO
import pandas as pd

#===============================================================================
#                      some basic sequence process
#===============================================================================

def vari3letter2vari1letter(vari):
    """
    This function change variants in the form of Arg359Lys, to one letter format
    return R359K
    
    * vari: input vairants
    """
    dic = {'Ala':'A', 'Arg':'R', 'Asn':'N', 'Asp':'D', 'Cys':'C',
           'Glu':'E', 'Gln':'Q', 'Gly':'G', 'His':'H', 'Ile':'I',
           'Leu':'L', 'Lys':'K', 'Met':'M', 'Phe':'F', 'Pro':'P',
           'Ser':'S', 'Thr':'T', 'Trp':'W', 'Tyr':'Y', 'Val':'V',
           'Ter':'*'}
    res = ''
    n = 0
    for i in range(len(vari)):
        if n == len(vari) - 1:
            res = res + vari[n]
            break
        else:
            if vari[i].isupper():
                res = res + dic[vari[n:n+3]]
                n = n + 3
            else:
                if i < n:
                    continue
                else:
                    res = res + vari[n]
                    n = n + 1
 
    return res

def DNA_complement_reverse(sequence):
    """
    This fucntion gets the complement reverse sequence of a sequence
    
    * sequence: a string of sequence in ATCG format
    """
    dic = {'A':'T','T':'A','C':'G','G':'C','N':'N'}
    result = ''
    for letter in sequence:
        result = dic[letter] + result
    return result

def RNA2Pr(RNA,CodonFile):
    """
    This function translate RNA to protein,
    return amino acid sequence
    
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
    res.close()
    # decode
    i = 0
    code = RNA[i:i+3]
    try:
        amino = codon[code]
    except:
        amino = 'X'
    AA = amino
    while amino != '':
        i = i + 3
        if i < len(RNA) - 2:
            try:
                code = RNA[i:i+3]
                amino = codon[code]
                AA = AA + amino            
            except:
                amino = amino + 'X'
                AA = AA + amino
        else:
            break
    return AA


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
    chr_seq = str(refDNA_dic[chrom].seq)
    # get CDS, the coding sequence
    cmd = ('grep Parent={trans_id} {gff} | awk -F \'\\t\' \'{ifcmd} {printrow}\'').format(trans_id=transcript,gff=gffFile,
                                                                                          ifcmd='$3==\"CDS\"',printrow='{print $4FS$5FS$7FS$9}')
    exon_loci = subprocess.check_output(cmd,shell=True)
    exon_loci_list = exon_loci[:-1].split('\n')
    # combine CDS together
    start_loci = []; stop_loci = [];strand = []; 
    for item in exon_loci_list:
        locis = item.split('\t')
        start_loci.append(locis[0])
        stop_loci.append(locis[1])
        strand.append(locis[2])
    transcript_seq = ''
    for i,j,k in zip(start_loci,stop_loci,strand):
        i = int(i); j = int(j)
        if k == '+': # this means the transcript is in positive strand
            transcript_seq = transcript_seq + chr_seq[i-1:j]
        elif k == '-':  # this transcript is in negtive strand
            transcript_seq = transcript_seq + DNA_complement_reverse(chr_seq[i-1:j])
            
    # get the amino acid sequence
    AA = RNA2Pr(transcript_seq,CodonFile)
    # get the online protein sequence
    index = exon_loci_list[0].index('protein_id=')
    protein_id = exon_loci_list[0][index + 11:]
    from Bio import Entrez
    Entrez.email = 'shl198@eng.ucsd.edu'
    handle = Entrez.efetch(db = 'protein',id=protein_id,rettype='fasta',retmode='text')
    record = handle.read()
    sequence = ''.join(record.split('\n')[1:-2])
    with open('LocalOnlineProteinDiff.fa','a') as f:
        if AA != sequence:
            f.write('>' + exon_loci_list[0].split('\t')[-1] + '\n' + sequence + '\n\n')
    return AA

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
#                 provean
#===============================================================================
def vcf2input4provean(filteredVCF,record_dict,gffFile,CodonFile):
    vcf_df = pd.read_csv(filteredVCF,sep='\t')
    if vcf_df.empty:
        return [[],[]]
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
                one_letter_vari = vari3letter2vari1letter(vari[2:])
                if impact == 'HIGH' or impact == 'MODERATE':
                    try:
                        pro_vari_dic[trid].append(one_letter_vari)
                    except:
                        pro_vari_dic[trid] = [one_letter_vari]
                    # add trid:gene to the dictionary
                    if trid not in trid_gene_dic:
                        trid_gene_dic[trid] = gene
                else:
                    continue
    
    protein_files = [] # stores input files for 
    variant_files = []
    for trid in pro_vari_dic:
        gene = trid_gene_dic[trid]
        protein_file = gene + '_' + trid + '.protein.fa'
        vari_file = gene + '_' + trid + '.variant.txt'
        
        protein_files.append(protein_file)
        variant_files.append(vari_file)
        
        protein_output = open(protein_file,'w')
        AA = get_AA_from_gff(record_dict,gffFile,chrom,trid,CodonFile)
        protein_output.write(('>{gene_trid} \n{AA}').format(gene_trid=gene+'_'+trid,AA=AA))
        protein_output.close()
        
        variant_output = open(vari_file,'w')
        for vari in pro_vari_dic[trid]:
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
    output.write('# VARIATION\tSCORE\n')
    for protein,vari in zip(proteinFiles,variantFiles):
        result = run_provean(provean,protein,vari,support_set_path,support_set,thread)
        lines = result[:-1].split('\n')
        for line in lines:
            if line.startswith('#') or line.startswith('['):
                continue
            elif line.startswith('No'):
                output.write(vari[:-12] + '\t' + 'No variantions entered' + '\n')
            else:
                item = line.split('\t')
                output.write(vari[:-12] + '_' + '\t'.join(item) + '\n')
        print protein,'provean analysis finishes\n'
    output.close()

# os.chdir('/data/shangzhong/VariantCall/CHOS_ToT_VS_chok1/filteredVcfResults/CHOS_TOT_A4')
# provean = '/home/shangzhong/Installation/provean-1.1.5/bin/provean.sh'
# proteinFiles = ['Nbn_rna17749.protein.fa']
# variantFiles = ['Nbn_rna17749.variant.txt']
# support_set_path = '/data/shangzhong/VariantCall/chok1_Provean_Supporting_Set'
# support_set = ['Nbn_rna17749.sss']
# outputFile = 'pro.txt'
# capture_provean_scores(outputFile,provean,proteinFiles,variantFiles,support_set_path,support_set,thread = 19)

def merge_all_provean_scores(previousFile='',*pathways):
    """
    This function combine all the provean results together in one table.
    if previousFile is provided, then the new result is appended to it.
    
    * previousFile: file that has provean results.
    * pathways: pathways that has folders for samples, with file named proveanScore.txt
    """
    columns = ['sample_name','gene_trID_vari','Score']
    if previousFile=='':
        df = pd.DataFrame(columns=columns)
    else:
        df = pd.read_csv(previousFile,sep='\t',header=None,
                         names=columns)
    for pathway in pathways:
        os.chdir(pathway)
        folders = [f for f in os.listdir(pathway) if os.path.isdir(os.path.join(pathway,f))]
        for folder in folders:
            f = folder + '/' + 'proveanScore.txt'
            data = pd.read_csv(f,sep='\t',header=0,names=['gene_trID_vari','Score'])
            sample = [folder]*len(data['Score'])
            data.insert(0,'sample_name',sample)
            df = df.append(data,ignore_index=True)
    
#merge_all_provean_scores('','/data/shangzhong/VariantCall/CHOS_ToT_VS_chok1filteredVcfResults')