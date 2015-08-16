import subprocess,os
from Bio import SeqIO,Seq
import pandas as pd
def extract_geneID(outputFile,gff3):
    """
    This function extracts gene ids from gff3 annotation file
    """
    res = open(gff3,'r')
    gene_list = []
    output = open(outputFile,'w')
    for line in res:
        try:
            index = line.index('GeneID:')
            index = index + 7
        except:
            continue
        geneid = ''
        while line[index] !=';' and line[index] != ',':
                geneid = geneid + line[index]
                index = index + 1
                
        gene_list.append(geneid)
    gene_list = list(set(gene_list))    
    for item in gene_list:
        output.write(item + '\n')
    output.close()


def extract_gff4cufflink_symbol2ID(outputFile,gffFile,source,transcript='no',protein='no'):
    """
    This function extracts gene id,gene symbol,chromosome from annotation file and stores these inoformation
    into a file with 3 colums. Can also store 4th column for transcript accession or protein accession.
    Notice that only one of transcript and protein can be set to yes.
    
    * outputFile: outputFile name
    * gffFile: gffFile name from ncbi or ensembl
    * source: str. annotation file source. eg: 'ncbi','ensembl'
    return a file with 3 columns: [id,symbol,chromosome] or 4 columsn: [id,symbol,chromosome,trid/prid]
    """
    handle = open(gffFile,'r')
    #columns = ['GeneID','GeneSymbol','GeneAccession']
    f = open('inter.txt','w')
    #f.write('\t'.join(columns)+'\n')
    for line in handle:
        if line.startswith('#'): continue
        item = line.split('\t')
        chm = item[0]  # chromosome
        des = item[-1] # description
        if source == 'ncbi':
            geneID = 'GeneID:'
            geneName = 'gene='
            trAc = 'transcript_id='
            prAc = 'protein_id='
        else:
            geneID = 'gene_id='
            geneName = 'gene_name='
        
        if transcript == 'yes':
            if item[2] != 'mRNA': continue
            # get transcript accession number
            try:
                index = des.index(trAc)
                index = index + len(trAc)
                trAc = ''
                while des[index] !=';' and des[index] != ',' and des[index]!='\n':
                    trAc = trAc + des[index]
                    index = index + 1
            except:
                continue
            # get transcript id
            try:
                index = des.index('ID=')
                index = index + len('ID=')
                trID = ''
                while des[index] !=';' and des[index] != ',' and des[index]!='\n':
                    trID = trID + des[index]
                    index = index + 1
            except:
                continue
            outline = [trID,trAc]
        if protein == 'yes':
            if item[2] != 'CDS': continue
        # get protein accession number
            try:
                index = des.index(prAc)
                index = index + len(prAc)
                prid = ''
                while des[index] !=';' and des[index] != ',' and des[index]!='\n':
                    prid = prid + des[index]
                    index = index + 1
            except:
                continue
        # get transcript id
            try:
                index = des.index('Parent=')
                index = index + len('Parent=')
                trID = ''
                while des[index] !=';' and des[index] != ',' and des[index]!='\n':
                    trID = trID + des[index]
                    index = index + 1
            except:
                continue
            outline = [trID,prid]
        # get gene id
        try:
            index = des.index(geneID)
            index = index + len(geneID)
            ids = ''
            while des[index] !=';' and des[index] != ',':
                    ids = ids + des[index]
                    index = index + 1
        except:
            continue
        # get gene symbol
        try:
            index = des.index(geneName)
            index = index + len(geneName)
            symbol = ''
            while des[index] !=';' and des[index] != ',' and des[index]!='\n':
                symbol = symbol + des[index]
                index = index + 1
        except:
            continue
        outline = [ids,symbol,chm] + outline
        f.write('\t'.join(outline) + '\n')
    f.close()
    cmd = ('sort inter.txt | uniq > {output}').format(output=outputFile)
    subprocess.call(cmd,shell=True)
    os.remove('inter.txt')
    return outputFile


def extractAllGffIDs(outputFile,gffFile,source):
    """
    This function extracts all the ids in gff file and output all result into 5 column file.
    [gene_id,gene_symbol,chromosome,rna_accession,protein_accession]
    
    * gffFile: str. annotation
    """
    file1 = extract_gff4cufflink_symbol2ID(outputFile,gffFile,source,transcript='yes',protein='no')
    df1 = pd.read_csv(file1,header=None,sep='\t',names=['GeneID','GeneSymbol','Chrom','TrID','TrAccess'])
    file2 = extract_gff4cufflink_symbol2ID(outputFile,gffFile,source,transcript='no',protein='yes')
    df2 = pd.read_csv(file2,header=None,sep='\t',names=['GeneID','GeneSymbol','Chrom','TrID','PrAccess'])
    df = pd.merge(df1,df2,how='outer',on=['GeneID','GeneSymbol','Chrom','TrID'])
    df.fillna('-',inplace=True)
    df[['GeneID','GeneSymbol','Chrom','TrAccess','PrAccess']].to_csv(outputFile,sep='\t',index=False)
# gff3 ='/data/shangzhong/RibosomeProfiling/Database/combined.gff'
# extractAllGffIDs('/data/shangzhong/RibosomeProfiling/combined_AllIDS.txt',gff3,'ncbi')

def removeIDVersion(ID):
    """
    This function remove version information in all kinds of id versions
    
    * ID: str. should be genome mRNA or protein accession number with version information. eg: NW0123.1
    """
    if '.' in ID:
        return ID[:ID.index('.')]
    else:
        return ID
        
def extractAllPr(gffFile,feature='CDS',out=''):
    """
    This function converts gff file to file for a list of genes in the input. Keep different proteins.
    Only extract exons or CDS of genes.
    
    * gffFile: str. Filename of gff file
    * genes: list. A list of gene IDs.
    """
    gff_handle = open(gffFile,'r')
    if out =='':
        out = gffFile[:-4] + '.txt'
    out_handle = open(out,'w')
    for line in gff_handle:
        if line.startswith('#'):
            continue
        item = line[:-1].split('\t')
        if item[2] != feature:
            continue
        try:
        # get protein id
            prid = ''
            index = line.index('protein_id=')
            index = index + 11
            while line[index] not in ['\n',';']:
                prid = prid + line[index]
                index = index + 1
        except:
            print line, 'dont have prid'
            continue
        try:
            # get gene id
            geneid = ''
            index = line.index('GeneID:')
            index = index + 7
            while line[index] !=';' and line[index] != ',':
                    geneid = geneid + line[index]
                    index = index + 1
        except:
            geneid = '-'
            print line, 'dont have geneid'
        outline = '\t'.join([item[0],str(int(item[3])),item[4],geneid,prid,item[6]]) + '\n'
        out_handle.write(outline)
    gff_handle.close()
    out_handle.close()
    df = pd.read_csv(out,sep='\t',header=None,names=['Chr','cds_start','cds_end','GeneID','Pr_Access','Strand'])
    df = df.drop_duplicates()
    df = df.sort(['Pr_Access','cds_start'])
    #df['Chr'] = df['Chr'].astype(str).apply(lambda x: removeIDVersion(x))
    #df['Pr_Access'] = df['Pr_Access'].astype(str).apply(lambda x: removeIDVersion(x))
    df['GeneID'] = df['GeneID'].astype(str)
    df.to_csv(out,index=False,sep='\t')
# gffFile = '/data/shangzhong/RibosomeProfiling/Database/combined.gff'
# extractAllPr(gffFile,feature='CDS',out='/data/shangzhong/RibosomeProfiling/cho_pr/10_pr_cds.txt')

#===============================================================================
#                      some basic sequence process
#===============================================================================

def vari3letter2vari1letter(vari):
    """
    This function changes variants in the form of Arg359Lys, to one letter format
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
    sequence = sequence.upper()
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
    # start decode
    i = 0
    code = RNA[i:i+3]
    try:
        amino = codon[code]
    except:
        amino = 'X'
    AA = amino  # get the 1st amino acids
    while amino != '':
        i = i + 3
        if i < len(RNA) - 2:
            try:
                code = RNA[i:i+3]
                amino = codon[code]
                AA = AA + amino            
            except:
                amino = 'X'
                AA = AA + amino
        else:
            break
    return AA


def get_AA_from_gff(refDNA_dic,gffFile,chrom,gene,transcript,outputFile=''):
    """
    This function extracts amino acid sequence from gff annotation file, given 
    reference sequence dictionary, chromosome name, gene name, transcript id
    and CodonFile. And then compare the amino acid seqeunce to online protein sequence.
    Will output 3 files, one stores inconsistant prs in gff file, one stores inconsistant
    prs online, one stores gene--transcript--protein ids.
    
    * refDNA_dic: dictionary get from command record_dict = SeqIO.index(fastaFile,'fasta')
    * gffFile: gff annotation file
    * chrom: chromosome name
    * gene: gene symbol
    * transcript: transcript id; eg: rna1234
    * outputFile: file that stores online refseq pr sequence that is different from those in annotation file
    
    """
    #-------- (1) get chromosome sequence ------
    try:
        chr_seq = str(refDNA_dic[chrom].seq)
    except:
        print chrom,'does not exit in fa file'
        return ''
    # get CDS, the coding sequence, different genes may have same transcript id, so here we need to grep gene after grep trid.
    cmd = ('grep \'Parent={trans_id}\' {gff} | grep \'{gene}\' | awk -F \'\\t\' \'{ifcmd} {printrow}\'').format(trans_id=transcript,gff=gffFile,gene=gene,
                                                                                          ifcmd='$3==\"CDS\"',printrow='{print $4FS$5FS$7FS$9}') 
    exon_loci = subprocess.check_output(cmd,shell=True)
    exon_loci_list = exon_loci[:-1].split('\n')
    
    #-------- (2) combine CDS together ---------
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
            transcript_seq = transcript_seq + chr_seq[i-1:j].upper()
        elif k == '-':  # this transcript is in negtive strand
            transcript_seq = transcript_seq + DNA_complement_reverse(chr_seq[i-1:j])
    #------- (3) get the amino acid sequence in gff file ---------- 
    AA = Seq.translate(transcript_seq)
    if AA.endswith('*'):
        AA = AA[:-1]
    #AA = RNA2Pr(transcript_seq,CodonFile)
    
    #------- (4) get the online protein sequence-----
    try:
        index = exon_loci_list[0].index('protein_id=')
        protein_id = exon_loci_list[0][index + 11:]
        if ';' in protein_id:
            index = protein_id.index(';')
            protein_id = protein_id[:index]
        from Bio import Entrez
        Entrez.email = 'shl198@eng.ucsd.edu'
        handle = Entrez.efetch(db = 'protein',id=protein_id,rettype='fasta',retmode='text')
        record = handle.read()
        sequence = ''.join(record.split('\n')[1:-2])
        if AA != sequence:
            if outputFile=='':
                filename = 'LocalOnlineProteinDiff.fa'
            else:
                filename = outputFile    
            with open(filename,'a') as f:
                f.write('>' + exon_loci_list[0].split('\t')[-1] + '\n' + sequence + '\n\n')
            with open(filename[:-2]+'gff.fa','a') as f_gff:
                f_gff.write('>'+gene+'_'+transcript+'_'+protein_id+'\n'+AA+'\n\n')
            with open(filename[:-2]+'id.fa','a') as f_id:
                f_id.write(gene+'\t'+transcript+'\t'+protein_id+'\n')
    except:
        print transcript,'do not have protein id in annotation file'
        
    return AA

# eg:
# record_dict = SeqIO.index('/opt/genome/cho/chok1.fa','fasta')
# get_AA_from_gff(record_dict,'/opt/genome/cho/chok1.gff3','NW_003615779.1','LOC103160018','gene23081',
#                 '/data/shangzhong/VariantCall/RNA_codon_table.txt')


def get_diff_pr_from_refseq(outputFile,ref_fa,ref_gff3,gene_list_file=''):
    """
    This function tests whether protein sequences encoded by the genes in the given list
    are different from refseq protein sequnce. Return inconsistant sequences
    
    * outputFile: file that stores the inconsistant sequences.
    * ref_fa: str. Reference fasta file.
    * ref_gffs: str. Reference annotation file 
    * gene_list: str. A list of gene names. If not provided, it will test all genes in gff file
    """
    # 1. generate gene list
    if gene_list_file != '':
        df = pd.read_csv(gene_list_file,header=0,names=['Gene'])
        genes = df['Gene'].tolist()
    else:
        # without gene list, it will process all the genes in annotation file
        genes = []
        with open(ref_gff3,'r') as gff3:
            for line in gff3:
                try:
                    index = line.index('gene=')
                    inter = line[index:]
                    if ';' in inter:
                        index = inter.index[';']
                        gene = inter[5:index]
                    else:
                        gene = inter[5:-1]
                except:
                    continue
                if gene not in genes:
                    genes.append(gene)
            genes.sort()
    print 'there are',len(genes),'genes in',ref_gff3
    # define three list
    query_gene=[];query_chrom=[];query_trid=[]
    # get all genes, chromosomes and trids
    for gene in genes:
        cmd = ('grep \'gene={gene};\' {gff3} | awk -F \'\\t\' \'{ifcmd} {printrow}\'').format(
                gene=gene,gff3=ref_gff3,ifcmd='$3==\"CDS\"',printrow='{print $1FS$9}')
        record = subprocess.check_output(cmd,shell=True)
        # no CDS exist
        if record == '':
            print gene,'does not have coding sequence'
            continue
        records = record[:-1].split('\n')
        for line in records:
            item = line.split('\t')
            # get trid
            index1=item[1].index('Parent=')
            index2=item[1].index('Dbxref=')
            trid=item[1][index1+7:index2-1]
            if (item[0] in query_chrom) and (trid in query_trid):
                continue
            else:
                query_gene.append(gene);query_chrom.append(item[0]);query_trid.append(trid)
    # loop for each trascript
    refDNA_dic = SeqIO.index(ref_fa,'fasta')
    for chrom,gene,trid in zip(query_chrom,query_gene,query_trid):
        print chrom,gene,trid,'start analysis'
        try:
            AA = get_AA_from_gff(refDNA_dic,ref_gff3,chrom,gene,trid,outputFile)
            if AA == '':
                print chrom,gene,trid,'not in reference fa file'
            else:
                print chrom,gene,trid,'finished'
        except:
            print 'fail to get AA and compare'
            raise
        
# outputFile = '/data/shangzhong/RibosomeProfiling/chok1_pr_different_from_refseq.fa'
#genelist = '/data/shangzhong/VariantCall/CHOk1/genelist.txt'
# get_diff_pr_from_refseq(outputFile,#'/opt/genome/cho/chok1.fa','/opt/genome/cho/chok1.gff3')
#                         '/data/shangzhong/RibosomeProfiling/Database/combined.fa',
#                         '/data/shangzhong/RibosomeProfiling/Database/combined.gff')

# outputFile = '/data/shangzhong/VariantCall/Hamster/hamster_pr_different_from_refseq.fa'
# genelist = '/data/shangzhong/VariantCall/Hamster/genelist.txt'
# get_diff_pr_from_refseq(outputFile,
#                         '/opt/genome/hamster/hamster.fa',
#                         '/opt/genome/hamster/hamster.gff3',genelist)
# print 'job done'

