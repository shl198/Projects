from Bio import SeqIO
import os,subprocess
import pandas as pd
from p05_ParseGff import extractAllGffIDs
from Bio.SeqRecord import SeqRecord
from Bio import Seq

def chunk(l,n):
    n = max(1,n)
    result = [l[i:i+n] for i in range(0,len(l),n)]
    return result

def change_ncbi_annotation_name(inputFa,inter='inter.fa'):
    """
    this function changes the reference name of fasta file to accession number
    eg: change '>gi|614415508|ref|NW_006834731.1| Cricetulus griseus unplaced genomic scaffold, 
    alternate assembly C_griseus_v1.0 C26700950, whole genome shotgun sequence' to 'NW_006834731.1'
    
    * inputFa: fa file with only name
    * inter: a temp file that stores the interium information, can be any name
    return fa file with accession nubmer as sequence name
    """
    reference = SeqIO.parse(open(inputFa,'rU'),'fasta')
    output = open(inter,'w')
    for item in reference:
        name =item.id
        start = name.index('|ref|')
        end = name.rfind('|')
        item.id = name[start+5:end]
        SeqIO.write(item,output,'fasta')
    output.close()
    outputFile = inputFa[:-2] + 'accession.fa'
    cmd = ("cut -d \' \' -f 1 {inter} > {out}").format(inter=inter,out=outputFile)
    subprocess.call(cmd,shell=True)
    os.system('rm ' + inter)


def remove_duplicate(outputfile,inputfile):
    """
    this function removes the duplicated references in a fasta file
    
    * outputfile: fa file with unique sequence
    * inputfile: fa file
    """
    reference = SeqIO.parse(open(inputfile,'rU'),'fasta')
    output = open(outputfile,'w')
    ref_name = []
    for item in reference:
        if item.id in ref_name:
            continue
        else:
            ref_name.append(item.id)
            SeqIO.write(item,output,'fasta')
    output.close()


def change_ensembl_annotation_name(outputfile,inputfile):
    """
    this function changes the reference name of fasta file from ensemble to chromosome number
    
    * outputfile:    fasta output file name
    * inputfile:    fasta input file name
    """
    reference = SeqIO.parse(open(inputfile,'rU'),'fasta')
    output = open(outputfile,'w')
    for item in reference:
        name = item.description
        start = name.index('chromosome')
        end = name.index(',')
        item.id = name[start + len('chromosome') + 1:end]
        item.description = ''
        SeqIO.write(item,output,'fasta')
    output.close


def extractRefseqPr(outputFile,refProteinFile,organism):
    """
    This function extracts the proteins of specific organism
    from refseq protein file which include all proteins.
    
    * refProtein: the refseq fasta file proteins
    * organism: organism name, must be the same as ncbi defines
    """
    res = SeqIO.parse(open(refProteinFile,'rU'),'fasta')
    output = open(outputFile,'w')
    for item in res:
        name = item.description
        if organism in name:
            SeqIO.write(item,output,'fasta')
    output.close()


def embl2fasta(in_embl,out_fa):
    """
    This function converts embl file format to fasta format
    
    * in_embl: input embl format sequence file
    * out_fa: output fasta format sequence file
    """
    count = SeqIO.convert(in_embl, 'embl', out_fa, 'fasta')
    return count


def fq2fa(fqs):
    """
    This function transfer fq files to fa files
    Return a list of fa files. ['f1.fa','f2.fa']
    
    * fqs: a list of fq files, for paired end files should be [['fq1.fq.gz','fq2.fq.gz'],...]
                               for single end files should be [['fq1.fq.gz']...]
                               can also be ['fq1.fq.gz',...] or 'fq1.fq.gz'
    return a list of fa files
    """
    if isinstance(fqs,str):
        fa = fqs[:-5] + 'fa'
        cmd = ('gunzip -c {fq_gz} | sed -n \'1~4s/^@/>/p;2~4p\' > {fa}').format(fq_gz=fqs,fa=fa)
        subprocess.call(cmd,shell=True)
        return fa
    fas = []
    for fq in fqs:
        if isinstance(fq,list):
            fa_Files = fq2fa(fq)
            fas.extend(fa_Files)
        else:
            fa = fq[:-6] + '.fa'
            cmd = ('gunzip -c {fq_gz} | sed -n \'1~4s/^@/>/p;2~4p\' > {fa}').format(fq_gz=fq,fa=fa)
            
            subprocess.call(cmd,shell=True)
            # Compress the files
            cmd = ('gzip {f}').format(f=fa)
            subprocess.call(cmd.split(' '))
            fas.append(fa+'.gz')

    return fas

def splitFa(inputFile,maxItem):
    """
    This function split fa file into several files.
    
    * inputFile: str. Fasta file name
    * maxItem: int. Number of items in each files
    """
    n = 0;suffix=1
    fa = SeqIO.parse(open(inputFile,'rU'),'fasta')
    out = inputFile[:-3] + str(suffix) + '.fa'
    outHandle = open(out,'w')
    for record in fa:
        if n == maxItem:
            n = 0
            outHandle.close()
            suffix = suffix + 1
            out = inputFile[:-3] + str(suffix) + '.fa'
            outHandle = open(out,'w')
        else:
            SeqIO.write(record,outHandle,'fasta')
            n = n + 1
    outHandle.close()
    

def chrLength(fa):
    """
    This functions generates two column file, 1st is chromosome name,
    2nd is chromosome length
    
    * fa: str. Fasta file
    """
    len_dic = {}
    handle = open(fa,'r')
    chrs = SeqIO.parse(handle,'fasta')
    for chrome in chrs:
        len_dic[chrome.id] = len(chrome.seq)
    df = pd.DataFrame(len_dic.items(),columns=['Chr','Length'])
    out = fa[:-3] + '_ChrLen.txt'
    df.to_csv(out,sep='\t',index=False)
# fa= '/data/shangzhong/RibosomeProfiling/Database/combined.fa'
# chrLength(fa)

def sub_scaffold(fa,gffFile,genes,sub_fa):
    """
    This function extract scaffolds on which locate the genes provided.
    * fa: str. File name of reference fasta file
    * gffFile: str. Annotation file
    * genes: list. A list of gene symbols.
    * outIDFile: str. File name that stores all ids from gff file. 
    """
    # 1. read fa file
    ref_dict = SeqIO.index(fa,'fasta')
    # 2. get chromosome of all genes
    ID_file = 'inter_ID.txt'
    extractAllGffIDs(ID_file,gffFile,'ncbi')
    ID_df = pd.read_csv(ID_file,sep='\t',header=0) # names=['GeneID','GeneSymbol','Chrom','TrAccess','PrAccess']
    gene_df = ID_df[ID_df['GeneSymbol'].isin(genes)]
    chroms = list(set(gene_df['Chrom'].tolist()))
    os.remove(ID_file)
    # 3. output the chromosome sequence to file
    out_handle = open(sub_fa,'w')
    for chrom in chroms:
        chr_seq = ref_dict[chrom]
        SeqIO.write(chr_seq,out_handle,'fasta')


def divide_scaffold_by_len(faFile,n_part,out_path):
    """
    This function sort the genome by length and then separate them in to n groups. Each group have similar amound the total nucleotide
    * faFiles: str. Fasta file.
    * n_part: int. Separate to several groups.
    """
    n_part = int(n_part)
    reference = SeqIO.parse(open(faFile,'rU'),'fasta')
    chrom_len = {}
    for item in reference:
        chrom = item.id
        if chrom in chrom_len:
            assert False,'repeated reference'
        chrom_len[chrom] = len(item.seq)
    # transfer to pandas table
    chrom_len_df = pd.DataFrame(chrom_len.items(),columns=['chr','len'])
    chrom_len_df = chrom_len_df.sort('len')
    chrs = chrom_len_df['chr'].tolist()
    res = [[] for i in range(n_part)]
    for i in range(len(chrs)):
        res[i%n_part].append(chrs[i])
    if not os.path.exists(out_path): os.mkdir(out_path)
    for i in range(len(res)):
        handle = open(out_path + '/L'+str(i)+'.list','w')
        handle.write('\n'.join(res[i]))
        handle.close()

def stitch_scaffolds(fa,outFile,len_limit=200000000,dist=500):
    """
    This function merge multiple scaffold together to form a longer sequence. 
    * fa: str. Reference fa file name
    * outFile: str. Filename output to the file.
    * len_limit: int. Maximum length of each merged scaffold.
    * dist: int. Distance between each scaffold.
    """
    in_handle = open(fa,'r')
    out_handle = open(outFile,'w')
    sequence = ''
    n = 1
    for record in SeqIO.parse(in_handle,'fasta'):
        sequence += str(record.seq)
        if len(sequence) >= len_limit:
            item = SeqRecord(Seq(sequence), id = 'chr'+str(n),description="")
            SeqIO.write(item,out_handle,'fasta')
            sequence = ''
            n += 1
        else:
            sequence += 'N'*500
    if sequence != '':
        item = SeqRecord(Seq(sequence[:-500]), id = 'chr'+str(n),description="")
        SeqIO.write(item,out_handle,'fasta')
    # output the last one
    handle = open(outFile)
    for record in SeqIO.parse(handle,'fasta'):
        print len(record.seq)