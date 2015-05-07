from Bio import SeqIO
import os,subprocess
def change_ncbi_annotation_name(outputfile,inputfile,inter):
    """
    this function changes the reference name of fasta file to accession number
    eg: change '>gi|614415508|ref|NW_006834731.1| Cricetulus griseus unplaced genomic scaffold, 
    alternate assembly C_griseus_v1.0 C26700950, whole genome shotgun sequence' to 'NW_006834731.1'
    
    * outputfile: fa file with accession nubmer as sequence name
    * inputfile: fa file with only name
    * inter: a temp file that stores the interium information, can be any name
    """
    reference = SeqIO.parse(open(inputfile,'rU'),'fasta')
    output = open(inter,'w')
    for item in reference:
        name =item.id
        start = name.index('|ref|')
        end = name.rfind('|')
        item.id = name[start+5:end]
        SeqIO.write(item,output,'fasta')
    output.close()
    os.system("""cut -d ' ' -f 1 %s > %s"""  % (inter,outputfile))
    os.system('rm ' + inter)

# input = '/data/shangzhong/DetectVirus/Database/150420virus.fa'
# out = input[:-2]+'gff.fa'; inter= 'test.txt'
# change_ncbi_annotation_name(out,input,inter)

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

# extractRefseqPr('/data/shangzhong/CHO2Human/mRNA/rna68_human.fna','/data/shangzhong/CHO2Human/mRNA/rna68.fna','Homo sapiens')
# extractRefseqPr('/data/shangzhong/CHO2Human/mRNA/rna68_cho.fna','/data/shangzhong/CHO2Human/mRNA/rna68.fna','Cricetulus griseus')
# extractRefseqPr('/data/shangzhong/CHO2Human/mRNA/rna68_mouse.fna','/data/shangzhong/CHO2Human/mRNA/rna68.fna','Mus musculus')

def embl2fasta(in_embl,out_fa):
    """
    This function convert embl file format to fasta format
    
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