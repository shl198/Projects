from Bio import SeqIO
import os 
def change_ncbi_annotation_name(outputfile,inputfile,inter):
    """this function can change the reference name of fasta file to accession number
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


def remove_duplicate(outputfile,inputfile):
    """this function can remove the duplicated references in a fasta file"""
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
    this function change the reference name of fasta file from ensemble to chromosome number
    outputfile:    fasta output file name
    inputfile:    fasta input file name
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

extractRefseqPr('/home/shangzhong/refseq/refseq68_human.faa','/home/shangzhong/refseq/ref.protein.faa','Homo sapiens')