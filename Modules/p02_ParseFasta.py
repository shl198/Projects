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
    outputfile:    the output file name
    inputfile:    input file name
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
inputfile = '/data/shangzhong/human/humanGenome/humangenome.fa'
outputfile = '/data/shangzhong/human/humanGenome/human.fa'
inter = '/opt/genome/hamster/inter.fa'
change_ensembl_annotation_name(outputfile,inputfile)