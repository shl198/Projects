from Bio import SeqIO
import os 
def change_name(outputfile,inputfile,inter):
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

inputfile = '/opt/genome/hamster/cgr_alt_C_griseus_v1.0_chrUn.fa'
outputfile = '/opt/genome/hamster/hamster.fa'
inter = '/opt/genome/hamster/inter.fa'
change_name(outputfile,inputfile,inter)