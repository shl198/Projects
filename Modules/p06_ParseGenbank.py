from BCBio import GFF
from Bio import SeqIO

def genbank2gff(inputfile):
    """
    This function transfer genbank file to gff file
    
    * inputfile: str. File name ends with gbk or gb
    
    return filename.gff
    """
    handle = open(inputfile,'r')
    if inputfile.endswith('gb'):
        out = inputfile[:-3] + '.gff'
    else:
        out = inputfile[:-4] + '.gff'
    out_handle = open(out,'w')
    result = SeqIO.parse(handle,'genbank')
    GFF.write(result,out_handle)
    
    handle.close()
    out_handle.close()
    return out

#genbank2gff('/data/shangzhong/RibosomeProfiling/Database/rRNAcho.gb')