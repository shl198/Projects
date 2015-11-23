import pysam
# this file does statistic analysis for the results
# from collections import Counter
#============ analyse the samfile output =========
#=================================================================================
def sam_result(filename):
    """
    This function does some statistic calculation of the input sam file.
    Including unique reads, uniqure reference, # of reads mapping to reference.
    """
    # define output filename
    result = open(filename,'r')
    data = [] # stores each line of result into tab separated
    for item in result:
        data.append(item.split("\t"))
    ref = [] # stores duplicated reference genome names in sam file
    #-----------------list unique reference genomes------------------------------
    for i in range(len(data)):
        ref.append(data[i][2])
    num_ref = list(set(ref))      # unique reference names
    #-----------calculate how many reads map to reference uniquely---------------
    reads = []  # stores number of reads mapping to each of num_ref
    uniq_read = [] # stores number of uniq reads mapping to each of num_ref
    uniq_seq = [] # store number of uniq seq mapping to each of num_ref
    
    output = open(filename+'.uniq_seq.fa','w') # alter
    for i in range(len(num_ref)):
        inter = []    
        inter_read = []  # read names
        inter_seq = []
        for j in range(len(ref)):
            if num_ref[i] == ref[j]:
                inter.append(num_ref[i])
                inter_read.append(data[j][0])
                inter_seq.append(data[j][9])
        reads.append(len(inter_read))
        uniq_read.append(len(list(set(inter_read))))
        uniq_seq.append(len(list(set(inter_seq))))
        
        """ Line 25 and the next few lines output the unique sequence 
            into a fasta file.
        """
        uniqueSequence = list(set(inter_seq))
        for i in range(len(uniqueSequence)):
            output.write('>' + num_ref[i] +'_'+ str(i) + '\n')
            output.write(uniqueSequence[i] + '\n')
    output.close()  # alternate
    # output better in format [ref,stat] = the return value
    return [num_ref,[reads,uniq_read,uniq_seq]]
#======================================================================================
##====================================================================================
#------------get gene name and type---------------------------------
from Bio import Entrez,SeqIO
Entrez.email = "shl198@eng.ucsd.edu"
def get_name(refname,datatype):
    """
    this function can search full name given the accession number.
    Usually the first 2 columns of blast output, 3rd column of sam file.
    
    refname: eg. 'gi|625180359|ref|XP_007605693.1|'
    
    datatype: 'nucleotide' or 'protein'
    """
    
    index = refname.replace('|','x',1).index('|')
    ncid = refname[3:index]   
    handle = Entrez.efetch(db=datatype,id=ncid,rettype='gb')
    record = handle.read()
    handle.close()
    # try to find the record on line
    if len(record)>0:
        gene_result = record.split()
    else:
        gene_result = 'NA'
    if ("DEFINITION" in gene_result) & ("ACCESSION" in gene_result):
        defination = gene_result.index("DEFINITION")
        accession = gene_result.index("ACCESSION")
        gene_name = " ".join(gene_result[defination+1:accession])
        if datatype == 'nucleotide':
            gene_type = gene_result[4]
        else:
            if 'accession' in gene_result:
                nt_accession = gene_result.index('accession')
                ncid = gene_result[nt_accession + 1]
                handle = Entrez.efetch(db='nucleotide',id=ncid, rettype='gb')
                record = handle.read()
                gene_result = record.split()
                gene_type = gene_result[4]
            else:
                gene_type ='linear'
    else:
        gene_name = 'NA'
        gene_type = 'linear'
    return [gene_name,gene_type]
##============================================================================
##---integrate genome and type------------
def integrate_name_type(result,column,datatype):
    # result is the array that has the statistic resuls. column indicate which column of result stores the name information
    # each item in result represent a column in the matrix
    ref_name = []
    ref_type = []
    name = result[column]
    for i in range(len(name)):
        [gene_name,gene_type]=get_name(name[i],datatype)
        ref_name.append(gene_name)
        ref_type.append(gene_type)
    result.insert(1,ref_name)
    result.insert(1,ref_type)
    
def write_file(outputfile,result):
    stat_result = open(outputfile,'w')
    tresult = map(list,zip(*result))
    for line in tresult:
        inter = [str(f) for f in line]
        stat_result.write('\t'.join(inter) + '\n')
    stat_result.close()
#=============  pipeline for read and output sam result  =======
def outputSam(samfiles):
    for filename in samfiles:
        #----- (1) stats of sam ------
        [num_ref,stat]= sam_result(filename)
        sam = [num_ref,stat[0],stat[1],stat[2]]
        #----- (2) add annotation ----
        integrate_name_type(sam,0,'nucleotide')
        #----- (3) write to file -----
        outputfile = filename[:-3] + '.txt'
        write_file(outputfile,sam)
        print  filename + ' done'


class bam_parse(object):
    """
    This is for single end data analysis.
    Input is a bamHandle read by pysam.AlignmentFile(bamFile,'rb')
    """
    def __init__(self,bamHandle):
        self.bamHandle = bamHandle
        
    def bam_fwd_rev_count(self,max_len=37,seq_len=50):  # get fwd and reverse count data
        """
        Return two dictioanry: {name\tpos:count} for fwd and reverse alignment
        """
        count_dic = {}
        for align in self.bamHandle.fetch():
            if align.alen > seq_len:
                length = align.qlen   # query length: length of the nt in the read that maps
            else:
                length = align.alen   # alignment length
            if length > max_len: continue
            tid = align.reference_id
            ref_name = self.bamHandle.getrname(tid)
            align5 = align.pos + 1
            align3 = align.aend
            if align.is_reverse:
                key = '\t'.join([ref_name,str(align5),str(align3),'-',str(length)])
                if key in count_dic:
                    count_dic[key] = count_dic[key] + 1
                else:
                    count_dic[key] = 1
            else:  # forward strand
                key = '\t'.join([ref_name,str(align5),str(align3),'+',str(length)])
                if key in count_dic:
                    count_dic[key] = count_dic[key] + 1
                else:
                    count_dic[key] = 1
        return count_dic
    
    def align_len_distribution(self,seq_len=50):
        """
        This function returns a dictioanry of alignment length distribution. eg: {29:30%}
        """
        len_dist = {}
        for align in self.bamHandle.fetch():
            if align.alen > seq_len:
                length = align.qlen
            else:
                length = align.alen
            if length in len_dist:
                len_dist[length] = len_dist[length]+1
            else:
                len_dist[length] = 1
        total = self.bamHandle.mapped
        for key in len_dist:
            len_dist[key] = len_dist[key]/float(total)
        return len_dist
        
    def get_gene_cov(self,chrom,pos,end='end5'):
        """get covreage at gene position"""
        