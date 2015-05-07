# this file does statistic analysis for the results
# from collections import Counter
from Bio.Blast import NCBIXML
from Bio import SeqIO
import os
from combine_files import *
#============analyse the output sam file : notchok1_virus.final.sort.sam==========
#=================================================================================
def sam_result(filename):
    result = file(filename).readlines()
    data = [] # stores each line of result into tab separated
    for i in range(len(result)):
        data.append(result[i].split("\t"))
    ref = [] # stores duplicated reference genome names in sam file
    #-----------------list unique reference genomes------------------------------
    for i in range(len(data)):
        ref.append(data[i][2])
    num_ref = list(set(ref))      # unique reference names
    #---------------calculate duplication of each genome-------------------------
#     dupli = [] # stores how many reads mapping to unique reference genome
#     stat = Counter(ref) # generate a dictionary count the number
#     for reference in num_ref:
#         dupli.append(stat[reference])
    #-----------calculate how many reads map to reference uniquely---------------
    reads = []  # stores number of reads mapping to each of num_ref
    uniq_read = [] # stores number of uniq reads mapping to each of num_ref
    uniq_seq = [] # store number of uniq seq mapping to each of num_ref
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
    return [num_ref,[reads,uniq_read,uniq_seq]]
##==============================================================================
##=====analyse the blast output tabular file: neither_chok1_virus1.txt======
##==============================================================================
def single_blast_result(filename,orig_file):
    #filename is blast tabular file
    #orig_file is neitherchok1_virus.fa
    result = file(filename).readlines() 
    data = []   # stores each tab separated line of result
    for i in range(len(result)):
        data.append(result[i].split("\t"))
    ref = []    # stores all duplicated reference
    #-------------------------unique reference genome----------------------------
    for i in range(len(data)):
        ref.append(data[i][1])
    num_ref = list(set(ref))        # stores unique reference
    #--------------------calculate duplication of each genome--------------------
#     dupli = []
#     stat = Counter(ref)
#     for reference in num_ref:
#         dupli.append(stat[reference])
    #---------------calculate how many reads map to reference uniquely-----------
    reads = []
    reads_name = []
    uniq_read = [] # store number of uniq reads mapping to each of num_ref
    uniq_read_name = []
    for i in range(len(num_ref)):   
        inter_read = []
        for j in range(len(ref)):
            if num_ref[i] == ref[j]:
                inter_read.append(data[j][0])
        reads.append(len(inter_read))
        reads_name.append(inter_read)      
        uniq_read.append(len(list(set(inter_read))))
        uniq_read_name.append(list(set(inter_read)))
    #--------------find number of unique sequences mapping to each reference----------------
    result = SeqIO.parse(open(orig_file,"rU"),"fasta")
    orig_read = []
    orig_seq = []
    uniq_seq = []
    #-------------store orig_file sequence and read name into two variables
    for item in result:
        orig_read.append(item.id)
        orig_seq.append(str(item.seq))
    
    for i in range(len(num_ref)):
        inter_seq = []
        for j in range(len(uniq_read_name[i])):
            index = orig_read.index(uniq_read_name[i][j])
            inter_seq.append(orig_seq[index])
        uniq_seq.append(len(list(set(inter_seq))))
    return [num_ref,[reads,uniq_read,uniq_seq]]
#------------analyze two blast result----------------------------------------
# def pair_blast_result(path,file1,file2,origion1,origion2):
#     #file are two blast tablar result
#     #origion are two neitherchok1_virus.fa 
#     combine_files(path + '/output.txt', file1,file2)
#     combine_files(path + '/origion.fa',origion1,origion2)
#     [num_ref,[dupli,uniq_read,uniq_seq]] = single_blast_result(path + '/output.txt', path + '/origion.fa')
#     os.remove(path + '/origion.fa')
#     os.remove(path + '/output.txt')
#     return [num_ref,[dupli,uniq_read,uniq_seq]]    
#============================================================================
#===============combine sam result and blastresult===========================
def combine_results(fir_ref,fir_stat,sec_ref,sec_stat):
    # first_ref = a list of reference names. first_dupli = the namber of 
    # reads mapped to each of the reference. They have same length.    
    ref = fir_ref
    sam = fir_stat
    blast_dupli = sec_stat[0]
    blast_read =sec_stat[1]
    blast_seq = sec_stat[2]
    dupli = [0]*len(ref)         # these three stores the statistic results
    read = [0] * len(ref)
    sequence = [0] * len(ref)
    
    for reference in sec_ref:
        if reference in fir_ref:
            index = fir_ref.index(reference)  # index of blast reference in sam reference
            n = sec_ref.index(reference)  # index of blast reference in blast reference        
            dupli[index] = blast_dupli[n]
            read[index] = blast_read[n]
            sequence[index] = blast_seq[n]
        else:
            ref.append(reference)
            for i in range(len(sam)):
                sam[i].append(0)
            index = sec_ref.index(reference)
            dupli.append(blast_dupli[index])
            read.append(blast_read[index])
            sequence.append(blast_seq[index])
    sam.append(dupli); sam.append(read);sam.append(sequence)
    sam.insert(0,ref)
    return sam

##========================get final results============================================
##=====================================================================================
def get_final_results(output,*paths):
    #this function combines all the results into one variable
    sam_file = 'notchok1_virus.final.sort.sam'
    blast_file1 = 'neither_chok1_virus1.txt'
    blast_file2 = 'neither_chok1_virus2.txt'
    orig1 = 'neither_chok1_virus1.fa'
    orig2 = 'neither_chok1_virus2.fa'
    [sam_ref,sam_stat] = sam_result(paths[0] + '/' + sam_file)
    print 'sam result done'
    [blast_ref,blast_stat] = pair_blast_result(paths[0],paths[0] + '/' + blast_file1,paths[0] + '/' + blast_file2, paths[0]+'/'+ orig1,paths[0]+'/'+orig2)
    print 'blast result done'
    result = combine_results(sam_ref,sam_stat,blast_ref,blast_stat)    
    for i in range(1,len(paths)):
        [sam_ref,sam_stat] = sam_result(paths[i] + '/' + sam_file)
        print 'sam result done'
        result = combine_results(result[0],result[1:],sam_ref,sam_stat)
        print 'sam combine done'
        [blast_ref,blast_stat] = pair_blast_result(paths[i],paths[i] + '/' + blast_file1,paths[i] + '/' + blast_file2, paths[i]+'/'+ orig1,paths[i]+'/'+orig2)
        print 'blast result done'
        result = combine_results(result[0],result[1:],blast_ref,blast_stat)
        print 'blast combine done'    
    # get name
    integrate_name_type(result,0)
    # get results for simplified version
    simp_result = result[:3]
    n = len(result) 
    iteration = (n - 3)/6
    for i in range(iteration):
        simp_dupli = [result[i*6 + 3 ][j] + result[i*6 + 6][j] for j in range(len(result[0]))]
        simp_read = [result[i*6 + 4][j] + result[i*6 + 7][j] for j in range(len(result[0]))]
        simp_seq = [result[i*6 + 5][j] + result[i*6 +8][j] for j in range(len(result[0]))]
        simp_result.append(simp_dupli)
        simp_result.append(simp_read)
        simp_result.append(simp_seq)
    write2file(output,result,simp_result)
    return [result,simp_result]
#======================================================================================
##====================================================================================


#------------get gene name and type---------------------------------
from Bio import Entrez
Entrez.email = "shangzhong0619@gmail.com"
def get_name(refname,datatype):
    index = refname.replace('|','x',1).index('|')
    ncid = refname[3:index]   
    handle = Entrez.efetch(db=datatype,id=ncid,rettype='gb')
    record = handle.read()
    gene_result = record.split()
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
    return [gene_name,gene_type]
# #============================================================================
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
#------output the result------------
def write2file(output,result,simp_result):
    stat_result = open(output + '/' + "stat_result.txt",'w')
    for i in range(len(result[0])):
        for j in range(len(result)):    
            stat_result.write(str(result[j][i])+"\t")
        stat_result.write("\n")
    stat_result.close()
    stat_result = open(output + '/' + "simp_result.txt",'w')
    for i in range(len(simp_result[0])):
        for j in range(len(simp_result)):    
            stat_result.write(str(simp_result[j][i])+"\t")
        stat_result.write("\n")
    stat_result.close()
    
def write_file(outputfile,result):
    stat_result = open(outputfile,'w')
    for i in range(len(result[0])):
        for j in range(len(result)):    
            stat_result.write(str('\t' + result[j][i]))
        
    stat_result.close()
#============================================================================
#             get stat of reads used for assembly                   =========
#============================================================================
# the sam result would be the same, the blast result is different.
def stat_assem(sam1,sam2,blast1,blast2):
    # these two files are final blast paired reads. blast_reads1.fa, blast_reads2.fa
    sam_stat1 = list(SeqIO.parse(open(sam1,"rU"),"fasta"))
    sam_stat2 = list(SeqIO.parse(open(sam1,"rU"),"fasta"))
    sam_sequence = []
    sam_reads = len(sam_stat1)   # number of uniq pairs
    for item in sam_stat1:
        sam_sequence.append(str(item.seq))
    for item in sam_stat2:
        sam_sequence.append(str(item.seq))
    sam_uniq_seq = len(list(set(sam_sequence)))  # number of uniq sequence
    # blast results    
    blast_stat1 = list(SeqIO.parse(open(blast1,"rU"),"fasta"))
    blast_stat2 = list(SeqIO.parse(open(blast2,"rU"),"fasta"))
    blast_seq = []
    blast_reads = len(blast_stat1)
    for item in blast_stat1:
        blast_seq.append(item.seq)
    for item in blast_stat2:
        blast_seq.append(item.seq)
    blast_uniq_seq = len(list(set(blast_seq)))
    uniq_seq = sam_uniq_seq + blast_uniq_seq
    uniq_reads = sam_reads + blast_reads
    return [uniq_reads, uniq_seq]
def stat4allassem(output,*paths):
    n = len(paths)
    sam1 = 'notchok1_virus1.fa'
    sam2 = 'notchok1_virus2.fa'
    blast1 = 'blast_reads1.fa'
    blast2 = 'blast_reads2.fa'    
    result = []
    for i in range(n):
        [uniq_reads,uniq_seq] = stat_assem(paths[i] + '/' + sam1, paths[i] + '/' + sam2, paths[i] + '/' + blast1, paths[i] + '/' + blast2)
        result.append(uniq_reads)
        result.append(uniq_seq)
    return result