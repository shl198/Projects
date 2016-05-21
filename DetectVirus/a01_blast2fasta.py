##==============transfer blast tabular file to fasta file====================
from Bio import SeqIO
import pandas as pd
import subprocess,os
#--------------read the 1st file -----------------
def blast2faPaired(blast1, blast2,origin1,origin2):
    """
    This function transfers blast tab demilinated txt files back to fasta files for paired end read
    
    * blast1: str. filename of blast result from mapping 1st read. eg: filename_1.blast.txt
    * blast2: str. filename of blast result from mapping 2nd read
    * origin1: str. filename of original 1st end fa file. eg: filename_1.fa
    * origin2: str. filename of original 2nd end fa file
    """
    df = pd.read_csv(blast1,sep='\t',names=['ReadName'],usecols=[0],header=None)
    df['ReadName'] = df['ReadName'].str[:-2].astype(str)
    name1 = list(set(df['ReadName'].tolist()))
    df = pd.read_csv(blast1,sep='\t',names=['ReadName'],usecols=[0],header=None)
    df['ReadName'] = df['ReadName'].str[:-2].astype(str)
    name2 = list(set(df['ReadName'].tolist()))
    # get the intersection name
    pair = list(set(name1).intersection(set(name2)))
    # read the original files
    handle1 = SeqIO.parse(open(origin1,'rU'),'fasta')
    handle2 = SeqIO.parse(open(origin2,'rU'),'fasta')
    # define output file name
    out1 = blast1[:-4] + '.fa';  out2=blast2[:-4] + '.fa'
    outhandle1 = open(out1,'w')
    for item in handle1:
        if item.id[:-2] in pair:
            SeqIO.write(item,outhandle1,'fasta')
    handle1.close()
    outhandle1.close()
    outhandle2 = open(out2,'w')
    for item in handle2:
        if item.id[:-2] in pair:
            SeqIO.write(item,outhandle2,'fasta')
    handle2.close()
    outhandle2.close()
    
    return [out1,out2]
#     # this command generates two files: blast_reads1.fa, blast_reads2.fa
#     #blast1 = 'blast_reads1.fa'
#     #blast2 = 'blast_reads2.fa'
#     blast1 = out1
#     blast2 = out2
#     # outputpath is in form: /media/shl/CHOS/A4
#     # origion files have reads neither mapped to chok1 nor virus after gsnap 
#     #alignment. File format is fasta.
#     # filename files have reads that map to virus, not chok1 after blast map origion
#     # files to virus. File format is tabular txt file.
#     #result1 = file("/home/shl/analyze_result/CHOS/Tot_CHOS_A4/neither_chok1_virus1.txt").readlines()
#     result1 = open(filename1).readlines()
#     data1 = []
#     for line in result1:
#         data1.append(line.split("\t")[0])
#     uniq_data1 = list(set(data1))
#     uniq_data1.sort()  # list all unique reads.
#     
#     #-------------read the 2nd file-----------------
#     result2 = open(filename2).readlines()   
#     data2 = []
#     for line in result2:
#         data2.append(line.split("\t")[0])
#     uniq_data2 = list(set(data2))
#     uniq_data2.sort()
#     
#     #-------------read the original file 1st and 2nd (neither_chok1_virus)-------
#     # In the file, it has reads not mapping to chok1 nor virus
#     records1 = SeqIO.parse(open(origion1,"rU"),"fasta")
#     records2 = SeqIO.parse(open(origion2,"rU"),"fasta")
#     #-------------find the matched reads and output-----------
#     normalized_name1 = []    # these two lists will store the reads names with out /1 or /2
#     normalized_name2 = []
#     for i in uniq_data1:
#         normalized_name1.append(i[:-2])
#     for i in uniq_data2:
#         normalized_name2.append(i[:-2])    # remove last two characters /1 and /2  
#     pair = list(set(normalized_name1).intersection(set(normalized_name2))) # find the paired reads in two files
#     pair1 = []
#     pair2 = []
#     for i in pair:
#         pair1.append(i + "/1")
#         pair2.append(i + "/2")
#     #------------output file into two fasta files--------------------
#     output_handle1 = open(outputpath + '/' + blast1,"w")
#     for item in records1:
#         if item.id in pair1:
#             SeqIO.write(item, output_handle1,"fasta") # write to file.
#     output_handle1.close()        
#     output_handle2 = open(outputpath + '/' + blast2,"w")
#     for item in records2:
#         if item.id in pair2:
#             SeqIO.write(item, output_handle2,"fasta") # write to file.
#     output_handle2.close()


def blast2faSingle(blastResultFiles,OriginFaFiles):
    """
    This function transfers blast tab demilinated txt files back to fasta files for single end read
    
    * blastResultFiles: list.  a list of blast tab demilited (format 7) files. eg: ['f1.txt','f2.txt',...]
    * OriginFaFiles:    list.  a list of fasta files that are queiries when doing blast. eg: ['f1.fa','f2.fa',...]
    return a fasta file for reads in balst results
    """
    # loop for each blast and origin
    res = []
    for blast,origin in zip(blastResultFiles,OriginFaFiles):
        df = pd.read_csv(blast,header=None,sep='\t',usecols=[0],names=['readName'])
        readName = df['readName'].tolist()
        handle = SeqIO.parse(open(origin,'rU'),'fasta') # read the origin file
        outFile = blast[:-3] + 'fa'
        outputHandle = open(outFile,'w')
        for item in handle:
            if item.id in readName:
                SeqIO.write(item,outputHandle,'fasta')
        outputHandle.close()
        res.append(outFile)
    
    return res


def blast2fa(blastResultFiles,OriginFaFiles,seqType):
    """
    This function combines the single end and pair end blast2fa file functions
    
    * blastResultFiles: list.  a list of blast result txt file in format 6 (tab delimited). eg: ['f1.txt','f2.txt']
    * OriginFaFiles:    list.  a list of fa.gz files that are used as queries in blastn. eg: ['f1.fa.gz','f2.fa.gz']
    
    return a list of fa files. eg: ['f1.blast.fa']
    """
    if OriginFaFiles[0].endswith('.gz'):
        for f in OriginFaFiles:
            cmd = ('gunzip {file}').format(file=f)
            subprocess.call(cmd.split(' '))
        OriginFaFiles = [f[:-3] for f in OriginFaFiles]
    #=== (1) deal with the single end samples
    if seqType == 'single':
        res = blast2faSingle(blastResultFiles,OriginFaFiles)
        # compress the results
#         for f in res:
#             cmd = ('gzip {file}').format(file=f)
#             subprocess.call(cmd.split(' '))
#         res = [f+'.gz' for f in res]
    #=== (2) deal with the paired end samples
    if seqType == 'pair':
        res = []
        blast_left=[f for f in blastResultFiles if f.endswith('1.blast.txt')]
        blast_right=[f for f in blastResultFiles if f.endswith('2.blast.txt')]
        origin_left=[f for f in OriginFaFiles if f.endswith('1.fa')] 
        origin_right=[f for f in OriginFaFiles if f.endswith('2.fa')]
        for blast1,blast2,origin1,origin2 in zip(blast_left,blast_right,origin_left,origin_right):
            pair = blast2faPaired(blast1, blast2,origin1,origin2)
            res.extend(pair)
    for f in OriginFaFiles: os.remove(f)
    return res


def merge_fa(seqType,*files):
    """
    This function merges fq or fa files in either single end or paired end
    
    * seqType: str. paired or single
    * files: list.  can be any number of list of files
    
    return a list with length of 1 or 2
    """
    res = []
    final_list = []
    for f in files:
            final_list.extend(f)
    # unzip the gzipped files
    for i in range(len(final_list)):
        if final_list[i].endswith('.gz'):
            cmd = ('gunzip {f}').format(f=final_list[i])
            print cmd
            subprocess.call(cmd.split(' '))
            final_list[i] = final_list[i][:-3]
    # deal with the single end samples
    if seqType == 'single':
        outFile = 'merged4Trinity.fa'
        res.append(outFile)
        cmd = ('cat {f_list} > {out}').format(f_list=' '.join(final_list),out=outFile)
        print cmd
        subprocess.call(cmd,shell=True)
    # deal with paired end samples
    if seqType == 'pair':
        left = [f for f in final_list if '_1.' in f]
        right = [f for f in final_list if '_2.' in f]
        outFile1 = 'merged4Trinity1.fa'
        cmd = ('cat {f_list} > {out}').format(f_list=' '.join(left),out=outFile1)
        print cmd
        subprocess.call(cmd,shell=True)
        outFile2 = 'merged4Trinity2.fa'
        cmd = ('cat {f_list} > {out}').format(f_list=' '.join(right),out=outFile2)
        print cmd
        subprocess.call(cmd,shell=True)
        res = [outFile1,outFile2]
    return res