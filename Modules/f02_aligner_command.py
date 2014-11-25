import subprocess
import os
from os import listdir
#==========  gsnap alignment  ==========================
def gsnap(fastqFiles,db_path, db_name,annotation,thread=1):
    """
    This function run gsnap command to map. Fastq files is a list, each item
    includes either paired fastq files or single ended files.
    db_path is the pathway. db_name is the prefix of database
    """
    map_result = []
    cmd = ''
    for fastq in fastqFiles:
        #------- define output sam file name ----
        if fastq[0].endswith(".fastq.gz"):
            output = fastq[0][:-8] + 'sam'
        else:
            output = fastq[0][:-5] + 'sam'
        map_result.append(output)
        #-------- map without annotation ---------
        if annotation == '':
            if len(fastq) == 2:
                cmd = cmd + ('gsnap --input-buffer-size=5000 -D {db_path} -d {db_name} --gunzip -A sam '
                             '-B 4 -i 2 -t {thread} '
                            '-N 1 {fastq1} {fastq2} > {output} && ').format(
                            db_path=db_path,db_name=db_name,thread=thread,
                            fastq1=fastq[0],fastq2=fastq[1],output=output)
            else:
                cmd = cmd + ('gsnap --input-buffer-size=5000 -D {db_path} -d {db_name} -A sam '
                             '--gunzip -B 4 -i 2 -t {thread} '
                            '-N 1 {fastq} > {output} && ').format(db_path=db_path,
                            db_name=db_name,thread=thread,
                            fastq=fastq[0],output=output)
        #-------- map with annotation ------------
        else:
            if len(fastq) == 2:
                cmd = cmd + ('gsnap --input-buffer-size=5000 -D {db_path} -d {db_name} --gunzip -A sam '
                             '-B 4 -i 2 -t {thread} '
                            '-N 1 -s {annotation} {fastq1} {fastq2} --force-xs-dir > {output} && ').format(
                            db_path=db_path,db_name=db_name,thread=thread,annotation=annotation,
                            fastq1=fastq[0],fastq2=fastq[1],output=output)
            else:
                cmd = cmd + ('gsnap --input-buffer-size=5000 -D {db_path} -d {db_name} -A sam --gunzip '
                             '-B 4 -i 2 -t {thread} '
                            '-N 1 -s {annotation} {fastq} --force-xs-dir > {output} && ').format(db_path=db_path
                            ,db_name=db_name,thread=thread,annotation=annotation,
                            fastq=fastq[0],output=output)
    subprocess.check_call(cmd[:-3],shell=True)
    return map_result

#=========  bowtie2 alignment  ========================= 
def bowtie2(fastqFiles,database,thread=1):
    map_result = []
    cmd = ''
    for fastq in fastqFiles:
        if fastq[0].endswith(".fastq.gz"):
            output = fastq[0][:-8] + 'sam'
        else:
            output = fastq[0][:-5] + 'sam'
        map_result.append(output)
        if len(fastq) == 2:
            cmd = cmd + ('bowtie2 -x {database} -p {thread} -1 {fastq1} -2 {fastq2} '
            '-S {outputfile} & ').format(database=database,thread = thread, 
                             fastq1 = fastq[0], fastq2 = fastq[1],outputfile=output)
        else:
            cmd = cmd + ('bowtie2 -x {database} -p {thread} -U {fastq} '
                         '-S {outputfile} && ').format(database=database,
                        thread = thread, fastq = fastq[0], outputfile=output)
    subprocess.check_call(cmd[:-3],shell=True)
    return map_result
#============  Tophat Alignment  =======================            
def tophat(fastqFiles,database,annotation,thread=1):
    """
    This function will map the fastq files to the reference
    And will returen the path of the bam file results,
    In each folder, there is a file 'accepted_hits.bam' which is important
    """
    map_result = []
    cmd = ''
    for fastq in fastqFiles:
        if len(fastq) == 2:
            output = fastq[0][:-8] + '_map'  # define the output folder
            map_result.append(output)
            cmd = cmd + ('tophat -p {thread} -G {gff} --b2-very-sensitive -o {output} ' 
            '{reference} {fastq1} {fastq2} & ').format(thread = thread, gff = annotation, 
            output = output,reference = database, fastq1 = fastq[0],fastq2 = fastq[1])
        else:
            output = fastq[0][:-8] + '_map'
            map_result.append(output)
            cmd = cmd + ('tophat -p {thread} -G {gff} --b2-very-sensitive -o {output} '
            '{reference} {fastq} && ').format(thread = thread, gff = annotation, 
            output = output,reference = database, fastq = fastq[0])
    subprocess.check_call(cmd[:-3],shell=True)
    return map_result

#============  Blast Alignment  ===========================
def blastn(faFiles,database,thread = 1):
    """
    this function run blastn,default model is running through internet 
    """
    map_result = []
    for fa in faFiles:
        output = fa[:-2] + 'blast.txt'
        map_result.append(output)
        blastn = ('blastn -query {input} -task megablast -out {output} '
                  '-db {db} -evalue 1e-10 -word_size 10 -outfmt 6 '
                  '-num_alignments 1 -num_threads {thread} ').format(
                  input=fa,output=output,db=database,thread=thread)
        subprocess.check_call(blastn,shell=True)
    return map_result
#============ bwa alignment  ===============================
def bwa_vari(readgroup,fqFiles,database,thread=1):
    """
    this function run bwa, whose downstream analysis is to find variant
    """
    map_result = []
    bwaCmd = ''
    for fastq,rg in zip(fqFiles,readgroup):
        if fastq[0].endswith(".fastq.gz"):
            output = fastq[0][:-8] + 'sam'
        else:
            output = fastq[0][:-5] + 'sam'
        map_result.append(output)
        if len(fastq) == 2:
            bwaCmd = bwaCmd + ('bwa mem -t {thread} -M -R {readgroup} '
                      '{database} {fq1} {fq2} > {output} && ').format(
                        thread=thread,readgroup=rg,database=database,
                        fq1=fastq[0],fq2=fastq[1],output=output)
        else:
            bwaCmd = bwaCmd + ('bwa mem -t {thread} -M -R {readgroup} '
                       '{database} {fq} > {output} && ').format(thread=thread,
                        readgroup=rg,database=database,fq=fastq[0])
    subprocess.check_call(bwaCmd[:-3],shell=True)
    return map_result        
#============  STAR alignment  ===============================
def STAR_Db(db_path,ref_fa,thread=1,annotation = ''):
    """
    This function generates database for alignment using STAR
    """
    if listdir(db_path) == []:
        cmd = ('STAR --runMode genomeGenerate --genomeDir {db_path} '
               '--genomeFastaFiles {ref_fa} --runThreadN {thread} '
               '--limitGenomeGenerateRAM 310000000000 ').format(
                db_path=db_path,ref_fa=ref_fa,thread=thread)
        if annotation != '':
            cmd = cmd + ('--sjdbGTFfile {gff3} --sjdbGTFtagExonParentTranscript Parent '
                         '--sjdbOverhang 100').format(gff3=annotation)
    subprocess.check_call(cmd,shell=True)
    
def STAR(fastqFiles,db_path,thread=1,otherParameters=['']):
    """
    STAR are more proper for aligning RNA seq
    otherParameters: a list of added parameters
    """
    map_results = []
    cmd = ''
    for fastq in fastqFiles:
        #------- define output sam file name ----
        if fastq[0].endswith(".fastq.gz"):
            output = fastq[0][:-8] 
        else:
            output = fastq[0][:-5]
        map_results.append(output + 'Aligned.out.sam')
        #-------- map without annotation ---------
        if len(fastq) == 2:
            starCmd = ('STAR --genomeDir {ref} '
                         '--readFilesCommand zcat '
                         '--readFilesIn {fq1} {fq2} --runThreadN '
                         '{thread} --outFileNamePrefix {output} '
                         '--outSAMunmapped Within').format(
                        ref=db_path,fq1=fastq[0],fq2=fastq[1],
                        thread=thread,output=output)
            cmd = cmd + starCmd + ' ' + ' '.join(otherParameters) + ' && '
        else:
            starCmd = ('STAR --genomeDir {ref} '
                         '--readFilesCommand zcat '
                         '--readFilesIn {fq1} --runThreadN '
                         '{thread} --outFileNamePrefix {output} '
                         '--outSAMunmapped Within').format(
                        ref=db_path,fq1=fastq[0],
                        thread=thread,output=output)
            cmd = cmd + starCmd + ' ' + ' '.join(otherParameters) + ' && '
    subprocess.check_call(cmd[:-3],shell=True)
    final_name = []
    for sam in map_results:
        new_name = sam[:-15] + 'sam'
        rename = ('mv {star_result} {modified_name}').format(star_result=sam,
                  modified_name=new_name)
        final_name.append(new_name)
        subprocess.check_call(rename,shell=True)
    return final_name

#============  STAR 2 pass alignment  ===============================
def STAR2Pass(fastqFiles,starDb,ref_fa,thread=1):
    # this function run 2 pass of mapping
    map_sams = STAR(fastqFiles,starDb,thread)
    if not os.path.exists("starDb2Pass"):
        subprocess.check_call('mkdir starDb2Pass',shell=True)
    for sam,fastq in zip(map_sams,fastqFiles):
        SJFile = sam[:-3] + 'SJ.out.tab'
        cmd = ('STAR --runMode genomeGenerate --genomeDir starDb2Pass '
               '--genomeFastaFiles {ref_fa} --sjdbFileChrStartEnd {SJ} '
               '--sjdbOverhang 100 --runThreadN {thread} '
               '--limitGenomeGenerateRAM 310000000000').format(
                ref_fa=ref_fa,SJ=SJFile,thread=thread)
        subprocess.check_call(cmd,shell=True)
        result = STAR([fastq],'starDb2Pass',thread)
    subprocess.check_call('rm -r starDb2Pass',shell=True)
    return map_sams

