import subprocess
import os
from os import listdir
import copy
#==========  gsnap alignment  ==========================
def gsnap_Db(fa,db_path,db_name,annotation):
    """
    This function builds index for genome
    """
    cmd = ('gmap_build -D {db_path} -d {db_name} {fa}').format(db_path=db_path,db_name=db_name,fa=fa)
    subprocess.call(cmd,shell=True)
    if annotation != '':
        cmd = ('iit_store -G -o {db_name} {annotation}').format(db_name=db_name,annotation=annotation)
    subprocess.call(cmd,shell=True)
    # move the the target maps folder
    cmd = ('mv {iit_file} {folder}').format(iit_file=db_name+'.iit',folder=db_path+'/'+db_name+'/'+db_name+'.maps')
    subprocess.call(cmd,shell=True)
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
    print cmd[:-3]
    subprocess.call(cmd[:-3],shell=True)
    return map_result

#=========  bowtie2 alignment  ========================= 
def bowtie2(fastqFiles,database,thread=1,otherParameters=['']):
    """
    This function runs bowtie2

    * fastqFiles: list of fastqFiles. paired_end:[[f1.fq.gz,f2.fq.gz]]. single_end:[[f1.fq.gz]]
    * database: index of bowtie2
    * thread: number of thread to use
    """
    map_result = []
    nomap_fqs = []
    return_non = ''
    cmd = '' # final command
    parameters = copy.copy(otherParameters)
    for fastq in fastqFiles:
        if fastq[0].endswith(".fastq.gz"):
            output = fastq[0][:-8] + 'sam'
        else:
            output = fastq[0][:-5] + 'sam'
        map_result.append(output)
        # define otherParameters for extracting unmappd reads
        if '--un-gz' in otherParameters:
            nomap = output[:-3] + 'norna.fq.gz'
            return_non = 'yes'
        # define align command
        if len(fastq) == 2:
            cmd = cmd + ('bowtie2 -x {database} -p {thread} -1 {fastq1} -2 {fastq2} '
            '-S {outputfile} ').format(database=database,thread = thread, 
                             fastq1 = fastq[0], fastq2 = fastq[1],outputfile=output)
            # output non map fq.gz
            if return_non=='yes':
                nomap_fqs.append([nomap[:-6]+'_1.fq.gz',nomap[:-6]+'_2.fq.gz'])
                
        else:
            cmd = cmd + ('bowtie2 -x {database} -p {thread} -U {fastq} '
                         '-S {outputfile} ').format(database=database,
                        thread = thread, fastq = fastq[0], outputfile=output)
            # otuput non map fq.gz
            if return_non=='yes':
                nomap_fqs.append([nomap])
        otherParameters[otherParameters.index('--un-gz')] = '--un-gz ' + nomap
        cmd = cmd + ' '.join(otherParameters) + ' && '
        otherParameters = copy.copy(parameters)
    print cmd[:-3]
    subprocess.check_call(cmd[:-3],shell=True)
    if return_non=='yes':
        return nomap_fqs
    else:
        return map_result

def bowtie_rRNA(fastqFiles,database,thread=1):
    """
    This function runs bowtie to align reads to rRNA and output the reads that don't align to rRNA
    
    * fastqFiles: list of fastqFiles. paired_end:[[f1.fq.gz,f2.fq.gz]]. single_end:[[f1.fq.gz]]
    * database: index of bowtie2
    * thread: number of thread to use
    """
    map_result = []
    for fastq in fastqFiles:
        if fastq[0].endswith(".fastq.gz"):
            output = fastq[0][:-8] + 'norna.fq'
        else:
            output = fastq[0][:-5] + 'norna.fq'
        map_result.append([output+'.gz'])
        if len(fastq) == 1:
            cmd = ('bowtie -l 23 -p {thread} --un {out} -S {db} <(gunzip -c {input}) > temp.sam').format(
                           thread=thread,out=output,db=database,input=fastq[0])
            subprocess.call(cmd,shell=True,executable='/bin/bash')
            os.remove('temp.sam')
            subprocess.call('gzip {f}'.format(f=output),shell=True)
    
    return map_result
#============  Tophat Alignment  =======================            
def tophat(fastqFiles,database,annotation,thread=1,otherParameters=['']):
    """
    This function maps the fastq files to the reference
    
    * fastqFiles: list. List of fastq files. eg: [[f1.fq.gz],[f2,fq.gz],...] or [[f1.fq.gz,f2.fq.gz],...]
    * database: str. Name of index generated by bowtie or bowtie2.
    * annotation: str. Annotation file name. eg: name.gff
    
    And will returen the path of the bam file results,
    In each folder, there is a file 'accepted_hits.bam' which stores the results
    """
    map_result = []
    cmd = ''
    for fastq in fastqFiles:
        if len(fastq) == 2:
            output = fastq[0][:-5] + 'tophat'  # define the output folder
            map_result.append(output)
            tophat_cmd1 = ('tophat -p {thread} -G {gff} -o {output} ' 
            '--transcriptome-index tophatTranscriptIndex/trptIndex ').format(thread=thread,
                                                            gff=annotation,output=output)
            tophat_cmd2 = (' {reference} {fq1} {fq2}').format(reference=database,fq1=fastq[0],
                                                                 fq2=fastq[1])
            cmd = cmd + tophat_cmd1 + ' '.join(otherParameters) + tophat_cmd2 + ' && '
        else:
            output = fastq[0][:-5] + 'tophat'
            map_result.append(output)
            tophat_cmd1 = ('tophat -p {thread} -G {gff} -o {output} '
                           '--transcriptome-index tophatTranscriptIndex/trptIndex ').format(thread = thread,
                            gff = annotation,output = output)
            tophat_cmd2 = (' {reference} {fq}').format(reference=database,fq=fastq[0])
            cmd = cmd + tophat_cmd1 + ' '.join(otherParameters) + tophat_cmd2 + ' && '
    print cmd[:-3]
    subprocess.call(cmd[:-3],shell=True)
    return map_result

    
#============  Blast Alignment  ===========================
def makeblast(ref_fa,db_type,out):
    """
    This function make index for the reference fa file
    
    * ref_fa: reference fa file
    * db_type: database type. Default: 'nucl', Alternate: 'prot'
    * out: name for created database
    """
    cmd = ('makeblastdb -in {ref} -dbtype {type} -out {out} -title {title}').format(
            ref=ref_fa,type=db_type,out=out,title=out)
    subprocess.call(cmd.split(' '))
    

def blastn(faFiles,database,thread = 1):
    """
    this function run blastn,default model is running through internet 
    """
    map_result = []
    for fa in faFiles:
        if fa.endswith('.gz'):
            output = fa[:-5] + 'blast.txt'
            fa = '<(gunzip -c {fa})'.format(fa=fa)
        else:
            if fa.endswith('.fasta'):
                output = fa[:-5] + 'blast.txt'
            else:
                output = fa[:-2] + 'blast.txt'
        map_result.append(output)
        blastn = ('blastn -query {input} -task megablast -out {output} '
                  '-db {db} -evalue 1e-7 -word_size 10 -outfmt 6 '
                  '-num_alignments 1 -num_threads {thread} ').format(
                  input=fa,output=output,db=database,thread=thread)
        print blastn
        subprocess.check_call(blastn,shell=True,executable='/bin/bash')
    return map_result
#============ bwa alignment  ===============================
def bwa_vari(readgroup,fqFiles,database,thread=1):
    """
    this function run bwa, whose downstream analysis is to find variant
    """
    map_result = []
    bwaCmd = ''
    for fastq,rg in zip(fqFiles,readgroup):
        rg = rg + '\\tPL:illumina\\tLB:lib20000\\tPU:unit1'
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
    print bwaCmd[:-3]
    subprocess.call(bwaCmd[:-3],shell=True)
    return map_result        
#============  STAR alignment  ===============================
def STAR_Db(db_path,ref_fa,thread=1,annotation = ''):
    """
    This function generates database for alignment using STAR
    """
    if not os.path.exists(db_path): os.mkdir(db_path)
    if listdir(db_path) == []:
        cmd = ('STAR --runMode genomeGenerate --genomeDir {db_path} '
               '--genomeFastaFiles {ref_fa} --runThreadN {thread} '
               '--limitGenomeGenerateRAM 100000000000 ').format(
                db_path=db_path,ref_fa=ref_fa,thread=str(thread))
        if annotation != '':
            cmd = cmd + ('--sjdbGTFfile {gff3} --sjdbGTFtagExonParentTranscript Parent '
                         '--sjdbOverhang 100').format(gff3=annotation)   # for geneDb add --sjdbGTFfeatureExon CDS
    print cmd
    subprocess.check_call(cmd,shell=True)
    
def STAR(fastqFiles,db_path,thread=1,annotation='',otherParameters=['']):
    """
    STAR are more proper for aligning RNA seq
    otherParameters: a list of added parameters
    """
    map_results = []
    cmd = ''
    if annotation != '':
            otherParameters.extend(['--sjdbGTFfile {gff}'.format(gff=annotation), '--sjdbGTFtagExonParentTranscript Parent'])
    for fastq in fastqFiles:
        #------- define output sam file name ----
        if fastq[0].endswith(".fastq.gz"):
            output = fastq[0][:-8] 
        else:
            output = fastq[0][:-5]
        if '--outSAMtype BAM SortedByCoordinate' in otherParameters:
            map_results.append(output + 'Aligned.sortedByCoord.out.bam')
        else:
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
    print cmd[:-3]
    subprocess.check_call(cmd[:-3],shell=True)
    final_name = []
    for sam in map_results:
        new_name = sam.split('.')[0] + '.' + sam[-3:]
        rename = ('mv {star_result} {modified_name}').format(star_result=sam,
                  modified_name=new_name)
        final_name.append(new_name)
        print rename
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
               '--limitGenomeGenerateRAM 110000000000').format(
                ref_fa=ref_fa,SJ=SJFile,thread=thread)
        subprocess.check_call(cmd,shell=True)
        result = STAR([fastq],'starDb2Pass',thread)
    subprocess.check_call('rm -r starDb2Pass',shell=True)
    return map_sams


    
