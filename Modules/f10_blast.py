import subprocess,os

def makeblastdb(fastaFile,datatype,outputname):
    """
    this function build database given a fasta file
    
    * fastaFile: can be gzipped or not
    """
    if fastaFile.endswith('.gz'):
        cmd = ('gunzip -c {input} | makeblastdb -in - -dbtype {type} -title {title} '
               '-out {outputname}').format(input=fastaFile,
                    type=datatype,outputname=outputname,title=outputname)
    else:
        cmd = ('makeblastdb -in {input} -dbtype {type} -title {title} '
               '-out {outputname}').format(input=fastaFile,
                    type=datatype,outputname=outputname,title=outputname)
    subprocess.call(cmd,shell=True)
    

def blastp(query,database,outputFile,threads,evalue,fmt,mapnum):
    """
    This function run blastp
    
    * query: fasta file which you want to map
    
    * database: database path/name
    
    * outputFile: tabular blast result
    """
    if query.endswith('.gz'):
        cmd = ('gunzip -c {input} | blastp -query - -db {database} '
               '-out {outputFile} -evalue {evalue} -outfmt {format} '
               '-seg yes -num_threads {thread} -num_alignments {mapnum}').format(input=query,
                    database=database,outputFile=outputFile,evalue=evalue,
                    format=str(fmt),thread=str(threads),mapnum=mapnum)
    else:
        cmd = ('blastp -query {input} -db {database} -out {outputFile} '
               '-evalue {evalue} -outfmt {format} -seg yes '
               '-num_threads {thread} -num_alignments {mapnum}').format(input=query,
                    database=database,outputFile=outputFile,evalue=evalue,
                    format=str(fmt),thread=str(threads),mapnum=mapnum)
    subprocess.call(cmd,shell=True)
