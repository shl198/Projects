library(ShortRead)
file1 <- '/data/shangzhong/Diff_express/test/f1.fastq.gz'
f1 <- FastqSampler(file1, n=1e7)
set.seed(123L); p1 <- yield(f1)
writeFastq(p1,'/data/shangzhong/Diff_express/t1.fq.gz')

file2 <- '/data/shangzhong/Diff_express/test/f2.fastq.gz'
f2 <- FastqSampler(file2, n=1e7)
set.seed(123L); p2 <- yield(f2)
writeFastq(p1,'/data/shangzhong/Diff_express/t2.fq.gz')