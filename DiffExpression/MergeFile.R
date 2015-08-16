filepath = '/data/shangzhong/finalResult'
files = list.files(filepath,pattern='\\.txt$')
data <- read.table(files[1])
for (i in 2:length(files)) {
    data1 <- read.table(files[i])
    data <- merge(data,data1,by='V1')
}
logic <- rowSums(data[,-1]) > 0
res <- data[logic,]
name <- as.vector(res[,1])
write(name,file='name.txt')
