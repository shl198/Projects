dir <- getwd()
res <- list.files(dir,pattern='\\.csv$')
res <- res[order(res)]
res <- res[-1]
res

result <- list()
for (i in 1:length(res)){
    data <- read.csv(res[i])
    data <- data$X
    result[[length(result) + 1]] <- data
}


# statistics
num <- c()
for (i in 1:length(result)){
    num <- c(num,length(result[[i]]))
}
num

# GS VS Cosmc3C9 and Cosmc 4B7
length(result[[1]])
length(result[[5]])
length(intersect(result[[1]],result[[5]]))


# GS VS Cosmc 3C9
length(result[[2]])
length(result[[3]])
length(result[[4]])
length(intersect(intersect(result[[2]],result[[3]]),result[[4]]))

# GS VS Cosmc 4B7
length(result[[6]])
length(result[[7]])
length(result[[8]])
length(intersect(intersect(result[[6]],result[[7]]),result[[8]]))

# Cosmc 3C9 VS Cosmc 4B7
length(result[[9]])

# GS VS 91-1C8 and 91-2A6
length(result[[10]])
length(result[[14]])
length(intersect(result[[10]],result[[14]]))

# GS VS 91-1C8
length(result[[11]])
length(result[[12]])
length(result[[13]])
length(intersect(intersect(result[[11]],result[[12]]),result[[13]]))

# GS VS 91-2A6
length(result[[15]])
length(result[[16]])
length(result[[17]])
length(intersect(intersect(result[[15]],result[[16]]),result[[17]]))

# 91-1C8 VS 91-2A6
length(result[[18]])

install.packages('VennDiagram')
library(VennDiagram)

par(mfrow=c(2,1))
par(mar=c(5,10,5,10))
# plot GS_3C9_4B7
venn.plot <- draw.pairwise.venn(2165,2751,1646,
                   category=c('GS_VS_3c9','GS_VS_4B7'),
                   euler.d=T,scaled=T,inverted=F,ext.text=T,
                   ext.percent=rep(0.05,3),lwd=c(3,3),
                   col=rep('black',2),fill=c('red','green'),
                   alpha=c(0.7,0.7),cat.pos=c(-30,30),
                   cat.dist=rep(0.03,2),cat.cex=rep(1.1,2),
                   cat.col=c('red','green'),ext.pos=30)
grid.draw(venn.plot)

# plot GS_91-1C8_91-2A6
venn.plot <- draw.pairwise.venn(2165,2751,1646,
                                category=c('GS_VS_1C8','GS_VS_2A6'),
                                euler.d=T,scaled=T,inverted=F,ext.text=T,
                                ext.percent=rep(0.05,3),lwd=c(3,3),
                                col=rep('black',2),fill=c('red','green'),
                                alpha=c(0.7,0.7),cat.pos=c(-150,150),
                                cat.dist=rep(0.03,2),cat.cex=rep(1.1,2),
                                cat.col=c('red','green'),ext.pos=30)

# plot GS_VS_3c9.     7,8,9
paste('7 and 8:',length(intersect(result[[2]],result[[3]])))
paste('8 and 9:',length(intersect(result[[3]],result[[4]])))
paste('7 and 9:',length(intersect(result[[2]],result[[4]])))

venn.plot <- draw.triple.venn(990,1175,1029,826,847,843,754,
                  category=c('GS_VS_3C9-7','GS_VS_3C9-8','GS_VS_3C9-9'),
                  rotation=1,reverse=F,euler.d=T,scaled=T,lwd=rep(3,3),
                  col=rep('black',3),fill=c('red','green','blue'),
                  alpha=rep(0.7,3),
                  cat.pos=c(-30,30,180),cat.col=c('red','green','blue'),
                  cat.cex=rep(1.1,3),
                  rotation.degree=0)

# plot GS_VS_4B7 10,11,12
paste('10 and 11:',length(intersect(result[[6]],result[[7]])))
paste('11 and 12:',length(intersect(result[[7]],result[[8]])))
paste('10 and 12:',length(intersect(result[[6]],result[[8]])))

venn.plot <- draw.triple.venn(295,318,340,223,228,225,199,
                              category=c('GS_VS_4B7-10','GS_VS_4B7-11','GS_VS_4B7-12'),
                              rotation=1,reverse=F,euler.d=T,scaled=T,lwd=rep(3,3),
                              col=rep('black',3),fill=c('red','green','blue'),
                              alpha=rep(0.7,3),
                              cat.pos=c(-30,30,180),cat.col=c('red','green','blue'),
                              cat.cex=rep(1.1,3),
                              rotation.degree=0)


paste('13 and 14:',length(intersect(result[[11]],result[[12]])))
paste('14 and 15:',length(intersect(result[[12]],result[[13]])))
paste('13 and 15:',length(intersect(result[[11]],result[[13]])))

venn.plot <- draw.triple.venn(1593,1603,1900,1317,1189,1228,1116,
                              category=c('GS_VS_1C8-13','GS_VS_1C8-14','GS_VS_1C8-15'),
                              rotation=1,reverse=F,euler.d=T,scaled=T,lwd=rep(3,3),
                              col=rep('black',3),fill=c('red','green','blue'),
                              alpha=rep(0.7,3),
                              cat.pos=c(-30,30,180),cat.col=c('red','green','blue'),
                              cat.cex=rep(1.1,3),
                              rotation.degree=0)

paste('16 and 17:',length(intersect(result[[15]],result[[16]])))
paste('17 and 18:',length(intersect(result[[16]],result[[17]])))
paste('16 and 18:',length(intersect(result[[15]],result[[17]])))

venn.plot <- draw.triple.venn(2052,1973,1976,1707,1627,1568,1480,
                              category=c('GS_VS_2A6-16','GS_VS_2A6-17','GS_VS_2A6-18'),
                              rotation=1,reverse=F,euler.d=T,scaled=T,lwd=rep(3,3),
                              col=rep('black',3),fill=c('red','green','blue'),
                              alpha=rep(0.7,3),
                              cat.pos=c(-30,30,180),cat.col=c('red','green','blue'),
                              cat.cex=rep(1.1,3),
                              rotation.degree=0)

