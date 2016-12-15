require(ggplot2)
require(pcadapt)
source("~/Documents/thesis/git/FDR/fdrUtils.R")
source("~/Documents/thesis/git/pooled-samples/poolUtils.R")
source("~/Documents/thesis/git/pooled-samples/poolOld.R")

case <- "IM"
min.maf <- 0.05

if (case=="DM"){
  geno.matrix <- read.table(system.file("extdata","geno3pops",package="pcadapt"))
  gt <- 1:150
} else if (case=="IM"){
  geno.matrix <- read.table("~/Documents/thesis/MolEcolRes/IM/dems_3_neu_10_out_001_n150_p1000")
  gt <- read.table("~/Documents/thesis/MolEcolRes/IM/dems_3_neu_10_out_001_n150_p1000.gt")$V1
}

pop <- c(rep(1,50),rep(2,50),rep(3,50))
pool.matrix <- get.pool.matrix(data=geno.matrix,pop=pop,ploidy=2)
cover.matrix <- simulate.cover(nrow=nrow(pool.matrix),ncol=ncol(pool.matrix),
                               min.cover=c(100,1000),max.cover=c(200,2000),
                               high.cov.loc = 350:400)
pool.matrix.test <- cover.to.pool(data=geno.matrix,cover.matrix = cover.matrix,pop = pop)
new.geno <- sample.geno(pool.matrix = pool.matrix.test,cover.matrix,nINDperPOOL = c(100,100,100))
bino.geno <- sample.geno(pool.matrix.test,nINDperPOOL = c(100,100,100))
ind.geno <- sample.geno(pool.matrix.test,cover.matrix=cover.matrix,nINDperPOOL = c(100,100,100),method = "per.ind")

## Get the statistics
filename <- read.pcadapt(new.geno, type = "lfmm", local.env = TRUE)
x <- pcadapt(filename,K=2,min.maf = min.maf)
filename.bino <- read.pcadapt(bino.geno, type = "lfmm", local.env = TRUE)
x.bino <- pcadapt(filename.bino,K=2,min.maf = min.maf)
filename.ind <- read.pcadapt(ind.geno,type="lfmm",local.env=TRUE)
x.ind <- pcadapt(filename.bino,K=2,min.maf = min.maf)
y <- pool.old(pool.matrix.test,K=2)
z <- pool.old.corrected(pool.matrix.test,K=2,min.maf = min.maf,cover.matrix = cover.matrix)

stat.multiple <- 0
for (n.test in 1:10){
  multiple.geno <- sample.geno(pool.matrix = pool.matrix.test,cover.matrix,nINDperPOOL = c(20,20,20))
  filename.multiple <- read.pcadapt(multiple.geno, type = "lfmm", local.env = TRUE)
  x.multiple <- pcadapt(filename.multiple,K=2,min.maf = min.maf)
  stat.multiple <- stat.multiple + x.multiple$chi2.stat
}

## Order the statistics
rnk.1 <- sort(x$stat,decreasing = TRUE,index.return=TRUE)$ix
rnk.2 <- sort(y$stat,decreasing = TRUE,index.return=TRUE)$ix
rnk.3 <- sort(x.bino$stat,decreasing = TRUE,index.return=TRUE)$ix
rnk.4 <- sort(z$stat,decreasing = TRUE,index.return=TRUE)$ix
rnk.5 <- sort(x.ind$stat,decreasing = TRUE,index.return=TRUE)$ix
rnk.6 <- sort(stat.multiple,decreasing = TRUE,index.return=TRUE)$ix

## Create the data frames
df.1 <- create.fdr.pow(rnk.1,ground.truth = gt,soft.name = "Beta-binomial sampling",smooth = TRUE)
df.2 <- create.fdr.pow(rnk.2,ground.truth = gt,soft.name = "First version",smooth = TRUE)
df.3 <- create.fdr.pow(rnk.3,ground.truth = gt,soft.name = "Binomial sampling",smooth = TRUE)
df.4 <- create.fdr.pow(rnk.4,ground.truth = gt,soft.name = "Coverage correction",smooth = TRUE)
df.5 <- create.fdr.pow(rnk.5,ground.truth = gt,soft.name = "Beta-binomial per individual",smooth = TRUE)
df.6 <- create.fdr.pow(rnk.6,ground.truth = gt,soft.name = "Multiple runs",smooth = TRUE)

## Bind the data frames
df <- rbind(df.1,df.2,df.3,df.4,df.5,df.6)
df[,1] <- as.character(df[,1])
df[,2] <- as.numeric(as.character(df[,2]))
df[,3] <- as.numeric(as.character(df[,3]))
colnames(df) <- c("Software","FDR","Power")

## ggplot
p0 <- ggplot(data = df,aes(x=FDR,y=Power)) + 
  geom_line(aes(linetype=Software, color=Software),size=2,na.rm = TRUE) +
  xlim(0,1) + ylim(0,1) +
  theme_bw() +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"),
        title=element_text(size=15,face="bold"),
        legend.text=element_text(size=15),
        legend.key.height=unit(2,"line"),
        legend.key.width=unit(2,"line")
  )
print(p0)

