require(ggplot2)
source("~/Documents/thesis/git/FDR/fdrUtils.R")
source("~/Documents/thesis/git/pooled-samples/poolUtils.R")
source("~/Documents/thesis/git/pooled-samples/poolOld.R")

geno.matrix <- read.table("~/Documents/thesis/Datasets/geno3pops")
pop <- c(rep(1,50),rep(2,50),rep(3,50))
pool.matrix <- get.pool.matrix(data=geno.matrix,pop=pop,ploidy=2)
cover.matrix <- simulate.cover(nrow(pool.matrix),ncol(pool.matrix),min.cover=c(40,5000),max.cover=c(50,6000))
new.geno <- sample.geno(pool.matrix,cover.matrix,nINDperPOOL = c(200,200,200))
bino.geno <- sample.geno(pool.matrix,nINDperPOOL = c(200,200,200))

## Get the statistics
filename <- read.pcadapt(new.geno, type = "lfmm", local.env = TRUE)
x <- pcadapt(filename,K=2)
filename.bino <- read.pcadapt(bino.geno, type = "lfmm", local.env = TRUE)
x.bino <- pcadapt(filename.bino,K=2)
y <- pool.old(pool.matrix,K=2)
z <- create.pcadapt.pool(pool.matrix,K=2,min.maf=0.05,cover.matrix = cover.matrix)

## Order the statistics
rnk.1 <- sort(x$stat,decreasing = TRUE,index.return=TRUE)$ix
rnk.2 <- sort(y$stat,decreasing = TRUE,index.return=TRUE)$ix
rnk.3 <- sort(x.bino$stat,decreasing = TRUE,index.return=TRUE)$ix
rnk.4 <- sort(z$stat,decreasing = TRUE,index.return=TRUE)$ix

## Create the data frames
df.1 <- create.fdr.pow(rnk.1,ground.truth = 1:150,soft.name = "Beta-binomial sampling",lmax = 1000,smooth = TRUE)
df.2 <- create.fdr.pow(rnk.2,ground.truth = 1:150,soft.name = "First version",lmax = 1000,smooth = TRUE)
df.3 <- create.fdr.pow(rnk.3,ground.truth = 1:150,soft.name = "Binomial sampling",lmax = 1000,smooth = TRUE)
df.4 <- create.fdr.pow(rnk.4,ground.truth = 1:150,soft.name = "Coverage correction",lmax = 1000,smooth = TRUE)

## Bind the data frames
df <- rbind(df.1,df.2,df.3,df.4)
df[,1] <- as.character(df[,1])
df[,2] <- as.numeric(as.character(df[,2]))
df[,3] <- as.numeric(as.character(df[,3]))
colnames(df) <- c("Software","FDR","Power")

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

