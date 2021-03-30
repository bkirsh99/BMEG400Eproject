library(dplyr)
library(tidyr)
library(readxl)
library(reshape2)
library(purrr)
library(ggplot2)
library(ggvenn)
library(ggVennDiagram)

dir <- "C:/Users/bkirs/Documents/School/BMEG400E/BMEG400Eproject"
setwd(dir)

data <- read_excel("225SNPs26populations.xlsx")
colnames(data) <- data[2, ] # the second row will be the header
data <- data[-c(1,2), ]          # removing the first two rows
#The first three columns are SNP information. After that, each of the 26 populations has 5 columns:
#("population","effect allele number","other allele number","total allele number","effect allele frequency")
#The remaining columns (138-[3+26*5=133]=5) summarize global population data ("population","effect allele number","other allele number","effect allele frequency","GWAS P-value")

populations <- c("CEU","ASW","YRI")
cols <- c("population","effect allele number","other allele number","total allele number","effect allele frequency","GWAS P-value")

SNP.df <- data[,1:3]

MAF.df <- data.frame(row.names= SNP.df$`SNP ID`)

ALL.df <- data.frame(matrix(nrow=nrow(data),ncol=length(cols)))
colnames(ALL.df) <- cols
ALL.df$population <- "ALL"
ALL.df$`GWAS P-value` <- data$`GWAS P-value`
ALL.df[is.na(ALL.df)] <- 0

dfs <- paste(populations,"df",sep=".")
hyper.df <- data.frame()

for(i in seq(from=4, to=ncol(data), by=5)){
  population <- toString(data[1,i])
  if(population %in% populations){
    MAF.df[population] <- data[,i+4]
    assign(paste(population,"df",sep="."), data[,i:(i+4)])
    ALL.df["effect allele number"] <- as.numeric(ALL.df[["effect allele number"]]) + as.numeric(data[["effect allele number"]])
    ALL.df["other allele number"] <- as.numeric(ALL.df[["other allele number"]]) + as.numeric(data[["other allele number"]])
    ALL.df["total allele number"] <- as.numeric(ALL.df[["total allele number"]]) + as.numeric(data[["total allele number"]])
    ALL.df["effect allele frequency"] <- as.numeric(ALL.df[["effect allele frequency"]]) + as.numeric(data[["effect allele frequency"]])
  }
}

for(df in dfs){
  a <- c()
  info <- get(df)
  for(j in 1:nrow(info)){
    m <- as.numeric(ALL.df$`total allele number`[j]) #failure in population size
    n <- as.numeric(ALL.df$`effect allele number`[j]) #success in population
    k <- as.numeric(info$`total allele number`[j]) #sample size
    x <- as.numeric(info$`effect allele number`[j]) #success in sample
    a[j] <- dhyper(x = x, m = m, n = n, k = k)
  }
  hyper.df[df] <- a
}

for(i in 1:nrow(ALL.df)){
  a <- c()
  for (df in dfs){
    info <- get(df)
    m <- as.numeric(info$`total allele number`[i]) #effect alleles in population
    n <- as.numeric(ALL.df$`total allele number`[i])-m #effect alleles NOT in population
    k <- as.numeric(ALL.df$`effect allele number`[i]) #effect alleles hit
    x <- as.numeric(info$`effect allele number`[i])
    a <- cbind(a,dhyper(x=x,m=m,n=n,k=k))
  }
  hyper.df <- rbind(hyper.df,a)
}

colnames(hyper.df) <- populations
rownames(hyper.df) <- SNP.df$`SNP ID`

hyper_adj.df <- p.adjust(as.matrix(hyper.df))

cutoff <- 0.01 / (2*3*225) #To control a family-wise error rate (FWER) of 0.01, we used a raw p-value of 7.41E-6 as cutoff. 
pvalFiltered.df <- subset(hyper.df < cutoff)

enriched <- data.frame(hyper.df < cutoff)
table(rowSums(enriched))
v <- ggvenn(enriched)


MAF.df["id"] <- rownames(MAF.df)
MAF_molten.df <- melt(MAF.df, id.vars="id", value.name="MAF", variable.name="Population")
p <- ggplot(MAF_molten.df,aes(x=id,y=MAF)) + geom_bar(position="identity", stat="identity") + facet_wrap(~Population,nrow=3)
p
q <- ggplot(MAF_molten.df,aes(x=id,y=MAF,fill=Population)) + geom_dotplot(binaxis='y', stackdir='center')
q
