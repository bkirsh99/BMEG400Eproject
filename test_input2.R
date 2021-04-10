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

#Read data from the original publication ("Population differentiation in allele frequencies of obesity-associated SNPs")
data <- read_excel("225SNPs26populations.xlsx")
colnames(data) <- data[2, ] # the second row will be the header
data <- data[-c(1,2), ]          # removing the first two rows
#The first three columns are SNP information. After that, each of the 26 populations has 5 columns:
#("population","effect allele number","other allele number","total allele number","effect allele frequency")
#The remaining columns (138-[3+26*5=133]=5) summarize global population data ("population","effect allele number","other allele number","effect allele frequency","GWAS P-value")

#Determine what populations to investigate
populations <- c("CEU","ASW","YRI")
cols <- c("population","effect allele number","other allele number","total allele number","effect allele frequency","GWAS P-value")

SNP.df <- data[,1:3]

MAF.df <- data.frame(row.names= SNP.df$`SNP ID`)

#Create a background "global population" (ALL.df) to use for hypergeometric testing. It is the sum of each column for the populations being interrogated.
ALL.df <- data.frame(matrix(nrow=nrow(data),ncol=length(cols)))
colnames(ALL.df) <- cols
ALL.df$population <- "ALL"
ALL.df$`GWAS P-value` <- data$`GWAS P-value`
ALL.df[is.na(ALL.df)] <- 0

dfs <- paste(populations,"df",sep=".")

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

#OPTION 1:
hyper.df <- data.frame()
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

#OPTION 2:
hyper.df <- data.frame()
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

#OPTION 3:
#Calculating the hypergeometric distributions using phyper().
#A random variable X has a hypergeometric distribution if X=x is the number of successes in a sample size of size k (without replacement) when the population contains M successes and N non-successes.
#For example, the probability of selecting 14 red marbles from a sample of 20 taken from an urn containing 70 red marbles and 30 green marbles is phyper(x=14,m=70,n=30,k=20).
#Thus, for each SNP, the parameters for this function are:
#q: number of effect alleles selected from the sample
#m: total number of effect alleles in the population
#n: total number of other alleles in the population
#k:	total number of alleles (effect + other) in the sample
#lower.tail=TRUE (depletion) or FALSE (enrichment)
enrichment.df <- data.frame()
depletion.df <- data.frame()
for(i in 1:nrow(ALL.df)){
  a <- c()
  b <- c()
  for (df in dfs){
    info <- get(df)
    m <- as.numeric(ALL.df$`effect allele number`[i]) #effect alleles in population
    n <- as.numeric(ALL.df$`other allele number`[i]) #other alleles in population
    k <- as.numeric(info$`total allele number`[i]) #total alleles in sample
    q <- as.numeric(info$`effect allele number`[i]) #effect alleles in sample
    a <- cbind(a,phyper(q=q-1,m=m,n=n,k=k,lower.tail=FALSE)) #When we set parameter lower.tail=FALSE in phyper the interpretation of the p-value is P[X > x] . But what we need to test is the null hypothesis P[X ≥ x], so we subtract x by 1.
    b <- cbind(b,phyper(q=q,m=m,n=n,k=k,lower.tail=TRUE))
  }
  enrichment.df <- rbind(enrichment.df,a)
  depletion.df <- rbind(depletion.df,b)
}

L <- c("enrichment.df","depletion.df")
#Assign new column and row names to each dataframe in "L"
for(df in L) {
  df.tmp <- get(df)
  colnames(df.tmp) <- populations
  rownames(df.tmp) <- SNP.df$`SNP ID`
  assign(df, df.tmp)
}

cutoff <- 0.01 / (2*3*225) #To control a family-wise error rate (FWER) of 0.01, we used a raw p-value of 7.41E-6 as cutoff. 
enrichedFiltered.df <- subset(enrichment.df < cutoff)
depleatedFiltered.df <- subset(depletion.df < cutoff)

enriched <- data.frame(enrichment.df < cutoff)
depleated <- data.frame(depletion.df < cutoff)

table(rowSums(enriched))
v1 <- ggvenn(enriched)
v1
table(rowSums(depleated))
v2 <- ggvenn(depleated)
v2


MAF.df["id"] <- rownames(MAF.df)
MAF_molten.df <- melt(MAF.df, id.vars="id", value.name="MAF", variable.name="Population")
p <- ggplot(MAF_molten.df,aes(x=id,y=MAF)) + geom_bar(position="identity", stat="identity") + facet_wrap(~Population,nrow=3)
p
q <- ggplot(MAF_molten.df,aes(x=id,y=MAF,fill=Population)) + geom_dotplot(binaxis='y', stackdir='center')
q
