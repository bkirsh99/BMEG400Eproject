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

#We calculate the significance of enrichment and depletion in the effect allele of each SNP using Fisher’s two-tail exact test (phyper in R).
#Then, we correct for multiple hypothesis testing using Bonferroni (FWER).
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

hyperRes <- data.frame()
L <- c("enrichment.df","depletion.df")
#Assign new column and row names to each dataframe in "L"
for(df in L) {
  df.tmp <- get(df)
  colnames(df.tmp) <- populations
  rownames(df.tmp) <- SNP.df$`SNP ID`
  df.tmp$type <- sub(".df","",df)
  df.tmp$`GWAS P-value` <- as.numeric(data$`GWAS P-value`[match(rownames(df.tmp),data$`SNP ID`)])
  assign(df, df.tmp)
  hyperRes <- rbind(hyperRes,df.tmp)
}

hyperRes_melted <- melt(hyperRes,id.vars=c(populations,"type","GWAS P-value"))

cutoff <- 0.01 / (2*3*225) #To control a family-wise error rate (FWER) of 0.01, we used a raw p-value of 7.41E-6 as cutoff. 
enrichedFiltered.df  <- data.frame(enrichment.df < cutoff)
depleatedFiltered.df <- data.frame(depletion.df < cutoff)

table(rowSums(enrichedFiltered.df ))
v1 <- ggvenn(enrichedFiltered.df )
v1
table(rowSums(depleatedFiltered.df))
v2 <- ggvenn(depleatedFiltered.df)
v2

#Combine enrichment and depletion p-scores into a single data structure, where values are log10 transformed.
#If the effect allele of an SNP is enriched in a population, then the negative of log10 of the enrichment p-value (a positive number) was used to represent the SNP in association with that population in a heatmap. On the other hand, if the allele of an SNP is depleted in a population, the value of log10 of the depletion p-value (a negative number) was used to represent the SNP for that population in the heatmap.
keep <- which(names(hyperRes_melted) %in% c("type", "GWAS P-value"))
hyperRes_sig <- hyperRes_melted[!rowSums(hyperRes_melted[,-keep] < cutoff)==0,,drop=FALSE]
hyperRes_log <- hyperRes_sig %>% mutate_at(-keep, log10)
hyperRes_log0 <- hyperRes_log
hyperRes_log[hyperRes_log$type == "enrichment",] <- hyperRes_log[hyperRes_log$type == "enrichment",] %>% mutate_at(-keep, ~ . * -1)

library(heatmap3)
log_enriched <- enrichment.df
log_enriched$CEU <- -log10(log_enriched$CEU)
log_enriched$ASW <- -log10(log_enriched$ASW)
log_enriched$YRI <- -log10(log_enriched$YRI)

enriched_hm <- heatmap(as.matrix(log_enriched), Colv = NA, Rowv = NA, hclustfun = function(x) hclust(x,method = 'centroid'), scale = "row")
enriched_hm

#Create a matrix from the FWER filtered, combined, log-transformed, and additive inverse-corrected p-values. Then, generate a clustermap.
hyperRes_matrix <- as.matrix(hyperRes_log[ , -keep])
col_fun = colorRampPalette(c("green", "black", "red"))(1024)
enriched_hm2 <- heatmap3(hyperRes_matrix, Colv = NA, Rowv = NA, hclustfun = function(x) hclust(x,method = 'centroid'), col = col_fun, scale = "row")
enriched_hm2

#Filter again for for a selected set of SNPs that have enrichment or depletion p-values of at least 10−E100 and have reached genome-wide significance (5 × 10−8) in GWA studies.
hyperRes_sig2 <- hyperRes_log[!rowSums(hyperRes_log[,-keep] < 10E-100)==0, , drop = FALSE]
hyperRes_sig2 <- hyperRes_sig2[hyperRes_sig2$`GWAS P-value` < 5E-08,]
filteredSNPs <- rownames(hyperRes_sig2)

#if in filteredSNPs
MAF.df["id"] <- rownames(MAF.df)
MAF_molten.df <- melt(MAF.df, id.vars="id", value.name="MAF", variable.name="Population")
MAF_filtered.df <- MAF_molten.df[MAF_molten.df$id %in% filteredSNPs ,]
# p <- ggplot(MAF_molten.df,aes(x=id,y=MAF)) + geom_bar(position="identity", stat="identity") + facet_wrap(~Population,nrow=3)
p
q <- ggplot(MAF_filtered.df,aes(x=id,y=MAF,fill=Population)) + geom_dotplot(binaxis='y', stackdir='center')
q
