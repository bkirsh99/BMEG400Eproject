library(dplyr)
library(tidyr)
library(readxl)
library(reshape2)
library(purrr)
library(ggplot2)
library(ggvenn)
library(ggVennDiagram)
library(heatmap3)

dir <- "C:/Users/bkirs/Documents/School/BMEG400E/BMEG400Eproject"
setwd(dir)

data <- read_excel("225SNPs26populations.xlsx") # read the data from the original publication
colnames(data) <- data[2, ] # the second row is now the header
data <- data[-c(1,2), ] # and the first two rows can be removed
SNP_15 <- c("rs6499640","rs11075990","rs9939609","rs7202116","rs7185735","rs9940128","rs1121980","rs17817449","rs8043757","rs8050136","rs1421085","rs1558902","rs12149832","rs62033400","rs17817964")
data <- data[data$`SNP ID` %in% SNP_15,,drop = FALSE] #filtering by SNPs of interest
# the first 3 columns contain general SNP information ("SNP ID", "effect allele", and "other allele")
# the remaining 130 columns contain population-specific information
# each of the 26 populations takes up 5 subsequent columns ("population", "effect allele number", "other allele number", "total allele number", and "effect allele frequency")
# the final 5 columns (138-[3+26*5=133]=5) summarize global population data ("population", "effect allele number", "other allele number", "effect allele frequency", and "GWAS P-value")

populations <- c("CEU","ASW","YRI") # determine what populations to investigate
cols <- c("population","effect allele number","other allele number","total allele number",
          "effect allele frequency") # determine relevant parameters to report in the data frames created from here on

varInfo <- data[,1:3] # this data frame is a subset of the original data and contains only general SNP information
popInfo <- data.frame(row.names= data$`SNP ID`) # this data frame is a subset of the original data and contains only population-specific information
bkgdInfo <- data.frame(matrix(nrow=nrow(data),ncol=length(cols))) # this data frame is a subset of the original data and contains only global population information                                            

colnames(bkgdInfo) <- cols
rownames(bkgdInfo) <- data$`SNP ID`
bkgdInfo$population <- "ALL"
bkgdInfo$`GWAS P-value` <- data$`GWAS P-value`
bkgdInfo[is.na(bkgdInfo)] <- 0

hyperRes <- data.frame() # this data frame holds the results of the enrichment analysis

popDfs <- paste(populations,"Info",sep="") # list of data frames to be created (one per user-defined population)
hyperDfs <- c("enrDf","depDf") # list of data frames to be created (one per hypergeometric test: enrichment or depletion)

renameDfs <- function(dfs){ # this function re-names data frames, modifying rows to SNP IDs and columns to field names containing a space
  # argument: dfs: a list of data frames (strings) to be modified
  # return: NA (modifies the data frames in "dfs" in the global environment)
  for(df in dfs) {
    df.tmp <- data.frame(get(df))
    rownames(df.tmp) <- data$`SNP ID`
    colnames(df.tmp) <- cols
    assign(df, df.tmp, envir = .GlobalEnv)
  }
}

renameColsRows <- function(dfs){ # this function re-names data frames, modifying rows to SNP IDs and columns to population IDs 
  # argument: dfs: a list of data frames (strings) to be modified
  # return: NA (modifies the data frames in "dfs" in the global environment)
  for(df in dfs) {
    df.tmp <- data.frame(get(df))
    colnames(df.tmp) <- populations
    rownames(df.tmp) <- data$`SNP ID`
    assign(df, df.tmp, envir = .GlobalEnv)
  }
}

populateDfs <- function(pop){ # this function populates data frames using the original data for user-specified populations
  # argument: pop: a list of user-defined populations to analyze
  # return: NA (modifies the "popInfo" and "bkgdInfo" data frames in the global environment)
  
  for(i in seq(from=4, to=ncol(data), by=5)){
    population <- toString(data[1,i])
    if(population %in% populations){
      popInfo[population] <<- data[,i+4]
      assign(paste(population,"Info",sep=""), data[,i:(i+4)], envir = .GlobalEnv)
      bkgdInfo["effect allele number"] <<- as.numeric(bkgdInfo[["effect allele number"]]) + as.numeric(data[["effect allele number"]])
      bkgdInfo["other allele number"] <<- as.numeric(bkgdInfo[["other allele number"]]) + as.numeric(data[["other allele number"]])
      bkgdInfo["total allele number"] <<- as.numeric(bkgdInfo[["total allele number"]]) + as.numeric(data[["total allele number"]])
      bkgdInfo["effect allele frequency"] <<- as.numeric(bkgdInfo[["effect allele frequency"]]) + as.numeric(data[["effect allele frequency"]])
    }
  }
  
  renameDfs(popDfs)
  
}

populateDfs(populations) # run the pipeline for the set of user-defined populations

combineHyper <- function(dfs){ # this function combines two data frames, namely the ones produced by hyperTest (enrichment and depletion)
  for(df in dfs) {
    df.tmp <- get(df)
    colnames(df.tmp) <- populations
    rownames(df.tmp) <- data$`SNP ID`
    df.tmp$type <- sub("Df","",df)
    df.tmp$`GWAS P-value` <- as.numeric(data$`GWAS P-value`[match(rownames(df.tmp),data$`SNP ID`)])
    assign(df, df.tmp)
    hyperRes <- rbind(hyperRes,df.tmp)
  }
  return(hyperRes)
  
}

hyperTest <- function(pop){ # this function performs the enrichment analysis using Fisher’s two-tail exact test
  # two hypergeometric tests, one for enrichment and another for depletion, are conducted per SNP
  # argument: pop: a list of user-defined populations to analyze
  # return: df1, df2: the data frames containing the results of the enrichment analysis (one for enrichment and another for depletion)
  
  enrichmentDf <- data.frame()
  depletionDf <- data.frame()
  
  for(i in 1:nrow(bkgdInfo)){
    a <- c()
    b <- c()
    for (df in popDfs){
      info <- get(df)
      m <- as.numeric(bkgdInfo$`effect allele number`[i]) # total number of effect alleles in the population
      n <- as.numeric(bkgdInfo$`other allele number`[i]) # total number of other alleles in the population
      k <- as.numeric(info$`total allele number`[i]) # total number of alleles (effect + other) in the sample
      q <- as.numeric(info$`effect allele number`[i]) # number of effect alleles selected from the sample
      
      # phyper(q,m,n,k) = phyper(success-in-sample, success-in-bkgd, failure-in-bkgd, sample-size)
      a <- cbind(a,phyper(q=q-1,m=m,n=n,k=k,lower.tail=FALSE)) # over-representation (enrichment)
      b <- cbind(b,phyper(q=q,m=m,n=n,k=k,lower.tail=TRUE)) # under-representation (depletion)
    }
    enrichmentDf <- rbind(enrichmentDf,a)
    depletionDf <- rbind(depletionDf,b)
  }
  return(list(enrichmentDf,depletionDf))
}

testRes <- hyperTest(populations)
enrDf <- testRes[[1]]
depDf <- testRes[[2]]

renameColsRows(hyperDfs)
hyperRes <- combineHyper(hyperDfs)

makeVenn <- function(df,pop,sig){ # this function creates a Venn diagram that illustrates the differences and similarities between the chosen populations
                              # arguments: df: a data frame of p-values from the enrichment analysis
                              #            pop: the populations to compare
                              #            sig: the significance level (alpha)
                              # return: a Venn diagram
  
  cutoff <- sig / (2*length(pop)*nrow(data)) # control the family-wise error rate (FWER) from MHT using Bonferroni correction
  # in the case of 3 populations and 225 SNPs, we use a raw p-value of 7.41E-6 as the cutoff for alpha = 0.01
  
  filteredDf <- data.frame(df < cutoff)
  summary <- table(rowSums(filteredDf))
  plot <- ggvenn(filteredDf)
  return(plot)
}

enrVenn <- makeVenn(enrDf,populations,0.01)
depVenn <- makeVenn(depDf,populations,0.01)

makePlots <- function(df,pop,sig){
  #dev.off()
  cutoff <- sig / (2*length(pop)*nrow(data))
  
  meltedDf <- melt(df,id.vars=c(pop,"type","GWAS P-value"))
  keep <- which(names(meltedDf) %in% c("type", "GWAS P-value")) # columns that should not be modified by mutate functions
  sigDf <- meltedDf[!rowSums(meltedDf[,-keep] < cutoff)==0,,drop=FALSE] # remove rows whose raw p-value is not significant at alpha 
  logDf <- sigDf %>% mutate_at(-keep, log10) # log10 transform p-values
  heatmapDf <- logDf
  heatmapDf[heatmapDf$type == "enr",] <- heatmapDf[heatmapDf$type == "enr",] %>% mutate_at(-keep, ~ . * -1)
  # the negative of log10 of the p-value (a positive number) is used to represent enriched effect allele of a SNP for a population in the heatmap
  # the actual value of log10 of the p-value (a negative number) is used to represent depleted effect allele of a SNP for a population in the heatmap
  
  heatmapMtx <- as.matrix(heatmapDf[ , -keep]) # create a matrix from the FWER filtered, combined, log-transformed, and additive inverse-corrected p-values
  col_fun = colorRampPalette(c("green", "black", "red"))(1024)
  clusterHm <- heatmap3(heatmapMtx, Colv = NA, Rowv = NA, hclustfun = function(x) hclust(x,method = 'centroid'), col = col_fun, scale = "row")
  
  # filter again for for a selected set of SNPs that have enrichment or depletion p-values of at least 10−E100 and have reached genome-wide significance (5 × 10−8) in GWA studies.
  freqDf <- logDf[!rowSums(abs(logDf[,-keep]) < 10E-100)==0,,drop=FALSE] # remove rows whose log10 transformed p-value is not significant at 10E-100
  freqDf <- freqDf[freqDf$`GWAS P-value` < 5E-08,]
  filteredSNP <- rownames(freqDf)
  
  popInfo["id"] <- rownames(popInfo)
  meltedPopDf <- melt(popInfo, id.vars="id", value.name="MAF", variable.name="Population")
  filteredPopDf <- meltedPopDf[meltedPopDf$id %in% filteredSNP ,]
  freqDiff <- ggplot(filteredPopDf,aes(x=id,y=MAF,fill=Population)) + geom_dotplot(binaxis='y', stackdir='center')
  return(list(clusterHm,freqDiff))
}

plotRes <- makePlots(hyperRes,populations,0.01)
