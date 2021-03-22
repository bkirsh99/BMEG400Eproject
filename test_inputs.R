library(dplyr)
library(tidyr)
library(sim1000G)
library(vcfR)

dir <- "C:/Users/bkirs/Documents/School/BIT/gwas-trial"
setwd(dir)

#Import GWAS data:
associations <- read.csv("coronary_heart_disease-associations.csv",header=TRUE)
studies <- read.csv("coronary_heart_disease-studies.csv",header=TRUE)
colnames(associations) <- c("variant_and_risk_allele","p-value","p-value_annotation",
                            "risk_allele_frequency","odds_ratio","effect_size", 
                            "confidence_interval", "mapped_gene","reported_trait","traits",
                            "study_accession","location")
colnames(studies) <- c("first_author","study_accession","publication_date","journal",
                       "title","reported_trait","traits","discovery_sample_number_and_ancestry",
                       "replication_sample_number_and_ancestry","association_count",
                       "summary_statistics")

#Import 1000 Genomes data:
vcf <- read.vcfR("ALL.chr1.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz")
vcf_field_names(vcf)
queryMETA(vcf)
gt <- extract.gt(vcf)
info <- getINFO(vcf)
fix <- getFIX(vcf)
fix.df <- as.data.frame(fix)
fix.df <- mutate(fix.df,location=paste(CHROM,POS,sep=":"))
fix.df$POS <- as.character(fix.df$POS)
Z <- vcfR2tidy(vcf) 

#Define variables to use:
chromosomes <- data.frame("chr"=c(seq(1,22),"X","Y","XY","MT"),"code"=seq(1,26))

#Define data frames to work with:
#1) Only information about GWAS studies:
summary <- as.data.frame(distinct(studies,study_accession,.keep_all=TRUE)[,"study_accession"])
#eliminate duplicate entries for the same study 

#2) Only information about GWAS associations:
data <- associations

#3) Information about GWAS associations and their corresponding studies:
#data <- merge(associations,studies) 

#Work with data frames:
#1) Extract relevant columns to investigate the populations sampled in each study: 
samples <- distinct(studies,study_accession,.keep_all=TRUE)[,c("discovery_sample_number_and_ancestry",
                                                               "replication_sample_number_and_ancestry")]
#NOTE: Possible formats of the "sample_number_and_ancestry" columns:
  #I) Single population: "9618 European"
  #II) Multiple populations with different sample sizes:
    #II.a) Comma-separated:
      #II.a.1) "18 African unspecified, Asian unspecified1794 European"
      #II.a.2) "340799 East Asian, European51442 East Asian"
    #II.b) Non-comma-separated:
      #II.b.1) "3197 NR3557 European"
      #II.b.2) 11323 East Asian25557 South Asian2268 Greater Middle Eastern (Middle Eastern, North African or Persian)4095 Hispanic
  #III) No record: "'-" or "NR"

#2) Retain only ancestry information, regardless of sample numbers: 
discovery_ancestry<- strsplit(gsub("[[:digit:]]",",",samples$discovery_sample_number_and_ancestry), "[,]+")
replication_ancestry<- strsplit(gsub("[[:digit:]]",",",samples$replication_sample_number_and_ancestry), "[,]+")
study_ancestry <- Map(c,discovery_ancestry,replication_ancestry)
  #combine discovery and replication populations 
populations <- list() 
  #initialize empty list to store all possible populations

#3) Create a "count" matrix of sampled populations for each study. Here, a 0 means that the population was sampled:
for (i in 1:nrow(summary)){
  vector <- study_ancestry[i]
  ancestries <- unlist(vector)
  ancestries <- trimws(unique(ancestries[ancestries != ""],"l"))
  for (j in 1:length(ancestries)){
    population <- ancestries[j]
    if (population == "'-"){
      population <- "NR"
    }
    summary[i,population] <- 0
  }
}

#4) Convert p-values from string format (5 x 10-9) to numeric values in scientific notation (5e-09):
data$`p-value` <- strsplit(gsub("x[^x]*-","-",gsub(" ","",data$`p-value`)),"-")
for (i in 1:nrow(data)){
  pval <- data$`p-value`[i]
  base <- as.numeric(unlist(pval)[1])
  exponent <- as.numeric(unlist(pval)[2])
  data$`p-value`[i] <- base * 10 ** (-1 * exponent) 
}

#5) Remove SNP-SNP interactions and multi-SNP haplotypes:
data <- data[!grepl("\\|", data$location),]
  #this analysis is simplified by excluding SNP=SNP interactions and multi-SNP haplotypes
  #however, they can be included by making the following formatting changes:
#data <- separate_rows(data,location,sep="\\|")
#data <- separate_rows(data,variant_and_risk_allele,sep="([x,])+")

#6) Split the "variant_and_risk_allele" column into two columns, one with the SNP id and the other with the effect allele:
#NOTE: Possible formats of the "variant_and_risk_allele" column:
  #I) SNP: rs140607780-<b>G</b>
  #II) SNP-SNP interaction: rs1165668-<b>G</b> x rs1165669-<b>C</b> (mapped gene separated by "x" and location separated by "|")
  #III) Multi-SNP haplotypes: rs7767084-<b>T</b>, rs10755578-<b>G</b>, rs3127599-<b>T</b>, rs2048327-<b>C</b> (mapped gene separated by ";" and location separated by "|")
data$variant_and_risk_allele <- gsub("<[^>]*>","",data$variant_and_risk_allele)
data <- separate(data,variant_and_risk_allele,c("variant","risk_allele"),"-")

#7) Subset data frame to retain only interesting information:
#basic <- c("location", "variant","risk_allele","p-value","risk_allele_frequency","odds_ratio",
#          "effect_size","confidence_interval","mapped_gene","study_accession")
extra <- c("p-value_annotation","reported_trait","traits")
keep <- data[,!(names(data) %in% extra)]

#7.a) Create ".map" file (chr, variant, pos, coord) for PLINK applications:
mapCol <- c("location","variant")
map.df <- data[,names(data) %in% mapCol]
map.df <- separate(map.df,location,c("chr","coord"),":")
map.df$pos <- 0
map.df <- drop_na(map.df[c("chr","variant","pos","coord")])
map.df$chr <- chromosomes$code[match(map.df$chr,chromosomes$chr)]

table(map.df$chr) #chr16 contains the largest number of CHD-associated SNPs (209)

#7. b) Create a ".txt" file for filtering based on snp ids:
snps <- write.table(data$variant, file = "snp_ids.txt", quote = FALSE, row.names = FALSE)

#8) Filter the SNPs and remove duplicate entries:
#NOTE: Pre-processing criteria:
  #I) reported risk allele
  #II) p-value < 9e−06
  #III) odds ratio > 1.0
filter <- keep[!(keep$risk_allele=="?" | keep$`p-value`>9e-06 | keep$odds_ratio<1.0),]
#byvar <- filter %>% group_by(variant, risk_allele)

info_snps <- merge(fix.df,filter,by=location)

#In GWAS, significant SNPs are always reported as the risk allele (RA), or effective allele, and can be presented either on the forward or the reverse strand of the genome.
#Typically, the reported RA for a SNP will be converted to the forward strand prior to publication.
#Data in the GWAS Catalog is currently mapped to genome assembly GRCh38.p13 and dbSNP Build 153.
