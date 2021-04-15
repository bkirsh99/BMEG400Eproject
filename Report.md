BMEG 400E/591E Project
================

-   [Introduction](#introduction)
    -   [Overview of the original
        experiment](#overview-of-the-original-experiment)
    -   [Goals of the re-analysis](#goals-of-the-re-analysis)
-   [Methods](#methods)
    -   [Quality Control](#quality-control)
    -   [Analysis](#analysis)
-   [Results](#results)
-   [Conclusions](#conclusions)

# Introduction

## Overview of the original experiment

Large research efforts have been made to identify genomic regions,
variants in candidate genes, and environmental factors that contribute
to disease-associated health disparities across populations. In a study
by Mao et al., worldwide differences in the effect allele frequency of
225 obesity-associated SNPs were investigated. The original analysis was
conducted to identify significantly enriched or depleted effect alleles
across 26 populations, which were later interrogated with respect to the
relevance and underlying biological mechanisms of their associated
genes. Additionally, researchers inspected the population-based patterns
that emerged from differences in allele frequency and genetic risk
scores. The results suggest that over 85% of effect alleles exhibit
significant allele frequency differences and that population-level
differences in genetic risk scores are correlated with obesity
prevalence.

## Goals of the re-analysis

In GWAS, the over-representation of European populations is frequently
coupled with inadequate sampling of African genomes in SNP arrays and
imputation reference panels. As a result, ascertainment bias commonly
confounds the associations between genetic ancestry and disease
susceptibility. The main goal of this re-analysis is to explore
obesity-related variants among a subset of samples, namely those
associated with the African American (ASW), African (YRI), and European
American (CEU) populations. As such, this project will interrogate the
relationship between disease susceptibility, ancestry, and admixture to
elucidate the genetic underpinnings of health disparities across
populations. Furthermore, this re-analysis will evaluate the
transferability of Genome-Wide Association Studies (GWAS) by considering
the ancestral diversity within and between study cohorts. Ultimately,
this re-analysis will illustrate the spectrum of genetic diversity
between geographically distant genomes to reveal ancestral biological
signatures and systematic differences underlying the obesity epidemic.
In doing so, we hope to encourage diversity and inclusion in genomics
research, as well as urge for the investigation of genetic variants with
extreme allele frequency differences (EAFD) between populations.

``` bash
#
wget https://static-content.springer.com/esm/art%3A10.1186%2Fs12864-017-4262-9/MediaObjects/12864_2017_4262_MOESM1_ESM.xlsx
```

``` r
rm(list = ls()) # clean up the environment

################################################################################

# Load dependencies
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(tidyr)
```

    ## Warning: package 'tidyr' was built under R version 4.0.4

``` r
library(readxl)
library(reshape2)
```

    ## 
    ## Attaching package: 'reshape2'

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     smiths

``` r
library(purrr)
library(ggplot2)
library(ggvenn)
```

    ## Loading required package: grid

``` r
library(ggVennDiagram)
```

    ## Warning: package 'ggVennDiagram' was built under R version 4.0.4

``` r
library(heatmap3)
```

    ## Warning: package 'heatmap3' was built under R version 4.0.5

# Methods

## Quality Control

## Analysis

``` bash
```

# Results

# Conclusions
