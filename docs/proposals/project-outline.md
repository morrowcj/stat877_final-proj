---
title: "Stat 877 Project"
authors: "Clay Morrow & Ting-Fung Ma"
date: 04-March-2019
---

# Problem:
In the field of community genetics, the concept of "extended phenotype" has been
used to identify elements of a biological community that are indirectly affected
by the genetics of another organism. This project aims to identify genetic
regions of DNA, from a non-model tree species, that indirectly influence associated
herbivorous insect communities.

# System Description:
Our project will use data collected from a common-garden experimental plantation
of the clonal tree species *Populus tremuloides* (aspen).
The plot contains 1600 individual trees, each from one of 500 genetically
distinct clones (genets).

## Insect Data (Response)
In both 2016 and 2017, 2 summer insect surveys were
conducted in which insect presence and abundance were measured for 18 common
species of herbivorous insects on each tree.  

## Genetic Data (Predictor)
DNA sequencing was done for each aspen genet and over 5000 coding region SNPs
were identified for the population. Data for major and minor alleles for each SNP
are available, as well as each genet's genotype at each SNP site.

## Tree Traits (Co-variate Predictors)
Tree trait data, for 10 traits, corresponding to the individual trees during
the insect survey periods are also available. Such traits include tree size
(height, diameter, volume), phytochemical concentrations (foliar N, secondary
  metabolites), physical leaf traits, phenological traits, and gender.

## Environmental Factors (more Co-variate predictors)  
There are also a number of environmental variables to consider (weather, soil).

# Goal
The goal of this project is to investigate the relationships between genotype,
at various SNP regions, and insect presence and/or abundance and to identify
potential quantitative trait loci (QTLs) associated with any of the common
insect species. Marker regression will be used to identify QTLs from the SNP
genotype markers (possibly a subset w/even distribution of genotypes)
and any significant covariates. This will likely be a computationally expensive
process so only 1 insect species will be considered as a response at first. If
time permits, we may move to a more complex multi-variate model.

# Responsibilities
- Clay: Will provide the data, insight into the system, and contribute to
building statistical models and R code.

- Ting-Fung: Will provide statistical expertise and will lead in the design
of statistical models and code for computation.
