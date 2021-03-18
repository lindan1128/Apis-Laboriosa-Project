---
title: "Calculate genetic distances"
author: "LIN Dan"
date: "12/3/2020"
---
  
# load packages
library(optparse)
library(ape)
library(phangorn)

# create parser
option_list <- list(
  make_option(c("-d", "--gene_family_directory"),type="character",help="Path to gene families folder. Please provide the absolute path.")
)
args <- parse_args(OptionParser(option_list=option_list))

# set working directory
setwd(args$gene_family_directory)

# go through all gene familyfiles
for (i in c(0:10640)){
  files = paste(i,'.conserved_gene_family_for_genetic_distance.phy', sep = '')
  if (file.exists(files) == TRUE){
    # calcute p distance
    ex.dna <- read.dna(files, format = "sequential", as.character=TRUE)
    dat <- phyDat(ex.dna, "USER", levels=unique(as.vector(ex.dna)))
    p_dist <- dist.p(dat)[1]
    len <- length(ex.dna)/2
    cat(p_dist, file = "p_distance.txt", append = T)
    cat('\t', file = "p_distance.txt", append = T)
    cat(len, file = "p_distance.txt", append = T)
    cat('\n', file = "p_distance.txt", append = T)
    
    # calculate genetic distance
    dat <- as.DNAbin(ex.dna)
    jc_dist <- dist.dna(dat, model = "JC69", gamma = TRUE)[1]
    k_dist <- dist.dna(dat, model = "K80", gamma = TRUE)[1]
    f_dist<- dist.dna(dat, model = "F84", gamma = TRUE)[1]
    logdet_dist <- dist.dna(dat, model = "logdet", gamma = TRUE)[1]
    cat(jc_dist, file = "genetic_distance.txt", append = T)
    cat('\t', file = "genetic_distance.txt", append = T)
    cat(k_dist, file = "genetic_distance.txt", append = T)
    cat('\t', file = "genetic_distance.txt", append = T)
    cat(f_dist, file = "genetic_distance.txt", append = T)
    cat('\t', file = "genetic_distance.txt", append = T)
    cat(logdet_dist, file = "genetic_distance.txt", append = T)
    cat('\n', file = "genetic_distance.txt", append = T)
  }
}
    