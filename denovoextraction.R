#!/usr/bin/env Rscript 

#Load command arguments
message("Start")
message("Reading arguments")
input <- commandArgs(trailingOnly = TRUE)

library_location <- input[1]
combined_mut_mat_location <- input[2]
perform_estimate <- input[3]
estimate_location <- input[4]
perform_NMF <- input[5]
nmf_location <- input[6]
number_signatures <- as.numeric(input[7])
print(number_signatures)

#Load packages
library(rngtools, lib.loc=library_location)
library(registry, lib.loc=library_location)
library(pkgmaker, lib.loc=library_location)
library(NMF, lib.loc=library_location)
library(MutationalPatterns, lib.loc=library_location)
library(tidyverse, lib.loc=library_location)
library(cowplot, lib.loc=library_location)

#Load combined mutation matrix
assign("mutation_matrix", get(load(combined_mut_mat_location)))

#Perform NMF estimate
mutation_matrix <- mutation_matrix + 0.0001
if (perform_estimate){
  print("Starting estimation")
  estimate <- NMF::nmf(mutation_matrix, rank=2:15, method="brunet", nrun=10, seed=123456)
  pdf(estimate_location)
  print(plot(estimate))
  dev.off()
}

#Perform NMF
if (perform_NMF){
  print("Starting NMF")
  nmf_res <- extract_signatures(mutation_matrix,
                                rank = number_signatures,
                                nrun = 200,
				single_core = TRUE)
  save(nmf_res,
     file = nmf_location)
}

message("Ready")
