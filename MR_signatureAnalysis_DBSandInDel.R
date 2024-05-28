#!/usr/bin/env Rscript 

# Clear data and set seed
rm(list = ls())
set.seed(1234)

# Load libraries
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)
library("MutationalPatterns")
library("ggplot2")
library("logging")
library("bedr")
library("ensemblVEP")
library("TVTB")
library("readxl")
library("stringr")
library("gtools")
library("tibble")
library("dplyr")
library("openxlsx")

# Load custom functions
FUNCTIONS_SCRIPT <- list.files(path = "/path/to/code/dir/", 
                               pattern = "createdFunctions_HM-ALL.R",
                               full.names = TRUE, recursive = TRUE)
source(FUNCTIONS_SCRIPT)

# Set global parameters
# All input files should be located somewhere in PROJECT_DIR
PROJECT_DIR <- "/path/to/project/dir/"
RESULTS_DIR <- paste0(PROJECT_DIR, "/path/to/result/dir/")
RANDOM_ID_INFO <- "MR_cohortDescription.xlsx"
CENTROMERES_FILE <- "cytoBand.hg38.centromeresOnly.txt"
MNP_EXT <- "-merged-mutect2Calls.passFiltered.mnp.vepAnnotated.germlinefiltered.vcf.gz"
INDEL_EXT <- "-merged-mutect2Calls.passFiltered.indel-WGS.vepAnnotated.vcf.gz"
if (!dir.exists(RESULTS_DIR)){dir.create(RESULTS_DIR)}
setwd(RESULTS_DIR)

# Read random ID info 
random_id_info <- as.data.frame(read_excel(list.files(path = PROJECT_DIR,
                                                      pattern = RANDOM_ID_INFO,
                                                      recursive = TRUE,
                                                      full.names = TRUE))[1:29, c(1,2)])
colnames(random_id_info) <- c("skion_nr", "random_id")

# List input files
centromeres <- list.files(path = PROJECT_DIR,
                          pattern = CENTROMERES_FILE,
                          full.names = TRUE,
                          recursive = TRUE)
mnp_files <- list.files(path = PROJECT_DIR,
                        pattern = MNP_EXT,
                        recursive = TRUE,
                        full.names = TRUE)
indel_files <- list.files(path = PROJECT_DIR,
                          pattern = INDEL_EXT,
                          recursive = TRUE,
                          full.names = TRUE)

# Load COSMIC signatures
dbs_signatures <- get_known_signatures(muttype = "dbs")
indel_signatures <- get_known_signatures(muttype = "indel")



#------------------------------ Variant filtering -----------------------------#

## DBS filtering
# Define sample IDs
sample_ids <- gsub(MNP_EXT, "", basename(mnp_files))
loginfo(paste("found:", length(sample_ids), "MNP VCF files"))

# Match sample IDs to random IDs
random_ids <- random_id_info$random_id[match(sample_ids, random_id_info$skion_nr)]

# Filter each MNP file
dbs_mut_mat <- matrix(nrow = 78)
for (i in (1:length(sample_ids))){
  random_id <- random_ids[i]
  sample_vcf <- loadVcfFile(mnp_files[i],
                            ref_genome,
                            "chromosomes")
  GenomeInfoDb::genome(sample_vcf)  <- "hg38"

  # Rename colnames
  for (sample_id in colnames(sample_vcf)){
    if (endsWith(sample_id, sample_ids[i])){
      time <- str_split(sample_id, "-")[[1]][1]
      new_name <- paste0(sample_ids[i], time)
    } else if (startsWith(sample_id, "P_")){
      time <- str_split(gsub("_F$", "_F1", sample_id), "_")[[1]][3]
      new_name <- paste0(sample_ids[i], time)
    } else {
      new_name <- sample_id
    }
    if (endsWith(new_name, "D")){
      new_name <- paste0(new_name, "x")
    }
    colnames(sample_vcf)[colnames(sample_vcf) == sample_id] <- new_name
  }

  ## Filter VCF
  # Remove variants in centromeric regions
  sample_vcf_no_centromeric_variants <- excludeVariantsInCentromericRegions(sample_vcf, centromeres)
  
  # Apply popMax filters
  sample_vcf_gnomad_filtered <- filterOnGnomadAlleleFrequency(sample_vcf_no_centromeric_variants,
                                                              0.01)
  sample_vcf_gonl_filtered <- filterOnGonlAlleleFrequency(sample_vcf_gnomad_filtered,
                                                          0.01)
  sample_gonl_filtered_df <- convertVcfObjectToDataFrame(sample_vcf_gonl_filtered)
  
  #Look at separate timepoints
  timepoints <- mixedsort(unique(grep("D|([^C]R)", colnames(sample_vcf_gonl_filtered), value = TRUE)))
  timepoints <- timepoints[!grepl("^X", timepoints)]

  # Filter each time point seperately
  for (timepoint in timepoints){
    tumor_time <- gsub(sample_ids[i], "", timepoint)
    tumor_info <- paste(random_id, tumor_time, sep = "_")
    
    # Find ref, alt and vaf columns for each time point
    read_count_ref_idx <- tail(which(names(sample_gonl_filtered_df) %in% 
                                  paste0(timepoint, "_read_count_ref")), n = 1)
    read_count_alt_idx <- tail(which(names(sample_gonl_filtered_df) %in% 
                                  paste0(timepoint, "_read_count_alt")), n = 1)
    AF_idx <- tail(which(names(sample_gonl_filtered_df) %in% 
                      paste0(timepoint, "_read_based_AF")), n = 1)
    
    # Generate depth column
    sample_gonl_filtered_df[paste0(timepoint, "_total_read_count")] <- 
      sample_gonl_filtered_df[read_count_ref_idx] + 
      sample_gonl_filtered_df[read_count_alt_idx]
    total_read_count_idx <- which(names(sample_gonl_filtered_df) %in% 
                                    paste0(timepoint, "_total_read_count"))
    
    # Apply filters. Use custom AF filters per sample when necessary. 
    # Low tumor purity could result in lower AF values in general. 
    # For these samples the standard threshold might be too high.
    if (random_id == "P0105" & tumor_time == "R2") {
      AF_threshold <- 0.13
    } else if (random_id == "P0121") {
      custom_AFs <- c(0.25, 0.25, 0.16, 0.13)
      names(custom_AFs) <- c("Dx", "R1", "R2", "R3")
      AF_threshold <- as.numeric(custom_AFs[tumor_time])
    } else if (random_id == "P0140") {
      custom_AFs <- c(0.22, 0.24, 0.15)
      names(custom_AFs) <- c("Dx", "R1", "R3")
      AF_threshold <- as.numeric(custom_AFs[tumor_time])
    } else if (random_id == "P0160") {
      AF_threshold <- ifelse(tumor_time == "R2", 0.09, 0.21)
    } else if (random_id == "P0164") {
      custom_AFs <- c(0.25, 0.14, 0.2)
      names(custom_AFs) <- c("Dx", "R1", "R2")
      AF_threshold <- as.numeric(custom_AFs[tumor_time])
    } else if (random_id == "P0180") {
      custom_AFs <- c(0.25, 0.25, 0.08, 0.18)
      names(custom_AFs) <- c("Dx", "R1", "R2", "R3")
      AF_threshold <- as.numeric(custom_AFs[tumor_time])
    } else {
      AF_threshold <- 0.25
    }
    df_read_based_filtered <- sample_gonl_filtered_df[sample_gonl_filtered_df[,total_read_count_idx] >= 20 &
                                                        sample_gonl_filtered_df[,read_count_alt_idx] >= 5 &
                                                        sample_gonl_filtered_df[,AF_idx] >= AF_threshold,]
    
    # Subset VCF for filtered variants
    sample_filtered_vcf_object <- sample_vcf[rownames(sample_vcf) %in% 
                                               rownames(df_read_based_filtered), ]

    # Determine mutation types and obtain DBS contexts
    sample_dbs_grl <- get_mut_type(rowRanges(sample_filtered_vcf_object), 
                                   "dbs", 
                                   predefined_dbs_mbs = TRUE)
    sample_dbs_grl <- get_dbs_context(sample_dbs_grl)
    
    # Create mutation matrix and add to object
    sample_mut_mat <- count_dbs_contexts(sample_dbs_grl)
    colnames(sample_mut_mat) <- tumor_info
    dbs_mut_mat <- cbind(dbs_mut_mat, sample_mut_mat)
  }
}
dbs_mut_mat <- dbs_mut_mat[, 2:ncol(dbs_mut_mat)]
colSums(dbs_mut_mat)



## InDel filtering
# Define sample IDs
sample_ids <- gsub(INDEL_EXT, "", basename(indel_files))
loginfo(paste("found:", length(sample_ids), "InDel VCF files"))

# Match sample IDs to random IDs
random_ids <- random_id_info$random_id[match(sample_ids, random_id_info$skion_nr)]

# Filter each InDel file
indel_mut_mat <- matrix(nrow = 83)
for (i in (1:length(sample_ids))){
  random_id <- random_ids[i]
  sample_vcf <- loadVcfFile(indel_files[i],
                            ref_genome,
                            "chromosomes")
  GenomeInfoDb::genome(sample_vcf)  <- "hg38"
  
  # Rename colnames
  for (sample_id in colnames(sample_vcf)){
    if (endsWith(sample_id, sample_ids[i])){
      time <- str_split(sample_id, "-")[[1]][1]
      new_name <- paste0(sample_ids[i], time)
    } else if (startsWith(sample_id, "P_")){
      time <- str_split(sample_id, "_")[[1]][3]
      new_name <- paste0(sample_ids[i], time)
    } else {
      new_name <- sample_id
    }
    if (endsWith(new_name, "D")){
      new_name <- paste0(new_name, "x")
    }
    colnames(sample_vcf)[colnames(sample_vcf) == sample_id] <- new_name
  }
  
  ## Filter VCF
  # Remove variants in centromeric regions
  sample_vcf_no_centromeric_variants <- excludeVariantsInCentromericRegions(sample_vcf, centromeres)
  
  # Apply popMax filters
  sample_vcf_gnomad_filtered <- filterOnGnomadAlleleFrequency(sample_vcf_no_centromeric_variants,
                                                              0.01)
  sample_vcf_gonl_filtered <- filterOnGonlAlleleFrequency(sample_vcf_gnomad_filtered,
                                                          0.01)
  sample_gonl_filtered_df <- convertVcfObjectToDataFrame(sample_vcf_gonl_filtered)
  
  #Look at separate timepoints
  timepoints <- mixedsort(unique(grep("D|([^C]R)", colnames(sample_vcf_gonl_filtered), value = TRUE)))
  timepoints <- timepoints[!grepl("^X", timepoints)]

  # Filter each time point seperately
  for (timepoint in timepoints){
    tumor_time <- gsub(sample_ids[i], "", timepoint)
    tumor_info <- paste(random_id, tumor_time, sep = "_")
    
    # Find ref, alt and vaf columns for each time point
    read_count_ref_idx <- tail(which(names(sample_gonl_filtered_df) %in% 
                                       paste0(timepoint, "_read_count_ref")), n = 1)
    read_count_alt_idx <- tail(which(names(sample_gonl_filtered_df) %in% 
                                       paste0(timepoint, "_read_count_alt")), n = 1)
    AF_idx <- tail(which(names(sample_gonl_filtered_df) %in% 
                           paste0(timepoint, "_read_based_AF")), n = 1)
    
    # Generate depth column
    sample_gonl_filtered_df[paste0(timepoint, "_total_read_count")] <- 
      sample_gonl_filtered_df[read_count_ref_idx] + 
      sample_gonl_filtered_df[read_count_alt_idx]
    total_read_count_idx <- which(names(sample_gonl_filtered_df) %in% 
                                    paste0(timepoint, "_total_read_count"))
    
    # Apply filters. Use custom AF filters per sample when necessary. 
    # Low tumor purity could result in lower AF values in general. 
    # For these samples the standard threshold might be too high.
    if (random_id == "P0105" & tumor_time == "R2") {
      AF_threshold <- 0.13
    } else if (random_id == "P0121") {
      custom_AFs <- c(0.25, 0.25, 0.16, 0.13)
      names(custom_AFs) <- c("Dx", "R1", "R2", "R3")
      AF_threshold <- as.numeric(custom_AFs[tumor_time])
    } else if (random_id == "P0140") {
      custom_AFs <- c(0.22, 0.24, 0.15)
      names(custom_AFs) <- c("Dx", "R1", "R3")
      AF_threshold <- as.numeric(custom_AFs[tumor_time])
    } else if (random_id == "P0160") {
      AF_threshold <- ifelse(tumor_time == "R2", 0.09, 0.21)
    } else if (random_id == "P0164") {
      custom_AFs <- c(0.25, 0.14, 0.2)
      names(custom_AFs) <- c("Dx", "R1", "R2")
      AF_threshold <- as.numeric(custom_AFs[tumor_time])
    } else if (random_id == "P0180") {
      custom_AFs <- c(0.25, 0.25, 0.08, 0.18)
      names(custom_AFs) <- c("Dx", "R1", "R2", "R3")
      AF_threshold <- as.numeric(custom_AFs[tumor_time])
    } else {
      AF_threshold <- 0.25
    }
    df_read_based_filtered <- sample_gonl_filtered_df[sample_gonl_filtered_df[,total_read_count_idx] >= 20 &
                                                      sample_gonl_filtered_df[,read_count_alt_idx] >= 5 &
                                                      sample_gonl_filtered_df[,AF_idx] >= AF_threshold,]
    
    # Filter indels with a length >= 50
    indel_length <- function(indel_vector){
      return(max(c(nchar(indel_vector[4]), nchar(indel_vector[5])))-1)
    }
    df_read_based_filtered$indel_length <- apply(df_read_based_filtered, 1, indel_length)
    df_read_based_filtered <- dplyr::filter(df_read_based_filtered, indel_length <= 50)

    # Remove fake indels
    #fake_indels <- c()
    #for (indel_count in 1:(nrow(df_read_based_filtered)-1)){
    #  chromosome <- df_read_based_filtered$chr[indel_count]
    #  next_chromosome <- df_read_based_filtered$chr[indel_count + 1]
    #  position <- df_read_based_filtered$start[indel_count]
    #  next_position <- df_read_based_filtered$start[indel_count + 1]
    #  if (chromosome == next_chromosome & ((next_position - position) <=20)){
    #    fake_indels <- c(fake_indels, indel_count, indel_count+1)
    #  }
    #}
    
    # Subset VCF for real indels
    #real_indels <- (1:nrow(df_read_based_filtered))[!1:nrow(df_read_based_filtered) %in% fake_indels]
    #df_read_based_filtered <- df_read_based_filtered[real_indels,]
    sample_filtered_vcf_object <- sample_vcf[rownames(sample_vcf) %in% 
                                             rownames(df_read_based_filtered), ]

    # Filter on MMQ
    #INFO_realdf <- as.data.frame(info(sample_filtered_vcf_object))
    #df_read_based_filtered <- left_join(rownames_to_column(df_read_based_filtered), 
    #                                    rownames_to_column(INFO_realdf), by = "rowname")
    #high_MMQ <- c()
    #for (MMQ_values in df_read_based_filtered$MMQ){
    #  MMQ_values_60 <- MMQ_values == 60
    #  high_MMQ <- c(high_MMQ, (MMQ_values_60[1] & MMQ_values_60[2]))
    #}
    #df_read_based_filtered <- df_read_based_filtered[high_MMQ,]
    #rownames(df_read_based_filtered) <- df_read_based_filtered$rowname
    
    # Subset VCF for filtered variants
    #sample_filtered_vcf_object <- sample_vcf[rownames(sample_vcf) %in% 
    #                                         rownames(df_read_based_filtered), ]
    
    # Create an indel mutation matrix
    sample_indel_grl <- get_indel_context(rowRanges(sample_filtered_vcf_object), 
                                          ref_genome = ref_genome)
    sample_indel_counts <- count_indel_contexts(sample_indel_grl)
    colnames(sample_indel_counts) <- tumor_info
    indel_mut_mat <- cbind(indel_mut_mat, 
                           sample_indel_counts)
  }
}
indel_mut_mat <- indel_mut_mat[, 2:ncol(indel_mut_mat)]
colSums(indel_mut_mat)

#load("MR_mut_mat_DBS_profilePerTimepoint_160524.Rda")
#load("MR_mut_mat_InDel_profilePerTimepoint_160524.Rda")



#----------------------- Mutational signature analysis ------------------------#

# Function to calculate relative contributions
calculate_relative_contribution <- function(contribution_col){
  contributions <- contribution_col[1:(length(contribution_col) - 1)]
  total <- contribution_col[length(contribution_col)]
  relative_contributions <- contributions / total
  return(relative_contributions)
}

# Function to calculate absolute contributions
calculate_absolute_contribution <- function(contribution_col){
  contributions <- contribution_col[1:(length(contribution_col) - 2)]
  total <- contribution_col[length(contribution_col) - 1]
  originaltotal <- contribution_col[length(contribution_col)]
  absolute_contributions <- ((contributions / total) * originaltotal)
  return(absolute_contributions)
}



### DBS analysis
# Conduct a strict refit on the DBS matrix
dbs_fit_res <- fit_to_signatures(dbs_mut_mat, dbs_signatures)

# Calculate contribution totals based on refit
contribution_with_totals <- rbind(dbs_fit_res$contribution, 
                                  apply(dbs_fit_res$contribution, 2, sum))
rownames(contribution_with_totals)[nrow(contribution_with_totals)] <- "Total"

# Calculate relative contributions
dbs_rel_contributions <- apply(contribution_with_totals, 2, calculate_relative_contribution)

# Calculate contribution totals based on original data
contribution_with_totals <- rbind(contribution_with_totals, 
                                  apply(dbs_mut_mat, 2, sum))
rownames(contribution_with_totals)[nrow(contribution_with_totals)] <- "OriginalTotal"

# Calculate absolute contributions
dbs_abs_contributions <- apply(contribution_with_totals, 2, calculate_absolute_contribution)

# Plot profiles
pdf("dbs_profiles.pdf", width = 16, height = 120)
plot_dbs_contexts(dbs_mut_mat) +
  xlab("Sequence context") +
  theme(strip.text.y = element_text(angle = 0, size = 10),
        strip.text.x = element_text(size = 10))
dev.off()


## Write contributions to file
sample_info <- data.frame(Patient_ID = sapply(strsplit(colnames(dbs_abs_contributions), "_"), `[`, 1),
                          Timepoint = sapply(strsplit(colnames(dbs_abs_contributions), "_"), `[`, 2))

# Relative contributions
# Change contributions from NA to 0
dbs_rel_contributions[is.na(dbs_rel_contributions)] <- 0

# Combine contributions with the mutation matrix
dbs_compl_info_rel <- as.data.frame(t(rbind(dbs_mut_mat, dbs_rel_contributions)))

# Add sample info and total counts to data frame
dbs_compl_info_rel <- cbind(sample_info,
                            dbs_compl_info_rel, 
                            data.frame(Total_count_MutMat = colSums(dbs_mut_mat)))
write.xlsx(dbs_compl_info_rel, file = "MR_DBS_MutMatAndRelativeContributions_16052024.xlsx")

# Absolute contributions
# Change contributions from NA to 0 and round values
dbs_abs_contributions[is.na(dbs_abs_contributions)] <- 0
dbs_abs_contributions_rounded <- round(dbs_abs_contributions, digits = 0)

# Combine contributions with the mutation matrix
dbs_compl_info_abs <- as.data.frame(t(rbind(dbs_mut_mat, dbs_abs_contributions_rounded)))

# Add sample info and total counts to data frame
dbs_compl_info_abs <- cbind(dbs_compl_info_abs, 
                            data.frame(Total_count_MutMat = colSums(dbs_mut_mat),
                                       Total_count_AbsContr = colSums(dbs_abs_contributions_rounded)))
# Write to Excel
write.xlsx(dbs_compl_info, file = "MR_DBS_MutMatAndAbsoluteContributions_16052024.xlsx")



### InDel analysis
# Conduct a strict refit on the DBS matrix
indel_fit_res <- fit_to_signatures(indel_mut_mat, indel_signatures)

# Calculate contribution totals based on refit
contribution_with_totals <- rbind(indel_fit_res$contribution, 
                                  apply(indel_fit_res$contribution, 2, sum))
rownames(contribution_with_totals)[nrow(contribution_with_totals)] <- "Total"

# Calculate relative contributions
indel_rel_contributions <- apply(contribution_with_totals, 2, calculate_relative_contribution)

# Calculate contribution totals based on original data
contribution_with_totals <- rbind(contribution_with_totals, 
                                  apply(dbs_mut_mat, 2, sum))
rownames(contribution_with_totals)[nrow(contribution_with_totals)] <- "OriginalTotal"

# Calculate absolute contributions
indel_abs_contributions <- apply(contribution_with_totals, 2, calculate_absolute_contribution)

# Plot profiles
pdf("indel_profiles.pdf", width = 16, height = 120)
plot_indel_contexts(indel_mut_mat) +
  xlab("Sequence context") +
  theme(strip.text.y = element_text(angle = 0, size = 10),
        strip.text.x = element_text(size = 10))
dev.off()


## Write contributions to file
sample_info <- data.frame(Patient_ID = sapply(strsplit(colnames(indel_abs_contributions), "_"), `[`, 1),
                          Timepoint = sapply(strsplit(colnames(indel_abs_contributions), "_"), `[`, 2))

# Relative contributions
# Combine contributions with the mutation matrix
indel_compl_info_rel <- as.data.frame(t(rbind(indel_mut_mat, indel_rel_contributions)))

# Add sample info and total counts to data frame
indel_compl_info_rel <- cbind(sample_info,
                              indel_compl_info_rel, 
                              data.frame(Total_count_MutMat = colSums(indel_mut_mat)))
write.xlsx(indel_compl_info_rel, file = "MR_InDel_MutMatAndRelativeContributions_16052024.xlsx")

# Absolute contributions
# Round values
indel_abs_contributions_rounded <- round(indel_abs_contributions, digits = 0)

# Combine absolute contributions with the mutation matrix
indel_compl_info_abs <- as.data.frame(t(rbind(indel_mut_mat, indel_abs_contributions_rounded)))

# Add sample info and total counts to data frame
indel_compl_info_abs <- cbind(sample_info,
                              indel_compl_info_abs, 
                              data.frame(Total_count_MutMat = colSums(indel_mut_mat),
                                         Total_count_AbsContr = colSums(indel_abs_contributions_rounded)))

# Write to Excel
write.xlsx(indel_compl_info_abs, file = "MR_InDel_MutMatAndAbsoluteContributions_16052024.xlsx")




