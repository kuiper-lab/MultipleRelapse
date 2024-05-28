#!/usr/bin/env Rscript 

# Clear data and set seed
rm(list = ls())
set.seed(1234)

# Load packages
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)
library(MutationalPatterns)
library(MatrixGenerics)
library(logging)
library(tidyverse)
library(grid)
library(gridExtra)
library(stringi)
library(gtools)
library(ggplot2)
library(cowplot)
library(readxl)
library(Hmisc)
library(RColorBrewer)
library(lsa)
library(plyr)
library(bedr)
library(ensemblVEP)
library(TVTB)

# Set global parameters
# All input files should be located somewhere in PROJECT_DIR
PROJECT_DIR <- "/Users/m.m.kleisman/Projects/HypermutatedALL_project/"
RESULTS_DIR <- paste0(PROJECT_DIR, "RESULTS/multipleRelapse/")
SKION_INFO <- "skion_gender.csv"
RANDOM_ID_INFO <- "MR_cohortDescription.xlsx"
CENTROMERES_FILE <- "cytoBand.hg38.centromeresOnly.txt"
MUTATION_MATRIX <- "combined_mut_mat_MRsamples_refcohort_05032024.rdata"
NMF_RES_FILE <- "nmf_res_MRsamples_refcohort_17032024_8signatures.rdata"
SBSBLOOD_SIGNATURE <- "41586_2022_5072_MOESM4_ESM.xlsx"
SNV_EXT <- "snp-WGS.vepAnnotated.vcf.gz"
CNV_EXT <- ".merged.tumor.cnv.txt"

setwd(RESULTS_DIR)

# Load custom functions
FUNCTIONS_SCRIPT <- list.files(path = "/Users/m.m.kleisman/Projects/git/pmc_kuiper_projects/", 
                               pattern = "createdFunctions_HM-ALL.R",
                               full.names = TRUE, recursive = TRUE)
source(FUNCTIONS_SCRIPT)

# Define filter parameters
REMOVE_GERMLINE_VARIANTS <- FALSE
PLOT_AF_AND_PROFILE_PLOT <- FALSE
split <- TRUE
countClusterSize <- FALSE
deleteSmallClusters <- TRUE
MIN_AF <- 0.25
MinClustSize <- 75
noise_threshold <- 0      

# Obtain paths to input files 
centromeres <- list.files(path = PROJECT_DIR,
                          pattern = CENTROMERES_FILE, 
                          full.names = TRUE, 
                          recursive = TRUE)
vcf_vector <- list.files(pattern = SNV_EXT, 
                         full.names = TRUE, 
                         recursive = TRUE)
cnv_files <- list.files(pattern = CNV_EXT, 
                        full.names = TRUE, 
                        recursive = TRUE)

# Read SKION and random ID info of the MR relapse cohort
skion_gender_df <- as.data.frame(read.csv(list.files(path = PROJECT_DIR,
                                                     pattern = SKION_INFO, 
                                                     full.names = TRUE, 
                                                     recursive = TRUE))) 
random_id_info <- as.data.frame(read_excel(list.files(path = PROJECT_DIR,
                                                      pattern = RANDOM_ID_INFO,
                                                      recursive = TRUE,
                                                      full.names = TRUE))[1:29, c(1,2)])
colnames(random_id_info) <- c("skion_nr", "random_id")

# Load and parse random ID info of the reference cohort
ref_random_id_info <- as.data.frame(read_excel(list.files(path = PROJECT_DIR,
                                                          pattern = "refCohort_sampleDescription_25032024.xlsx",
                                                          recursive = TRUE,
                                                          full.names = TRUE)))
rownames(ref_random_id_info) <- ref_random_id_info$SKION
ref_random_id_info <- ref_random_id_info[,-1]

# Combine random ID info
compl_random_id_info <- rbind(random_id_info,
                              data.frame(skion_nr = rownames(ref_random_id_info),
                                         random_id = ref_random_id_info$`Patient ID`))

# Load mutation matrix, NMF result and SBSblood
load(list.files(pattern = MUTATION_MATRIX, 
                full.names = TRUE, recursive = TRUE))
load(list.files(pattern = NMF_RES_FILE, 
                full.names = TRUE, recursive = TRUE))
SBSblood <- as.data.frame(read_excel(list.files(path = PROJECT_DIR, 
                                                pattern = SBSBLOOD_SIGNATURE, 
                                                full.names = TRUE, recursive = TRUE),
                                     sheet = 8, skip = 3))


# Define directory where to store output files
output_dir <- paste0(RESULTS_DIR, "ClusterAFDeNovoContributionPlots")
if (!dir.exists(output_dir)){dir.create(output_dir)}



#---------------------------- Function definitions ---------------------------#

# Define cluster order
tumor_types <- c("ClDx", do.call(paste0, expand.grid("ClR", 1:8)))
cluster_types <- as.character(sapply(tumor_types, function(x) 
  paste0(x, c("", "_A", "_B", "_sub", "_sub_A", "_sub_B", "_sub_C", "_fall", "_fall_A", 
              "_fall_B", "_fall_sub", "_fall_sub_A", "_fall_sub_B"))))
na_clusters <- paste0("NA_", c(paste0(0, 1:9), 10:11))
cluster_order <- c(cluster_types, na_clusters)

# Functions
splitVcfBasedOnCluster <- function (vcfObject, clusterDataframe){
  
  #Adapted so "cluster" does not get added to the cluster name
  cluster_list <- lapply(unique(clusterDataframe$cluster), 
                         function(index) {
                           master_list <- list()
                           cluster_specific_index <- rownames(clusterDataframe)[clusterDataframe$cluster == 
                                                                                  index]
                           cluster_specific_variants <- vcfObject[rownames(vcfObject) %in% 
                                                                    cluster_specific_index]
                           master_list[[index]] <- rowRanges(cluster_specific_variants)
                         })
  g_ranges_cluster_list <- GRangesList(cluster_list)
  header <- unique(clusterDataframe$cluster)
  names(g_ranges_cluster_list) <- header
  return(g_ranges_cluster_list)
}

ggplot_cnv_af  <- function (melted_df) {
  sample_af_merged_plots <- ggplot(data = melted_df, aes(x = time_point,
                                                         y = AF,
                                                         group = variant,
                                                         colour=as.factor(Type))) +
    scale_color_manual(name = "CNA Type", values = c("GAIN" = "#56B4E9", "LOSS" = "#CC79A7", 
                                                     "NO CNA" = "#000000", "CHR X" = "#009E73")) +
    geom_line(alpha = 0.02) +
    geom_point(alpha = 0.30) +
    facet_wrap(~cluster, ncol = 1) +
    scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.5)) +
    coord_cartesian(xlim = c(x_width, 1.45), ylim = c(0,1)) +
    theme_bw() +
    xlab(label = "Time point") +
    ylab(label = "Allele frequency") + 
    theme(axis.text = element_text(size = 10, color = "black"),
          axis.title = element_text(size = 11.5, color = "black"),
          strip.text = element_text(size = 10))
  
  return(sample_af_merged_plots)
}

renameCluster <- function (cluster_pattern, timepoints){
  #Name the clusters based on the cluster pattern eg 0-1-1
  if (length(cluster_pattern) == sum(cluster_pattern)) {
    return("ClDx")
  }
  first_occ_present <- match(1, cluster_pattern)
  
  #nonsense cluster
  if (is.na(first_occ_present)) {
    return("NA")
  }
  
  t <- str_replace(timepoints, skion_number, "")
  if (cluster_pattern[length(cluster_pattern)] == 0) {
    last_occ_present <- which(cluster_pattern == 1)
    first_occ_notpresent <- match(0, cluster_pattern)
    leftover_vec <- cluster_pattern[first_occ_present:last_occ_present[length(last_occ_present)]]
    exp_val_fall <- length(leftover_vec)
    if (exp_val_fall == sum(leftover_vec)) {
      if (first_occ_present == 1) {
        cluster_name <- paste0("Cl", "Dx", "_fall")
      } else {
        #cluster_name <- paste0("Cl", "R", (first_occ_present-1), "_fall")
        cluster_name <- paste0("Cl", t[first_occ_present], "_fall")
      }
      return(cluster_name)
    }
  }
  
  if (first_occ_present > 1) {
    leftover_vec <- cluster_pattern[first_occ_present:length(cluster_pattern)]
    exp_val_ris <- length(leftover_vec) 
    leftover_sum <- sum(leftover_vec)
    if (exp_val_ris == leftover_sum) {
      #cluster_name <- paste0("Cl", "R", (first_occ_present-1))
      cluster_name <- paste0("Cl", t[first_occ_present])
      return(cluster_name)
    }
  } 
  return("NA")
}

splitClusterOnTimePointSpecificAlleleFrequency <- function (cluster_df, time_point, cluster, AF_threshold = 0){ 
  #Split clusters, and name them
  if (!(cluster %in% cluster_df$cluster)) {
    print(paste0(cluster, " was not found amongst the existing clusters, check spelling?"))
    stop()
  }
  loginfo(paste0("splitting cluster: ", cluster, " based on timepoint: ", time_point))
  time_point <- grep(time_point, colnames(cluster_df), value = TRUE) 
  
  #check if the cluster has already been split before, if so give a unique postfix, otherwise the clusters will get merged
  if (grepl("_sub", cluster, fixed = TRUE) && cluster %in% cluster_df$cluster) {
    cluster_df$cluster <- gsub(cluster, paste0(cluster, "_A"), cluster_df$cluster)
    cluster_df$cluster[cluster_df$cluster == cluster] <- paste0(cluster, "_A")
    cluster_a_idx <- filterVariantsOnTimePointSpecificAlleleFrequency(cluster_df = cluster_df, 
                                                                      time_point = time_point, cluster = paste0(cluster, "_A"), 
                                                                      AF_threshold = AF_threshold, variants_to_keep = FALSE)
    cluster_df$cluster[rownames(cluster_df) %in% rownames(cluster_a_idx)] <- paste0(cluster, "_B")
    
  } else if (paste0(cluster, "_sub") %in% cluster_df$cluster) {
    cluster_df$cluster <- gsub(paste0(cluster, "_sub"), paste0(cluster, "_sub_A"), cluster_df$cluster)
    cluster_df$cluster[cluster_df$cluster == cluster] <- paste0(cluster, "_sub_B")
    cluster_a_idx <- filterVariantsOnTimePointSpecificAlleleFrequency(cluster_df = cluster_df, 
                                                                      time_point = time_point, cluster = paste0(cluster, "_sub_B"), 
                                                                      AF_threshold = AF_threshold, variants_to_keep = FALSE)
    cluster_df$cluster[rownames(cluster_df) %in% rownames(cluster_a_idx)] <- paste0(cluster)
  }
  
  else {
    cluster_df$cluster[cluster_df$cluster == cluster] <- paste0(cluster, "_sub")
    cluster_a_idx <- filterVariantsOnTimePointSpecificAlleleFrequency(cluster_df = cluster_df, 
                                                                      time_point = time_point, cluster = paste0(cluster, "_sub"), 
                                                                      AF_threshold = AF_threshold, variants_to_keep = FALSE)
    cluster_df$cluster[rownames(cluster_df) %in% rownames(cluster_a_idx)] <- paste0(cluster)
  }
  unique(cluster_df$cluster)
  return(cluster_df)
}

plotClusterSpecificAlleleFrequency <- function (alleleFrequencyDataFrame, minimal = FALSE){
  #Adapted so that "cluster" does not get added to the cluster name
  
  create_k_means_allele_frequency_plot_df <- function(alleleFrequencyDataFrame, 
                                                      sample_list) {
    n_samples <- length(sample_list)
    var = rep(x = row.names(alleleFrequencyDataFrame), times = n_samples)
    AF <- c()
    for (i in seq(sample_list)) {
      AF <- c(AF, as.vector(unlist(alleleFrequencyDataFrame[[i]])))
    }
    sample <- c()
    for (i in seq(sample_list)) {
      sample <- c(sample, rep(x = names(alleleFrequencyDataFrame)[i], 
                              times = nrow(alleleFrequencyDataFrame)))
    }
    cluster <- c(rep(x = alleleFrequencyDataFrame$cluster, 
                     times = n_samples))
    return(data.frame(var = var, AF = AF, sample = sample, 
                      cluster = cluster))
  }
  names(alleleFrequencyDataFrame) <- gsub("_read_based_AF", 
                                          "", names(alleleFrequencyDataFrame), perl = TRUE)
  ordered_vcf_sample_ids <- orderSampleIds(names(alleleFrequencyDataFrame))
  k_means_clustering_plot_df <- create_k_means_allele_frequency_plot_df(alleleFrequencyDataFrame, 
                                                                        ordered_vcf_sample_ids)
  k_means_clustering_plot_df$sample <- factor(k_means_clustering_plot_df$sample, 
                                              levels = ordered_vcf_sample_ids)
  tumor_k_means_clustering_plot_df <- k_means_clustering_plot_df
  if (minimal) {
    axis_text_x <- element_blank()
    axis_text_y <- element_blank()
    axis_title_x <- element_blank()
    axis_title_y <- element_blank()
  }
  else {
    axis_text_x <- element_text(angle = 45, hjust = 1)
    axis_text_y <- element_text(angle = 90)
    axis_title_x <- element_text()
    axis_title_y <- element_text()
  }
  cluster_allel_frequency_plots <- lapply(unique(alleleFrequencyDataFrame$cluster), 
                                          function(x) {
                                            cluster_specific_data_slice <- tumor_k_means_clustering_plot_df[tumor_k_means_clustering_plot_df$cluster == x, ]
                                            ggplot(data = cluster_specific_data_slice, aes(x = sample, 
                                                                                           y = AF, group = var)) + geom_point(alpha = 0.2) + 
                                              geom_line(alpha = 0.2) + facet_wrap(~cluster) + 
                                              facet_grid(cols = vars(cluster)) + ylim(0, 1) + 
                                              theme(panel.background = element_rect(fill = "white", colour = "grey50"), 
                                                    text = element_text(size = 16), 
                                                    legend.position = "none", panel.grid.major.y = element_line(colour = "black", size = 0.5), 
                                                    axis.text.x = axis_text_x, axis.title.x = axis_title_x, 
                                                    axis.text.y = axis_text_y, axis.title.y = axis_title_y)})
  names(cluster_allel_frequency_plots) <- unique(alleleFrequencyDataFrame$cluster)
  return(cluster_allel_frequency_plots)
}

plot_contribution_custom <- function (contribution, signatures, index = c(), coord_flip = FALSE, 
                                      mode = "relative", palette = c()){
  #customized plots for MR study
  if (!(mode == "relative" | mode == "absolute")) {
    stop("mode parameter should be either 'relative' or 'absolute'")
  }
  
  if (length(index > 0)) {
    contribution = contribution[, index]
  }
  Sample = NULL
  Contribution = NULL
  Signature = NULL
  
  if (mode == "relative") {
    m_contribution = melt(contribution)
    colnames(m_contribution) = c("Signature", "Sample", "Contribution")
    m_contribution$Sample <- factor(m_contribution$Sample, levels = c(colnames(contribution)))
    plot = ggplot(m_contribution, aes(x = factor(Sample), y = Contribution, 
                                      fill = factor(Signature), order = Sample)) + 
      geom_bar(position = "stack", stat = "identity", colour = "NA") + # here I changed colour into NA
      labs(x = "", y = "Relative contribution") + theme_bw() + 
      theme(panel.grid.minor.x = element_blank(), 
            panel.grid.major.x = element_blank(),
            panel.grid.minor.y = element_blank(), 
            panel.grid.major.y = element_blank(),
            axis.text.x = element_text(size = 10, color = "black"),
            axis.text.y = element_text(size = 10, color = "black"),
            axis.title = element_text(size = 11.5, color = "black"),
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 10))
  }
  else {
    if (missing(signatures)) {
      stop(paste("For contribution plotting in mode 'absolute':", 
                 "also provide signatures matrix"))
    }
    
    total_signatures = colSums(signatures)
    abs_contribution = contribution 
    m_contribution = reshape2::melt(abs_contribution)
    colnames(m_contribution) = c("Signature", "Sample", "Contribution")
    m_contribution$Sample <- factor(m_contribution$Sample, levels = c(colnames(contribution)))
    
    m_contribution[,"IsUnassigned"] <- NA
    for (i in 1:nrow(m_contribution)) {
      m_contribution[i, "IsUnassigned"] <- m_contribution[m_contribution$Sample == m_contribution$Sample[i] &
                                                            m_contribution$Signature == "Unassigned","Contribution"]
      if (m_contribution$IsUnassigned[i] == 0){
        m_contribution$IsUnassigned[i]<-NA
      }
    }
    m_contribution$Signature <- factor(m_contribution$Signature, levels = rev(unique(m_contribution$Signature)))
    plot = ggplot(m_contribution, aes(x = factor(Sample), y = Contribution, 
                                      fill = factor(Signature))) + 
      geom_bar(stat = "identity", position = "stack", linewidth = 0.15) + 
      labs(x = "", y = "Absolute contribution (no. mutations)") + 
      theme_bw() + theme(panel.grid.minor.x = element_blank(), 
                         panel.grid.major.x = element_blank(),
                         panel.grid.minor.y = element_blank(), 
                         panel.grid.major.y = element_blank(),
                         axis.text.x = element_text(size = 10, color = "black"),
                         axis.text.y = element_text(size = 10, color = "black"),
                         axis.title = element_text(size = 11.5, color = "black"),
                         legend.title = element_text(size = 12),
                         legend.text = element_text(size = 10)) 
  }
  if (length(palette) > 0) {
    plot = plot + scale_fill_manual(name = "Signature", values = palette, guide = guide_legend(reverse = TRUE))
  } else {
    plot = plot + scale_fill_discrete(name = "Signature", guide = guide_legend(reverse = TRUE))
  }
  
  if (coord_flip) {
    plot = plot + coord_flip() + xlim(rev(levels(factor(m_contribution$Sample)))) 
  } else {
    plot = plot + xlim(levels(factor(m_contribution$Sample)))
  }
  return(plot)
}



#--------------- Data processing and mutation matrix generation ---------------#

## Filter the vcf files and extract VAFs for each timepoint in each patient
filtered_vcf_objects_list <- rep(list(NA), length = length(vcf_vector))
filtered_variants_clusters_df_list <- rep(list(NA), length = length(vcf_vector))
for (i in 1:length(vcf_vector)){
  vcf_path <- vcf_vector[i]
  skion_number <- str_extract(vcf_path, "\\d\\d\\d\\d+")
  random_id <- random_id_info[random_id_info$skion_nr == skion_number, 2]
    
  # Read VCF and rename clusters
  sample_vcf <- loadVcfFile(vcf_path,
                            ref_genome,
                            "chromosomes")
  GenomeInfoDb::genome(sample_vcf)  <- "hg38"
  
  for (sample_id in colnames(sample_vcf)){
    if (endsWith(sample_id, skion_number)){
      time <- str_split(sample_id, "-")[[1]][1]
      new_name <- paste0(skion_number, time)
    } else if (startsWith(sample_id, "P_")){
      time <- str_split(sample_id, "_")[[1]][3]
      new_name <- paste0(skion_number, time)
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
  
  # Remove germline variants if requested
  if (REMOVE_GERMLINE_VARIANTS){
    sample_vcf_no_reads_in_normal_fraction <- filterOnReadsInNormalSamples(vcfObject = sample_vcf_no_centromeric_variants, 
                                                                           allowedNumberOfReads = 0)
  } else {
    sample_vcf_no_reads_in_normal_fraction <- sample_vcf_no_centromeric_variants
  }
  
  # Apply popMax filters
  sample_vcf_gnomad_filtered <- filterOnGnomadAlleleFrequency(sample_vcf_no_reads_in_normal_fraction,
                                                              0.01)
  sample_vcf_gonl_filtered <- filterOnGonlAlleleFrequency(sample_vcf_gnomad_filtered,
                                                          0.01)
  sample_gonl_filtered_df <- convertVcfObjectToDataFrame(sample_vcf_gonl_filtered)
  
  # Make data frame with only tumor data information
  timepoints <- grep("D|([^C]R)", orderSampleIds(colnames(sample_vcf_gonl_filtered)), value = T)
  for (timepoint in 1:length(timepoints)){
    time_name <- timepoints[timepoint]
    time_df <- sample_gonl_filtered_df[,grep(time_name, colnames(sample_gonl_filtered_df))]
    time_df <- time_df[,!grepl("mutect", colnames(time_df))]
    time_df[, paste0(time_name, "_read_count_total")] <- time_df[, paste0(time_name, "_read_count_ref")] +
      time_df[, paste0(time_name, "_read_count_alt")]
    if (timepoint == 1){
      tumor_df <- time_df
    } else{
      tumor_df <- cbind(tumor_df, time_df)
    }
  }
  
  # Define read info columns to which filters need to be applied
  alt_cols <- grep("read_count_alt", colnames(tumor_df))
  total_cols <- grep("read_count_total", colnames(tumor_df))
  vaf_cols <- grep("read_based_AF", colnames(tumor_df))
  
  # Apply filters. Use custom AF filters per sample when necessary. 
  # Low tumor purity could result in lower AF values in general. 
  # For these samples the standard threshold might be too high.
  if (random_id == "P0105") {
    sample_filtered_df <- sample_gonl_filtered_df[(apply(tumor_df[,alt_cols], 1, max) >= 3) &
                                                    (apply(tumor_df[,total_cols], 1, max) >= 20) &
                                                    (tumor_df[,3] >= 0.25 | tumor_df[,7] >= 0.25 | tumor_df[,11] >= 0.13),]
  } else if (random_id == "P0121") {
    sample_filtered_df <- sample_gonl_filtered_df[(apply(tumor_df[,alt_cols], 1, max) >= 3) &
                                                    (apply(tumor_df[,total_cols], 1, max) >= 20) &
                                                    (tumor_df[,3] >= 0.25 | tumor_df[,7] >= 0.25 | tumor_df[,11] >= 0.16 | tumor_df[,15] >= 0.13),]
  } else if (random_id == "P0140") {
    sample_filtered_df <- sample_gonl_filtered_df[(apply(tumor_df[,alt_cols], 1, max) >= 3) &
                                                    (apply(tumor_df[,total_cols], 1, max) >= 20) &
                                                    (tumor_df[,3] >= 0.22 | tumor_df[,7] >= 0.24 | tumor_df[,11] >= 0.15),]
  } else if (random_id == "P0160") {
    sample_filtered_df <- sample_gonl_filtered_df[(apply(tumor_df[,alt_cols], 1, max) >= 3) &
                                                    (apply(tumor_df[,total_cols], 1, max) >= 20) &
                                                    (tumor_df[,3] >= 0.21 | tumor_df[,7] >= 0.21 | tumor_df[,11] >= 0.09),]
  } else if (random_id == "P0164") {
    sample_filtered_df <- sample_gonl_filtered_df[(apply(tumor_df[,alt_cols], 1, max) >= 3) &
                                                    (apply(tumor_df[,total_cols], 1, max) >= 20) &
                                                    (tumor_df[,3] >= 0.25 | tumor_df[,7] >= 0.14 | tumor_df[,11] >= 0.2),]
  } else if (random_id == "P0180") {
    sample_filtered_df <- sample_gonl_filtered_df[(apply(tumor_df[,alt_cols], 1, max) >= 3) &
                                                    (apply(tumor_df[,total_cols], 1, max) >= 20) &
                                                    (tumor_df[,3] >= 0.25 | tumor_df[,7] >= 0.25 | tumor_df[,11] >= 0.08 | tumor_df[,15] >= 0.18),]
  } else {
    sample_filtered_df <- sample_gonl_filtered_df[(apply(tumor_df[,alt_cols], 1, max) >= 3) &
                                                    (apply(tumor_df[,total_cols], 1, max) >= 20) &
                                                    (apply(tumor_df[,vaf_cols], 1, max) >= MIN_AF),]
  }
  sample_hg38_filtered <- sample_vcf_gonl_filtered[rownames(sample_filtered_df)]
  sample_af_matrix <- extractAllelFrequenciesFromVcf(sample_hg38_filtered)
  sample_af_df <- sample_af_matrix[, (seq(ncol(sample_af_matrix)/4) * 4) - 1]
  loginfo(paste0("Using the read based allele frequencies of the vcf-file to cluster variants"))
  
  sample_tumour_ids <- grep("D|([^C]R)", colnames(sample_af_df))
  sample_af_df_tumours <- sample_af_df[,sample_tumour_ids]
  sample_transformed_tumours <- sample_af_df_tumours
  colnames(sample_transformed_tumours) <- orderSampleIds(colnames(sample_transformed_tumours))
  
  for (colname in colnames(sample_transformed_tumours)) {
    sample_transformed_tumours[,colname] = ifelse(sample_af_df_tumours[,colname] > noise_threshold, 1, 0)
  }
  
  sample_countdf <- plyr::count(sample_transformed_tumours)
  sample_significant_patterns <- sample_countdf
  new_rownames <- c()
  NA_counter <- 0
  
  for (pattern_row in 1:nrow(sample_significant_patterns)) {
    pattern <- sample_significant_patterns[pattern_row,1:(ncol(sample_significant_patterns)-1)]
    cluster_name <- renameCluster(pattern,timepoints)
    
    #check if a falling cluster already exists, if so give them a postfix in order to make the name unique
    if (cluster_name %in% new_rownames) {
      new_rownames <- gsub(cluster_name, paste0(cluster_name, "_A"), new_rownames)
      cluster_name <- paste0(cluster_name, "_B")
    } else if (paste0(cluster_name, "_G") %in% new_rownames) {
      cluster_name <- paste0(cluster_name, "_H")
    } else if (paste0(cluster_name, "_F") %in% new_rownames) {
      cluster_name <- paste0(cluster_name, "_G")
    } else if (paste0(cluster_name, "_E") %in% new_rownames) {
      cluster_name <- paste0(cluster_name, "_F")
    } else if (paste0(cluster_name, "_D") %in% new_rownames) {
      cluster_name <- paste0(cluster_name, "_E")
    } else if (paste0(cluster_name, "_C") %in% new_rownames) {
      cluster_name <- paste0(cluster_name, "_D")
    } else if (paste0(cluster_name, "_B") %in% new_rownames) {
      cluster_name <- paste0(cluster_name, "_C")
    }
    
    if (cluster_name == "NA") {
      NA_counter = NA_counter + 1
      if (NA_counter < 10) {
        cluster_name <- paste0("NA_", "0", NA_counter)
      } else {
        cluster_name <- paste0("NA_", NA_counter)
      }
    } 
    new_rownames <- append(new_rownames, cluster_name)
  }
  rownames(sample_significant_patterns) <- new_rownames
  
  # Assign clustername to each mutation
  clusters <- c()
  for (AF_row in 1:nrow(sample_transformed_tumours)) {
    for (pattern_row in 1:nrow(sample_significant_patterns)) {
      if (identical(as.numeric(sample_transformed_tumours[AF_row, ]), as.numeric(sample_significant_patterns[pattern_row, 1:ncol(sample_significant_patterns)-1]))) {
        clusters <- append(clusters, rownames(sample_significant_patterns[pattern_row,]))
        break
      }
    }
  }
  sample_tumours_clustered <- as.data.frame(sample_af_df_tumours)
  sample_tumours_clustered$cluster <- clusters
  sample_filtered_variants_clusters_df <- sample_tumours_clustered
  
  # Save variant info to objects
  filtered_vcf_objects_list[[i]] <- sample_hg38_filtered
  names(filtered_vcf_objects_list)[i] <- skion_number
  filtered_variants_clusters_df_list[[i]] <- sample_filtered_variants_clusters_df
  names(filtered_variants_clusters_df_list)[i] <- skion_number
}



## Create histograms of VAFs of each cluster at each timepoint
pdf(paste0(output_dir, "/clusterVAFdensities_smallClustersDeleted.pdf"), width = 27, height = 25)
for (i in 1:length(filtered_variants_clusters_df_list)){ 
  skion_number <- names(filtered_variants_clusters_df_list)[i]
  random_id <- random_id_info[random_id_info$skion_nr == skion_number, 2]
  
  # Parse and order column names for readability
  sample_variant_clusters <- filtered_variants_clusters_df_list[[i]]
  colnames(sample_variant_clusters) <- gsub(skion_number, "", colnames(sample_variant_clusters))
  colnames(sample_variant_clusters) <- gsub("_read_based_AF", "", colnames(sample_variant_clusters))
  
  # Order data frame by tumor
  sample_variant_clusters <- sample_variant_clusters[, c(mixedsort(colnames(sample_variant_clusters)[1:(ncol(sample_variant_clusters) - 1)]), "cluster")]
  
  # Delete small clusters
  small_clusters <- names(which(table(sample_variant_clusters$cluster) < MinClustSize))
  sample_variant_clusters <- sample_variant_clusters[!sample_variant_clusters$cluster %in% small_clusters, ]
  
  ## Create function of required plotting commands
  generate_vaf_overviews <- function(variant_clusters_df,
                                     plot_title = random_id){
    # Melt data frame
    sample_variant_clusters_melt <- reshape2::melt(variant_clusters_df)
    sample_variant_clusters_melt["grouping"] <- rep(1:nrow(variant_clusters_df), 
                                                    n = (ncol(variant_clusters_df)-1))
    
    # Define general plot parameters
    plot_params <- theme(axis.text.x = element_text(angle = -90, size = 12),
                         axis.text.y = element_text(size = 12),
                         strip.text = element_text(size = 20))
    
    # Generate VAF overview of each cluster
    cluster_vafs <- ggplot(sample_variant_clusters_melt, 
                           aes(x = variable, y = value, group = grouping)) +
      geom_line(alpha = 0.3) +
      geom_point() +
      facet_wrap(. ~ cluster, ncol = 1) +
      scale_y_continuous(limits = c(0, 1)) +
      xlab("Time point") + ylab("VAF") + theme_bw() + 
      plot_params
    
    # Generate histogram of VAF frequency at each time point
    tp_vafs <- ggplot(sample_variant_clusters_melt, aes(x = value)) +
      geom_histogram(binwidth = 0.01) +
      scale_x_continuous(breaks = seq(0, 1, 0.05), limits = c(-0.01, 1.01), 
                         expand = c(0.001, 0.001)) +
      scale_y_continuous(expand = c(0.05, 0)) +
      facet_wrap(cluster ~ variable, scales = "free", 
                 ncol = length(unique(sample_variant_clusters_melt$variable))) +
      xlab("VAF") + theme_bw() + 
      plot_params
    
    # Combine plots
    grid.arrange(cluster_vafs, tp_vafs, 
                 ncol = 2, widths = c(1, 4),
                 top = textGrob(plot_title, gp = gpar(fontsize = 25, font = 2)))
  }
  
  ## Check if sample contains > 5 NA clusters
  #  If so, divide figures over 2 pages
  if (sum(grepl("NA", unique(sample_variant_clusters$cluster))) > 5){
    main_sample_variant_clusters <- sample_variant_clusters[!grepl("NA", sample_variant_clusters$cluster),]
    na_sample_variant_clusters <- sample_variant_clusters[grepl("NA", sample_variant_clusters$cluster),]
    generate_vaf_overviews(variant_clusters_df = main_sample_variant_clusters,
                           plot_title = paste0(random_id, "_mainClusters"))
    generate_vaf_overviews(variant_clusters_df = na_sample_variant_clusters,
                           plot_title = paste0(random_id, "_NAclusters"))
    
  # Else, just call function to generate VAF overview plots
  } else {
    generate_vaf_overviews(variant_clusters_df = sample_variant_clusters)
  }
}
dev.off()

#load(paste0(output_dir, "/filtered_variants_clusters_df_list.Rda"))
#load(paste0(output_dir, "/filtered_vcf_objects_list.Rda"))

## Create a mutation matrix for each cluster
filtered_variants_clusters_df_list_split <- rep(list(NA), length = length(vcf_vector))
MR_mut_mat <- matrix(nrow = 96)
for (i in 1:length(filtered_variants_clusters_df_list)){ 
  skion_number <- names(filtered_variants_clusters_df_list)[i]
  random_id <- random_id_info[random_id_info$skion_nr == skion_number, 2]
  sample_filtered_vcf <- filtered_vcf_objects_list[[i]]
  sample_filtered_variants_clusters_df <- filtered_variants_clusters_df_list[[i]]

  # Split clusters where necessary
  if (split == TRUE) {
    if (random_id == "P0164") {
      sample_filtered_variants_clusters_df <- splitClusterOnTimePointSpecificAlleleFrequency(cluster_df = sample_filtered_variants_clusters_df,
                                                                                             time_point = "Dx",
                                                                                             cluster = "ClDx",
                                                                                             AF_threshold = 0.13)
      sample_filtered_variants_clusters_df <- splitClusterOnTimePointSpecificAlleleFrequency(cluster_df = sample_filtered_variants_clusters_df,
                                                                                             time_point = "R1",
                                                                                             cluster = "ClR1",
                                                                                             AF_threshold = 0.05)
    } 
    if (random_id == "P0608") {
      sample_filtered_variants_clusters_df <- splitClusterOnTimePointSpecificAlleleFrequency(cluster_df = sample_filtered_variants_clusters_df,
                                                                                             time_point = "Dx",
                                                                                             cluster = "ClDx",
                                                                                             AF_threshold = 0.15)
    } 
    if (random_id == "P0614") {
      sample_filtered_variants_clusters_df <- splitClusterOnTimePointSpecificAlleleFrequency(cluster_df = sample_filtered_variants_clusters_df,
                                                                                             time_point = "Dx",
                                                                                             cluster = "ClDx",
                                                                                             AF_threshold = 0.11) 
      sample_filtered_variants_clusters_df <- splitClusterOnTimePointSpecificAlleleFrequency(cluster_df = sample_filtered_variants_clusters_df,
                                                                                             time_point = "R1",
                                                                                             cluster = "ClDx_sub",
                                                                                             AF_threshold = 0.18)
      sample_filtered_variants_clusters_df <- splitClusterOnTimePointSpecificAlleleFrequency(cluster_df = sample_filtered_variants_clusters_df,
                                                                                             time_point = "R1",
                                                                                             cluster = "ClR1",
                                                                                             AF_threshold = 0.21)
      sample_filtered_variants_clusters_df <- splitClusterOnTimePointSpecificAlleleFrequency(cluster_df = sample_filtered_variants_clusters_df,
                                                                                             time_point = "R8",
                                                                                             cluster = "ClR1",
                                                                                             AF_threshold = 0.05)
    }
    if (random_id == "P0180") {
      sample_filtered_variants_clusters_df <- splitClusterOnTimePointSpecificAlleleFrequency(cluster_df = sample_filtered_variants_clusters_df,
                                                                                             time_point = "Dx",
                                                                                             cluster = "ClDx",
                                                                                             AF_threshold = 0.25)
    } 
    if (random_id == "P0616") {
      sample_filtered_variants_clusters_df <- splitClusterOnTimePointSpecificAlleleFrequency(cluster_df = sample_filtered_variants_clusters_df,
                                                                                             time_point = "R1",
                                                                                             cluster = "ClR1",
                                                                                             AF_threshold = 0.13)
    } 
    if (random_id == "P0066") {
      sample_filtered_variants_clusters_df <- splitClusterOnTimePointSpecificAlleleFrequency(cluster_df = sample_filtered_variants_clusters_df,
                                                                                             time_point = "Dx",
                                                                                             cluster = "ClDx",
                                                                                             AF_threshold = 0.2)
    }
    if (random_id == "P0077") {
      sample_filtered_variants_clusters_df <- splitClusterOnTimePointSpecificAlleleFrequency(cluster_df = sample_filtered_variants_clusters_df,
                                                                                             time_point = "Dx",
                                                                                             cluster = "ClDx",
                                                                                             AF_threshold = 0.2)
      sample_filtered_variants_clusters_df <- splitClusterOnTimePointSpecificAlleleFrequency(cluster_df = sample_filtered_variants_clusters_df,
                                                                                             time_point = "R1",
                                                                                             cluster = "ClR1",
                                                                                             AF_threshold = 0.24)
      sample_filtered_variants_clusters_df <- splitClusterOnTimePointSpecificAlleleFrequency(cluster_df = sample_filtered_variants_clusters_df,
                                                                                             time_point = "R2",
                                                                                             cluster = "ClR1",
                                                                                             AF_threshold = 0.22)
    } 
    if (random_id == "P0122") { 
      sample_filtered_variants_clusters_df <- splitClusterOnTimePointSpecificAlleleFrequency(cluster_df = sample_filtered_variants_clusters_df,
                                                                                             time_point = "Dx",
                                                                                             cluster = "ClDx",
                                                                                             AF_threshold = 0.1) 
      sample_filtered_variants_clusters_df <- splitClusterOnTimePointSpecificAlleleFrequency(cluster_df = sample_filtered_variants_clusters_df,
                                                                                             time_point = "R1",
                                                                                             cluster = "ClR1",
                                                                                             AF_threshold = 0.15)
    } 
    if (random_id == "P0144") { 
      sample_filtered_variants_clusters_df <- splitClusterOnTimePointSpecificAlleleFrequency(cluster_df = sample_filtered_variants_clusters_df,
                                                                                             time_point = "R1",
                                                                                             cluster = "ClR1",
                                                                                             AF_threshold = 0.11)
      sample_filtered_variants_clusters_df <- splitClusterOnTimePointSpecificAlleleFrequency(cluster_df = sample_filtered_variants_clusters_df,
                                                                                             time_point = "R3",
                                                                                             cluster = "ClR3",
                                                                                             AF_threshold = 0.11)   
      sample_filtered_variants_clusters_df <- splitClusterOnTimePointSpecificAlleleFrequency(cluster_df = sample_filtered_variants_clusters_df,
                                                                                             time_point = "R4",
                                                                                             cluster = "ClR4",
                                                                                             AF_threshold = 0.1) 
    }
    if (random_id == "P0160") {
      sample_filtered_variants_clusters_df <- splitClusterOnTimePointSpecificAlleleFrequency(cluster_df = sample_filtered_variants_clusters_df,
                                                                                             time_point = "R1",
                                                                                             cluster = "ClDx",
                                                                                             AF_threshold = 0.16)
      sample_filtered_variants_clusters_df <- splitClusterOnTimePointSpecificAlleleFrequency(cluster_df = sample_filtered_variants_clusters_df,
                                                                                             time_point = "R1",
                                                                                             cluster = "ClR1",
                                                                                             AF_threshold = 0.13)
    } 
    if (random_id == "P0613") {
      sample_filtered_variants_clusters_df <- splitClusterOnTimePointSpecificAlleleFrequency(cluster_df = sample_filtered_variants_clusters_df,
                                                                                             time_point = "Dx",
                                                                                             cluster = "ClDx",
                                                                                             AF_threshold = 0.15)
      sample_filtered_variants_clusters_df <- splitClusterOnTimePointSpecificAlleleFrequency(cluster_df = sample_filtered_variants_clusters_df,
                                                                                             time_point = "R2",
                                                                                             cluster = "ClDx",
                                                                                             AF_threshold = 0.2)
    } 
    if (random_id == "P0174") { 
      sample_filtered_variants_clusters_df <- splitClusterOnTimePointSpecificAlleleFrequency(cluster_df = sample_filtered_variants_clusters_df,
                                                                                             time_point = "R1",
                                                                                             cluster = "ClR1",
                                                                                             AF_threshold = 0.22)
    } 
    if (random_id == "P0618") { 
      sample_filtered_variants_clusters_df <- splitClusterOnTimePointSpecificAlleleFrequency(cluster_df = sample_filtered_variants_clusters_df,
                                                                                             time_point = "R2",
                                                                                             cluster = "ClR2",
                                                                                             AF_threshold = 0.11)
      sample_filtered_variants_clusters_df <- splitClusterOnTimePointSpecificAlleleFrequency(cluster_df = sample_filtered_variants_clusters_df,
                                                                                             time_point = "Dx",
                                                                                             cluster = "NA_03",
                                                                                             AF_threshold = 0.11)  
    } 
    if (random_id == "P0619") { 
      sample_filtered_variants_clusters_df <- splitClusterOnTimePointSpecificAlleleFrequency(cluster_df = sample_filtered_variants_clusters_df,
                                                                                             time_point = "R1",
                                                                                             cluster = "ClR1",
                                                                                             AF_threshold = 0.2)
    } 
    if (random_id == "P0081") { 
      sample_filtered_variants_clusters_df <- splitClusterOnTimePointSpecificAlleleFrequency(cluster_df = sample_filtered_variants_clusters_df,
                                                                                             time_point = "Dx",
                                                                                             cluster = "ClDx",
                                                                                             AF_threshold = 0.24)
      sample_filtered_variants_clusters_df <- splitClusterOnTimePointSpecificAlleleFrequency(cluster_df = sample_filtered_variants_clusters_df,
                                                                                             time_point = "R2",
                                                                                             cluster = "ClR2",
                                                                                             AF_threshold = 0.19)
    } 
    if (random_id == "P0098") { 
      sample_filtered_variants_clusters_df <- splitClusterOnTimePointSpecificAlleleFrequency(cluster_df = sample_filtered_variants_clusters_df,
                                                                                             time_point = "R1",
                                                                                             cluster = "ClR1",
                                                                                             AF_threshold = 0.2)
    } 
    if (random_id == "P0115") { 
      sample_filtered_variants_clusters_df <- splitClusterOnTimePointSpecificAlleleFrequency(cluster_df = sample_filtered_variants_clusters_df,
                                                                                             time_point = "R1",
                                                                                             cluster = "ClR1",
                                                                                             AF_threshold = 0.2)
    } 
    if (random_id == "P0121") { 
      sample_filtered_variants_clusters_df <- splitClusterOnTimePointSpecificAlleleFrequency(cluster_df = sample_filtered_variants_clusters_df,
                                                                                             time_point = "Dx",
                                                                                             cluster = "ClDx",
                                                                                             AF_threshold = 0.19) 
      sample_filtered_variants_clusters_df <- splitClusterOnTimePointSpecificAlleleFrequency(cluster_df = sample_filtered_variants_clusters_df,
                                                                                             time_point = "R2",
                                                                                             cluster = "ClR2",
                                                                                             AF_threshold = 0.13) 
    }
    if (random_id == "P0140") { 
      sample_filtered_variants_clusters_df <- splitClusterOnTimePointSpecificAlleleFrequency(cluster_df = sample_filtered_variants_clusters_df,
                                                                                             time_point = "Dx",
                                                                                             cluster = "ClDx",
                                                                                             AF_threshold = 0.13) 
      sample_filtered_variants_clusters_df <- splitClusterOnTimePointSpecificAlleleFrequency(cluster_df = sample_filtered_variants_clusters_df,
                                                                                             time_point = "R1",
                                                                                             cluster = "ClR1",
                                                                                             AF_threshold = 0.2) 
    }
  } 
  
  # Delete small clusters
  if (deleteSmallClusters) {
    small_clusters <- names(which(table(sample_filtered_variants_clusters_df$cluster) < MinClustSize))
    sample_filtered_variants_clusters_df <- sample_filtered_variants_clusters_df[!sample_filtered_variants_clusters_df$cluster %in% small_clusters, ]
    
    # Correct NA count after deleting the small ones
    NA_clusters_full <- sample_filtered_variants_clusters_df$cluster[grepl("NA", sample_filtered_variants_clusters_df$cluster)]
    NA_clusters <- mixedsort(unique(NA_clusters_full))
    if (length(NA_clusters) > 0){
      NA_clusters_renamed <- paste0("NA_0", 1:length(NA_clusters))
      names(NA_clusters_renamed) <- NA_clusters
      sample_filtered_variants_clusters_df$cluster[grepl("NA", sample_filtered_variants_clusters_df$cluster)] <-
        NA_clusters_renamed[match(NA_clusters_full, names(NA_clusters_renamed))]
    }
    
    # In case of the deletion of a small fall cluster, relabel the names
    fall_clusters_str <- paste0(c("ClDx", "ClR1", "ClR2", "ClR3", "ClR4", "ClR5"), "_fall")
    for (fall_cluster in fall_clusters_str) {
      fall_cluster_names <- unique(sample_filtered_variants_clusters_df$cluster[grep(fall_cluster, sample_filtered_variants_clusters_df$cluster)])
      
      # Rename fall clusters
      if (length(fall_cluster_names) == 1) {
        sample_filtered_variants_clusters_df$cluster <- gsub(fall_cluster_names, fall_cluster, sample_filtered_variants_clusters_df$cluster)
      } else if (length(fall_cluster_names) > 1){
        for (n_clust in (1:length(fall_cluster_names))) {
          newname <- paste0(fall_cluster, "_",  c("A", "B", "C", "D", "E", "F")[n_clust])
          sample_filtered_variants_clusters_df$cluster <- gsub(fall_cluster_names[n_clust], newname, 
                                                               sample_filtered_variants_clusters_df$cluster)
        }
      } else {
        next
      }
    }
    
    # Rename some clusters manually
    if (random_id == "P0614") {
      sample_filtered_variants_clusters_df$cluster <- gsub("NA_01", "ClDx_sub_C",
                                                           x = sample_filtered_variants_clusters_df$cluster) 
      sample_filtered_variants_clusters_df$cluster <- gsub("ClR1_sub_A", "ClR1_sub",
                                                           x = sample_filtered_variants_clusters_df$cluster) 
      
    } else if (random_id == "P0618") {
      sample_filtered_variants_clusters_df$cluster <- gsub("ClDx", "ClDx_A",
                                                           x = sample_filtered_variants_clusters_df$cluster) 
      sample_filtered_variants_clusters_df$cluster <- gsub("NA_01", "ClDx_B",
                                                           x = sample_filtered_variants_clusters_df$cluster) 
    }
  }

  # Order the cluster types
  orig_cluster_order <- unique(sample_filtered_variants_clusters_df$cluster)
  correct_cluster_order <- cluster_order[cluster_order %in% orig_cluster_order]
  sample_filtered_variants_clusters_df_ordered <- data.frame()
  for (cluster_type in correct_cluster_order) {
    sample_filtered_variants_clusters_df_ordered <- rbind(sample_filtered_variants_clusters_df_ordered, 
                                                          sample_filtered_variants_clusters_df[sample_filtered_variants_clusters_df$cluster == cluster_type,])
  }
  
  # Check that no clusters are missed
  if (!all(unique(sample_filtered_variants_clusters_df$cluster) %in% unique(sample_filtered_variants_clusters_df_ordered$cluster))){
    print(paste0(cluster, " was not found after ordering of the clusters. Check variable cluster_order if it contains this cluster."))
    stop()
  }
  
  # Split the vcf based on the clustering
  sample_filtered_variants_clusters_df <- sample_filtered_variants_clusters_df_ordered
  sample_filtered_variants_clusters_splitted_vcf <- splitVcfBasedOnCluster(sample_filtered_vcf,
                                                                           sample_filtered_variants_clusters_df)
  
  # Create a mutation matrix
  library(ref_genome, character.only = TRUE)
  sample_mut_mat <- mut_matrix(sample_filtered_variants_clusters_splitted_vcf,
                               ref_genome = ref_genome)
  colnames(sample_mut_mat) <- paste0(skion_number, "_", colnames(sample_mut_mat))
  
  # Add variants and mutation matrix to object
  filtered_variants_clusters_df_list_split[[i]] <- sample_filtered_variants_clusters_df
  names(filtered_variants_clusters_df_list_split)[i] <- skion_number
  MR_mut_mat <- cbind(MR_mut_mat, sample_mut_mat)
}
MR_mut_mat <- MR_mut_mat[, colnames(MR_mut_mat) != ""]
rm(filtered_vcf_objects_list)


## Create histograms of VAFs of each cluster at each timepoint
pdf(paste0(output_dir, "/clusterVAFdensities_smallClustersDeleted_afterSplitAndReordering.pdf"), width = 27, height = 25)
for (i in 1:length(filtered_variants_clusters_df_list)){ 
  skion_number <- names(filtered_variants_clusters_df_list)[i]
  random_id <- random_id_info[random_id_info$skion_nr == skion_number, 2]
  
  # Parse and order column names for readability
  sample_variant_clusters <- filtered_variants_clusters_df_list_split[[i]]
  colnames(sample_variant_clusters) <- gsub(skion_number, "", colnames(sample_variant_clusters))
  colnames(sample_variant_clusters) <- gsub("_read_based_AF", "", colnames(sample_variant_clusters))
  
  # Order data frame by tumor
  sample_variant_clusters <- sample_variant_clusters[, c(mixedsort(colnames(sample_variant_clusters)[1:(ncol(sample_variant_clusters) - 1)]), "cluster")]
  
  # Delete small clusters
  small_clusters <- names(which(table(sample_variant_clusters$cluster) < MinClustSize))
  sample_variant_clusters <- sample_variant_clusters[!sample_variant_clusters$cluster %in% small_clusters, ]
  
  # Set factor levels for desired plot order
  sample_variant_clusters$cluster <- factor(sample_variant_clusters$cluster,
                                            levels = unique(sample_variant_clusters$cluster))
  
  # Generate plots
  generate_vaf_overviews(variant_clusters_df = sample_variant_clusters)
}
dev.off()



#--------------------------------- Strict refit -------------------------------#

## Define input data
# Note: combined_mut_mat contains the same info as MR_mut_mat, which was 
# generated previously. Ref cohort samples are just added to it
combined_nmf_res <- nmf_res
combined_mut_mat <- MRsamples_ref_cohort_mut_mat

# Check if the generated mut matrix of MR samples is identical to the 
# previously loaded mut matrix object of MR samples
#! Should be TRUE, if not, check which samples are missing
identical(MR_mut_mat, combined_mut_mat[,colnames(MR_mut_mat)])

# Parse SBSblood
rownames(SBSblood) <- SBSblood$Signature
SBSblood <- SBSblood[, colnames(SBSblood) != "Signature"]
SBSblood <- t(SBSblood)[,"SBSblood", drop = FALSE]
rownames(SBSblood) <- paste0(substring(rownames(SBSblood), 1,1),
                             "[", substring(rownames(SBSblood), 2,2),
                             ">", substring(rownames(SBSblood), 6,6),
                             "]", substring(rownames(SBSblood), 7,7))

# Load COSMIC signatures and add SBSblood to the matrix
signatures <- get_known_signatures()
signatures <- cbind(signatures, SBSblood)

# Rename and reorder the de novo extracted signatures
renamed_nmf_res <- rename_nmf_signatures(combined_nmf_res, 
                                         signatures, 
                                         cutoff = 0.9)
renamed_nmf_res$signatures <- renamed_nmf_res$signatures[,c(4,2,5,8,1,6,3,7)]
renamed_nmf_res$contribution <- renamed_nmf_res$contribution[c(4,2,5,8,1,6,3,7),]

# Rename SBSA to SBSB and SBSB to SBSA to avoid confusion later on
colnames(renamed_nmf_res$signatures)[c(6,7)] <- c("SBSA", "SBSB")
rownames(renamed_nmf_res$contribution)[c(6,7)] <- c("SBSA", "SBSB")


# Show profile plot of de novo signatures
de_novo_signatures <- plot_96_profile(renamed_nmf_res$signatures,
                ymax = 0.3,
                condensed = TRUE) + 
  xlab("Sequence context") + 
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 11),
        strip.text = element_text(size = 12))
de_novo_signatures

# Scale the denovo signatures so that the sum is 1
signatures_denovo_scaled <- renamed_nmf_res$signatures
for (n_col in 1:ncol(signatures_denovo_scaled)){
  signatures_denovo_scaled[, n_col] <- renamed_nmf_res$signatures[, n_col] / 
                                       sum(renamed_nmf_res$signatures[, n_col])
}
renamed_nmf_res$signatures <- signatures_denovo_scaled

# Determine what the de novo signatures contain by doing a strict refit
strict_refit <- fit_to_signatures_strict(renamed_nmf_res$signatures, 
                                         signatures, 
                                         max_delta = 0.033)
fit_res_strict <- strict_refit$fit_res
pdf(paste0(output_dir, "/strictRefit_signatureContributions.pdf"), width = 7, height = 5)
plot_contribution(fit_res_strict$contribution[],
                  coord_flip = FALSE,
                  mode = "relative") + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.7))
dev.off()

# Determine cos sim between de novo signatures and reconstructed profile from refit with COSMIC
pdf(paste0(output_dir, "/strictRefit_reconstrCosSim.pdf"), width = 7, height = 5)
plot_original_vs_reconstructed(renamed_nmf_res$signatures[], 
                               fit_res_strict$reconstructed[],
                               y_intercept = 0.9) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.7))
dev.off()

# Print cosine similarities
cos_sim_reconstr <- as.data.frame(diag(cos_sim_matrix(renamed_nmf_res$signatures[], 
                                                      fit_res_strict$reconstructed[])))
cos_sim_cosmic <- cos_sim_matrix(renamed_nmf_res$signatures[],
                                 signatures)
cos_sim_cosmic[, colnames(cos_sim_cosmic)[apply(cos_sim_cosmic, 1, which.max)]]

## Select the collection of signatures. 
# Use matched COSMIC signature if the de novo signature has a >=90% cos sim with one or two COSMIC signatures (e.g. SBS2/SBS13)
# Conclusion: SBS1, SBS2, SBS13, SBS7a, SBS18, SBS86, SBS87 selected, still left with 2 de novo signatures (SBSA and SBSB)
# SBSA had good reconstructed vs original cos sim but contained too many signatures (7) 
# SBSB did not reach cos sim threshold of 0.9 to SBS3 

## Calculate relative contribution
# Function to calculate the relative contribution
calculate_relative_contribution <- function(contribution_col){
  contributions <- contribution_col[1:(length(contribution_col) - 1)]
  total <- contribution_col[length(contribution_col)]
  relative_contributions <- contributions / total
  return(relative_contributions)
}
contribution_with_totals <- rbind(renamed_nmf_res$contribution, 
                                  apply(renamed_nmf_res$contribution, 2, sum))
rownames(contribution_with_totals)[nrow(contribution_with_totals)] <- "Total"
relative_contributions <- apply(contribution_with_totals, 2, calculate_relative_contribution)

# Obtain samples with rel contr >=35% to SBSA and SBSB
rel_contr_data <- as.data.frame(t(relative_contributions[,colSums(combined_mut_mat) > 200]))
SBSA_samples <- rownames(rel_contr_data)[rel_contr_data$SBSA >= 0.35]
SBSB_samples <- rownames(rel_contr_data)[rel_contr_data$SBSB >= 0.35]
rownames(rel_contr_data) <- gsub("^P_", "", rownames(rel_contr_data))
for (i in 1:nrow(rel_contr_data)){
  rowname <- gsub("_Cl", "-Cl", rownames(rel_contr_data)[i])
  cluster <- unlist(strsplit(rowname, "-"))[2]
  sample <- compl_random_id_info$random_id[match(unlist(strsplit(rowname, "-"))[1], compl_random_id_info$skion_nr)]
  rownames(rel_contr_data)[i] <- ifelse(!is.na(cluster), paste0(sample, "_", cluster), sample)
}
rownames(rel_contr_data) <- gsub("_ClR1_sub", "\nClR1_sub", rownames(rel_contr_data))
random_SBSA_samples <- rownames(rel_contr_data)[rel_contr_data$SBSA >= 0.35]
random_SBSB_samples <- rownames(rel_contr_data)[rel_contr_data$SBSB >= 0.35]

# Generate plot for samples with > 200 muts
custom_cols <- c("#9BC0CD", "#F6EB13","#AA6AAC", "#0A0F25", "#5C6BC0", "#DADADA","#FF8FAB","#B12325")
sbsa_plot <- plot_contribution(t(rel_contr_data)[,random_SBSA_samples],
                  coord_flip = FALSE,
                  mode = "relative") +
  facet_grid(. ~ "SBSA") +
  scale_fill_manual(values = custom_cols) +
  theme(axis.text.x = element_text(angle = -90, vjust = 0.7, hjust = 0), legend.position = "none")
sbsb_plot <- plot_contribution(t(rel_contr_data)[,random_SBSB_samples],
                  coord_flip = FALSE,
                  mode = "relative") +
  scale_fill_manual(values = custom_cols) +
  facet_grid(. ~ "SBSB") +
  theme(axis.text.x = element_text(angle = -90, vjust = 0.7, hjust = 0),
        axis.title.y = element_blank())

# Save to PDF
pdf(paste0(output_dir, "/strictRefitSBSAandSBSB.pdf"), width = 7, height = 6)
cowplot::plot_grid(sbsa_plot, sbsb_plot, rel_widths = c(1,0.55))
dev.off()

# Redo strict refit on selected samples
SBSA_samples_strict_refit <- fit_to_signatures_strict(combined_mut_mat[,colnames(combined_mut_mat) %in% SBSA_samples], 
                                                      signatures, 
                                                      max_delta = 0.033)
SBSA_samples_strict_refit <- SBSA_samples_strict_refit$fit_res
colnames(SBSA_samples_strict_refit$contribution) <- random_SBSA_samples
SBSB_samples_strict_refit <- fit_to_signatures_strict(combined_mut_mat[,colnames(combined_mut_mat) %in% SBSB_samples], 
                                                      signatures, 
                                                      max_delta = 0.033)
SBSB_samples_strict_refit <- SBSB_samples_strict_refit$fit_res
colnames(SBSB_samples_strict_refit$contribution) <- random_SBSB_samples

# Show relative contributions of strict refit results
sbsa_plot <- plot_contribution(SBSA_samples_strict_refit$contribution[],
                  coord_flip = FALSE,
                  mode = "relative") + 
  facet_grid(. ~ "SBSA") +
  theme(axis.text.x = element_text(angle = -90, vjust = 0.7, hjust = 0))
sbsb_plot <- plot_contribution(SBSB_samples_strict_refit$contribution[],
                  coord_flip = FALSE,
                  mode = "relative") + 
  facet_grid(. ~ "SBSB") +
  theme(axis.text.x = element_text(angle = -90, vjust = 0.7, hjust = 0),
        axis.title.y = element_blank())

# Save to PDF
pdf(paste0(output_dir, "/strictRefitSBSAandSBSB_relativeContributions.pdf"), width = 8, height = 6)
cowplot::plot_grid(sbsa_plot, sbsb_plot, rel_widths = c(1,0.38))
dev.off()


## Conclusions about SBSA and SBSB
#
# SBSA found in multiple patients and is explained by multiple different signatures 
# in each patient, i.e. cannot be explained by COSMIC signatures. Add to data set as it is.
#
# SBSB found in 2 samples of 1 patient. Samples contain signatures SBS17a/b and SBS86/87. 
# SBS86/87 were previously selected, SBS17a/b will be added to the data set instead of SBSB

# Construct the final set of selected signatures
selected_signatures <- signatures[,c("SBS1","SBS2","SBS13","SBS7a","SBS17a","SBS17b","SBS18","SBS86","SBS87")]
selected_signatures <- cbind(selected_signatures, renamed_nmf_res$signatures[,"SBSA", drop = FALSE])

# Show profile plot of signature selection
final_signatures <- plot_96_profile(selected_signatures,
                                      ymax = 0.4,
                                      condensed = TRUE) + 
  xlab("Sequence context") + 
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 11),
        strip.text = element_text(size = 12))
final_signatures

# Save to PDF
pdf(paste0(output_dir, "/deNovoSignatures_vs_selectedSignatures_profilePlots_newNMF.pdf"), width = 12, height = 10)
grid.arrange(de_novo_signatures, final_signatures, ncol = 2)
dev.off()

# Final strict refit
strict_refit_final <- fit_to_signatures_strict(MR_mut_mat, 
                                               selected_signatures, 
                                               max_delta = 0.033)
fit_res_strict_final <- strict_refit_final$fit_res
fit_res_strict_final$signatures <- selected_signatures

# Visualize strict refit results
sample_order <- names(sort(colSums(MR_mut_mat), decreasing = TRUE))
contr_data <- fit_res_strict_final$contribution[,sample_order]
colnames(contr_data) <- paste0(colnames(contr_data), " (n = ", colSums(MR_mut_mat[,sample_order]), ")")
for (i in 1:ncol(contr_data)){
  colname <- gsub("_Cl", "-Cl", colnames(contr_data)[i])
  cluster <- unlist(strsplit(colname, "-"))[2]
  sample <- random_id_info$random_id[match(unlist(strsplit(colname, "-"))[1], random_id_info$skion_nr)]
  colnames(contr_data)[i] <- paste(sample, cluster)
}
reconstr_cossim_plot <- plot_original_vs_reconstructed(MR_mut_mat[,rev(sample_order)], 
                                                  fit_res_strict_final$reconstructed[,rev(sample_order)],
                                                  y_intercept = 0.85) + coord_flip() +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01),
                     limits = c(0,1),
                     expand = c(0.001, 0.001)) +
  ylab("Original vs reconstructed cos. sim.") + 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(size = 8))
sig_heatmap <- plot_contribution_heatmap(contr_data,
                                         cluster_samples = FALSE, 
                                         cluster_sigs = FALSE) +
  theme(legend.position = "left", 
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8),
        axis.text = element_text(size = 7.5)) 
layout_mat <- t(matrix(c(rep(1, 9), rep(2, 4))))
grid.arrange(sig_heatmap,
             reconstr_cossim_plot, 
             layout_matrix = layout_mat)



## Visualize strict refit result effects on SBSB samples
SBSB_samples_mut_mat <- combined_mut_mat[,SBSB_samples]

# Strict refit with SBSB instead of SBSB17a/b
SBSB_strict_refit <- fit_to_signatures_strict(combined_mut_mat[,SBSB_samples], 
                                              cbind(selected_signatures[,!colnames(selected_signatures) %in% c("SBS17a", "SBS17b")],
                                                    renamed_nmf_res$signatures[,"SBSB", drop = FALSE]), 
                                              max_delta = 0.033)$fit_res$reconstructed
# Strict refit with SBSB17a/b instead of SBSB
SBS17_strict_refit <- fit_to_signatures_strict(combined_mut_mat[,SBSB_samples], 
                                               selected_signatures, 
                                               max_delta = 0.033)$fit_res$reconstructed

# Determine original vs reconstructed cosine similarity
colnames(SBSB_strict_refit) <- random_SBSB_samples
colnames(SBS17_strict_refit) <- random_SBSB_samples
colnames(SBSB_samples_mut_mat) <- random_SBSB_samples
SBSB_plot <- plot_original_vs_reconstructed(SBSB_strict_refit, 
                                            SBSB_samples_mut_mat,
                                            y_intercept = 0.85) +
  labs(title = "Using SBSB") +
  #facet_grid(. ~ "Using SBSB") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.7),
        axis.title.y = element_text(size = 10.5),
        plot.title = element_text(size = 11, face = "bold", hjust = 0.5))
#strip.background = element_rect(fill = "grey85", color = "grey85"))
SBS17_plot <- plot_original_vs_reconstructed(SBS17_strict_refit, 
                                           SBSB_samples_mut_mat,
                                           y_intercept = 0.85) +
  labs(title = "Using SBS17a and SBS17b") +
  #facet_grid(. ~ "Using SBS17a/SBS17b") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.7),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 11, face = "bold", hjust = 0.5))
#strip.background = element_rect(fill = "grey85", color = "grey85"))

# Save to PDF
pdf(paste0(output_dir, "/strictRefitSBSBsamples_reconstrCosSim_v3.pdf"), width = 5, height = 6)
cowplot::plot_grid(SBSB_plot, SBS17_plot, rel_widths = c(1,0.89))
dev.off()





#-------------------------- Contribution calculations -------------------------#

## Add unassigned signature, clusters smaller than 200 will be assigned 'unassigned'
# Calculate the contribution of Unassigned clusters
fit_res_strict_final_unassigned <- fit_res_strict_final
Unassigned <- ifelse(colSums(MR_mut_mat) < 200, colSums(fit_res_strict_final_unassigned$contribution), 0)
fit_res_strict_final_unassigned$contribution <- rbind(fit_res_strict_final_unassigned$contribution,
                                                      Unassigned)

# Replace contr. of signatures to 0 if the cluster contains < 200 mutations
n_sign <- ncol(selected_signatures)
for (i in 1:ncol(fit_res_strict_final_unassigned$contribution)){
  if (fit_res_strict_final_unassigned$contribution["Unassigned", i] > 0){
    fit_res_strict_final_unassigned$contribution[1:n_sign, i] <- rep(0, n_sign)
  } else {
    next
  }
}
fit_res_strict_final <- fit_res_strict_final_unassigned

#load(list.files(pattern = "fit_res_strict_final_17032024.rdata", full.names = TRUE, recursive = TRUE))

## Calculate relative and absolute contributions
# Relative
contribution_with_totals <- rbind(fit_res_strict_final$contribution, 
                                  apply(fit_res_strict_final$contribution, 2, sum))
rownames(contribution_with_totals)[nrow(contribution_with_totals)] <- "Total"
relative_contributions <- apply(contribution_with_totals, 2, calculate_relative_contribution)

# Absolute
calculate_absolute_contribution <- function(contribution_col){
  contributions <- contribution_col[1:(length(contribution_col) - 2)]
  total <- contribution_col[length(contribution_col) - 1]
  originaltotal <- contribution_col[length(contribution_col)]
  absolute_contributions <- ((contributions / total) * originaltotal)
  return(absolute_contributions)
}
contribution_with_totals <- rbind(contribution_with_totals, 
                                  apply(MR_mut_mat, 2, sum))
rownames(contribution_with_totals)[nrow(contribution_with_totals)] <- "OriginalTotal"
absolute_contributions <- round(apply(contribution_with_totals, 2, calculate_absolute_contribution))

# Add reconstructed cosine similarity to data frame
Reconstruction <- diag(cos_sim_matrix(MR_mut_mat, 
                                      fit_res_strict_final$reconstructed[]))
absolute_contributions_withreconstr <- rbind(absolute_contributions, Reconstruction)

# Calculate bootstrap values
contri_bootstrap <- fit_to_signatures_bootstrapped(MR_mut_mat,
                                                   fit_res_strict_final$signatures,
                                                   n_boots = 100,
                                                   method = "strict",
                                                   max_delta = 0.033)

# Calculate bootstrap percentage for each signature
contributions_withbootstrap <- rbind(absolute_contributions_withreconstr, 
                                     matrix(ncol = ncol(absolute_contributions_withreconstr), 
                                            nrow = ncol(contri_bootstrap), 
                                            dimnames = list(paste0(colnames(contri_bootstrap), "_bootstrapPerc"), 
                                                            colnames(absolute_contributions_withreconstr))))
for (patient in colnames(fit_res_strict_final$contribution)){
  patient_bootstrap <- contri_bootstrap[str_remove(rownames(contri_bootstrap), "_\\d+$") == patient,]
  signature_percentage <- colSums(patient_bootstrap > 0)
  contributions_withbootstrap[grep("bootstrap", rownames(contributions_withbootstrap)), patient] <- signature_percentage
}

# Write info to out file, incl sample and cluster info
contributions_withbootstrap_df <- as.data.frame(t(contributions_withbootstrap))
rownames(contributions_withbootstrap_df) <- gsub("_Cl", "-Cl", rownames(contributions_withbootstrap_df))
sample_info <- data.frame(Patient_ID = sapply(strsplit(rownames(contributions_withbootstrap_df), "-"), `[`, 1),
                          Cluster = sapply(strsplit(rownames(contributions_withbootstrap_df), "-"), `[`, 2), 
                          row.names = rownames(contributions_withbootstrap_df))
sample_info$Patient_ID <- random_id_info$random_id[match(sample_info$Patient_ID, random_id_info$skion_nr)]
contributions_withbootstrap_df <- cbind(sample_info, contributions_withbootstrap_df)
write.table(contributions_withbootstrap_df, quote = FALSE, sep = "\t",
          paste0(output_dir, "/Refit_multiplerelapses_includingbootstrap_17032024.tsv"))


## Change signature contributions to the Unassigned signature if the bootstrap value is below 75%
for (sample in colnames(absolute_contributions)){
  sample_info <- data.frame(Contribution = absolute_contributions[colnames(selected_signatures), sample],
                            Bootstrap = as.data.frame(contributions_withbootstrap)[paste0(colnames(selected_signatures), "_bootstrapPerc"), sample])
  
  # Determine which samples have a signature contr. with a bootstrap value < 75
  sample_info_subset <- subset(sample_info, Contribution > 0 & Bootstrap < 75)
  subset_sign <- rownames(sample_info_subset)
  
  if (length(subset_sign) > 0){
    # Move selected signature contributions to Unassigned signature
    unassigned_contr <- colSums(absolute_contributions[subset_sign, sample, drop = FALSE])
    absolute_contributions["Unassigned", sample] <- unassigned_contr
    
    # Change signature contribution to 0
    absolute_contributions[subset_sign, sample] <- 0
  } else {
    next
  }
}
#load("./mutMatricesAndRefitResults/fit_res_strict_final_bootstrapCorrected_17032024.rdata")



#-------------------------------- Visualisations ------------------------------#

## Visualize VAFs, signatures and contributions in each cluster of each sample
for (i in 1:length(filtered_variants_clusters_df_list)){ 
  skion_number <- names(filtered_variants_clusters_df_list)[i]
  random_id <- random_id_info[random_id_info$skion_nr == skion_number, 2]
  sample_filtered_variants_clusters_df <- filtered_variants_clusters_df_list_split[[i]]
  sample_mut_mat <- MR_mut_mat[, grep(skion_number, colnames(MR_mut_mat)), drop = FALSE]
  sample_abs_contribution <- absolute_contributions[, grep(skion_number, colnames(absolute_contributions)), drop = FALSE]
  
  # Remove sample IDs from data frames
  colnames(sample_mut_mat) <- gsub(paste0(skion_number, "_"), "", colnames(sample_mut_mat))
  colnames(sample_abs_contribution) <- gsub(paste0(skion_number, "_"), "", colnames(sample_abs_contribution))
  
  
  
  ##### Create variant AF plot with CNV information #####
  # Melt variant info
  sample_filtered_variants_clusters_df$variant <- rownames(sample_filtered_variants_clusters_df)
  sample_filtered_variants_clusters_df_melt <- reshape2::melt(sample_filtered_variants_clusters_df, 
                                                              id = c("cluster", "variant"))
  names(sample_filtered_variants_clusters_df_melt) <- c("cluster", "variant", "time_point", "AF")
  
  # Parse info to increase readability
  sample_filtered_variants_clusters_df_melt$time_point <- gsub("_read_based_AF", "",
                                                               sample_filtered_variants_clusters_df_melt$time_point) 
  sample_filtered_variants_clusters_df_melt$time_point <- gsub(skion_number, "",
                                                               sample_filtered_variants_clusters_df_melt$time_point) 
  
  ## Annotate variants with CNV information
  # Read CNV file
  cnv_file <- cnv_files[grepl(skion_number, cnv_files)] 
  print(paste0("Opening CNV file: ", cnv_file))
  cnv_sample <- read.delim(cnv_file)
  CNVs <- dplyr::select(cnv_sample, all_of(c('CHR', 'CNV_START', 'CNV_STOP', 'CNV_CALL')))
  CNVs <- distinct(CNVs)

  # Check for each variant if it is located in a CNV
  for (variantnum in 1:nrow(sample_filtered_variants_clusters_df_melt)){
    loc_vec <- unlist(strsplit(sample_filtered_variants_clusters_df_melt$variant[variantnum],split = "[:_]+"))
    chromosome <- loc_vec[1]
    position <- as.integer(loc_vec[2])
    list_of_bools <- c()
    type <- "NO CNA"
    
    if (chromosome %in% unique(CNVs$CHR)) {
      chrom_CNV_subset <- CNVs[CNVs$CHR==chromosome,]

      for (cnv in 1:nrow(chrom_CNV_subset)) {
        in_cnv <- (position >= chrom_CNV_subset[cnv,]$CNV_START) & (position <= chrom_CNV_subset[cnv,]$CNV_STOP)
        list_of_bools <- append(list_of_bools, in_cnv)
        if (in_cnv == TRUE) {
          type <- chrom_CNV_subset[cnv,]$CNV_CALL
        }
      }
    } 
    max_AF_variant <- max(sample_filtered_variants_clusters_df_melt[sample_filtered_variants_clusters_df_melt$variant == 
                                                                    sample_filtered_variants_clusters_df_melt$variant[variantnum],]$AF)
    if (type == "NO CNA" && chromosome == "chrX" && 
        skion_gender_df[skion_gender_df$SKION == skion_number,]$Gender == "M" && 
        max_AF_variant > 0.7) {
      type <- "CHR X"
    }
    
    if(!(chromosome %in% unique(CNVs$CHR))){
      sample_filtered_variants_clusters_df_melt$In_CNV[variantnum] <- FALSE
    }
    else {
      sample_filtered_variants_clusters_df_melt$In_CNV[variantnum] <- sum(list_of_bools)>=1
    }
    sample_filtered_variants_clusters_df_melt$Type[variantnum] <- type
  }

  # Create AF plot with CNV info
  sample_filtered_variants_clusters_df_melt$cluster <- factor(sample_filtered_variants_clusters_df_melt$cluster,
                                                              levels = unique(sample_filtered_variants_clusters_df$cluster))
  x_width <- length(unique(sample_filtered_variants_clusters_df_melt$time_point))
  x_width <- x_width - 0.45
  sample_af_merged_plots <- ggplot_cnv_af(sample_filtered_variants_clusters_df_melt)
  
  
  
  ##### Create profile plot ######
  # Define colors used for the signatures
  CUSTOM_COLOUR_PALETTE <- c("#9BC0CD", "#CE2627", "#B12325", "#F6EB13", "#FFC2D1","#FF8FAB", 
                             "#AA6AAC", "#0A0F25", "#5C6BC0", "#DADADA", "#A9A9A9")
  names(CUSTOM_COLOUR_PALETTE) <- c("SBS1","SBS2","SBS13","SBS7a","SBS17a","SBS17b", 
                                    "SBS18","SBS86","SBS87", "SBSA", "Unassigned")
  
  # Adjust ymax when necessary
  if (random_id %in% c("P0077", "P0122", "P0174", "P0180", "P0609")) {
    ymax_custom <- 0.3
  } else {
    ymax_custom <- 0.2
  }
  
  # Generate profile plot
  colnames(sample_mut_mat) <- factor(colnames(sample_mut_mat), levels = colnames(sample_mut_mat))
  sample_profile_plots <- plot_96_profile(sample_mut_mat,
                                          ymax = ymax_custom,
                                          condensed = T) + 
    xlab("Sequence context") + 
    theme(axis.text.y = element_text(size = 10, color = "black"),
          axis.text.x = element_text(size = 6, color = c("black", "transparent")),
          axis.title = element_text(size = 11.5),
          strip.text.x = element_text(size = 10),
          strip.text.y = element_text(size = 10))
  

  
  ##### Create absolute contribution plot #####
  # Limit y-axis label size and set factor levels to increase readability
  colnames(sample_abs_contribution) <- gsub("_fall", "\nfall", colnames(sample_abs_contribution))
  colnames(sample_abs_contribution) <- gsub("_sub", "\nsub", colnames(sample_abs_contribution))
  colnames(sample_abs_contribution) <- factor(colnames(sample_abs_contribution), levels = colnames(sample_abs_contribution))
  
  # Generate absolute contribution plot
  sample_palette <- CUSTOM_COLOUR_PALETTE[names(CUSTOM_COLOUR_PALETTE) %in% rownames(sample_abs_contribution)]
  absolute_contribution_sample <-  plot_contribution_custom(contribution = sample_abs_contribution,
                                                            signatures = fit_res_strict_final$signatures,
                                                            mode = "absolute",
                                                            coord_flip = TRUE,
                                                            palette = sample_palette)
  
  # Combine all generated plots and save to PDF
  pdf_height <- ifelse(random_id == "P0615", 3.5, 8.5)
  fn <- paste0(output_dir, "/P_", skion_number, "-denovo_overview.pdf")
  pdf(file = fn, width = 18, height = pdf_height)
  gridExtra::grid.arrange(sample_af_merged_plots,
                          sample_profile_plots,
                          absolute_contribution_sample,
                          top = textGrob(random_id, gp = gpar(fontsize = 18, font = 2)),
                          layout_matrix = matrix(c(1,1,1,1,1,1,1,
                                                   2,2,2,2,2,2,
                                                   3,3,3,3,3,3,3,3), nrow = 1))
  dev.off()
}



#------------------------- SNV signature probabilities ------------------------#

# Read input file
snv_info <- as.data.frame(read_excel("2023_07_13_Overview of SNV signature probabilities.xlsx")[,c(1:35)])

# Rescale the contribution to [0,1] for each extracted signature
scaled_contributions <- reshape2::melt(fit_res_strict_final$contribution)
colnames(scaled_contributions) = c("Signature", "Sample", "Contribution")
samples <- unique(scaled_contributions$Sample)
scaled_contributions$Relative_contribution <- NA
for (i in 1:length(samples)){
  sample_specific_slice <- scaled_contributions[scaled_contributions$Sample == samples[i], ]
  sample_specific_slice_total_contribution <- sum(sample_specific_slice$Contribution)
  relative_contribution_vector <- sample_specific_slice$Contribution/sample_specific_slice_total_contribution
  scaled_contributions$Relative_contribution[scaled_contributions$Sample == samples[i]] <- relative_contribution_vector
}

# Reshape the matrix to increase readabiity
scaled_contributions_df <- reshape(scaled_contributions,
                                   idvar = "Sample",
                                   timevar = "Signature",
                                   direction = "wide")
scaled_contributions_df <- scaled_contributions_df[,!grepl("Unassigned", colnames(scaled_contributions_df))]

# convert the counts of each trinucleotide to relative contributions per extracted signature
relative_contribution_trinucleotide_per_signature <- apply(fit_res_strict_final$signatures, 2, function(x) x/sum(x))

# extract the number of groups we have to perform the operation on
number_of_cluster <- nrow(scaled_contributions_df)

# create the empty list that will contain all the data.frames
relative_contribution_list <- c(rep(list(NA), number_of_cluster))
names(relative_contribution_list) <- scaled_contributions_df[,1]

# for reach row (= sample/cluster) we create the relative contribution vector
for (i in 1:number_of_cluster){
  relative_contribution_vector <- as.numeric(scaled_contributions_df[i,seq(3, ncol(scaled_contributions_df), 2)]) 
  
  # multiply vector * matrix
  multiplied_matrix <- sweep(relative_contribution_trinucleotide_per_signature,
                             MARGIN=2,
                             relative_contribution_vector,
                             `*`)
  
  # store the matrix in the "master list"
  relative_contribution_list[[i]] <- multiplied_matrix
}

# create the empty list that will contain all the data.frames
normalized_relative_contribution_list <- c(rep(list(NA), number_of_cluster))
names(normalized_relative_contribution_list) <- scaled_contributions_df[,1]

# for reach row (= sample/cluster) we create the normalized relative contribution vector
for (i in 1:number_of_cluster){
  normalized_relative_contribution_df <- relative_contribution_list[[i]] 
  
  # normalize the values of the data frame
  for (j in 1:nrow(normalized_relative_contribution_df)){
    vector_total <- sum(normalized_relative_contribution_df[j,])
    rescaled_vector <- as.numeric(normalized_relative_contribution_df[j, ])/vector_total
    normalized_relative_contribution_df[j, ] <- rescaled_vector
  }
  
  # save the normalized data frame
  normalized_relative_contribution_list[[i]] <- normalized_relative_contribution_df
}

# Function to reverse complement DNA sequence
reverse_complement_dna_sequence <- function(dnaSequence){
  return(reverse(chartr("ATGC", "TACG", dnaSequence)))
}

# Get probability for each trinucleotide context of selected SNVs
probability_df <- data.frame()
for (i in 1:nrow(snv_info)){ 

  # Extract SNV info
  snv_line <- snv_info[i,]
  snv <- paste0(snv_line$chr, ":", snv_line$start, "_", snv_line$ref, "/", snv_line$alt)
  
  # Get corresponding sample and cluster info
  skion_number <- snv_info$SKION[i]
  sample_snvs <- filtered_variants_clusters_df_list_split[[which(names(filtered_variants_clusters_df_list_split) == skion_number)]]
  snv_cluster <- sample_snvs[snv,]$cluster

  # Get normalized relative contributions belonging to the cluster
  sample_cluster_probabilities <- normalized_relative_contribution_list[[which(names(normalized_relative_contribution_list) == 
                                                                               paste0(skion_number, "_", snv_cluster))]]
  
  # Check which signatures contributions are assigned to Unassigned after bootstrap correction
  present_signatures <- names(which(colSums(sample_cluster_probabilities) > 0))
  corrected_signatures <- names(which(absolute_contributions[,paste0(skion_number, "_", snv_cluster)] > 0))
  remove_signatures <- setdiff(present_signatures, corrected_signatures)
  if (length(remove_signatures) > 0){
    sample_cluster_probabilities[,remove_signatures] <- 0
  }
  
  # Get mutation type and context for the mutation
  mutation_type <- paste0(snv_line$ref, ">", snv_line$alt)
  grl <- makeGRangesFromDataFrame(snv_line[,c("chr", "start", "end", "ref", "alt")],
                                  keep.extra.columns = T)
  GenomeInfoDb::genome(grl) = 'hg38'
  mut_context <- mut_context(grl, ref_genome = ref_genome)
  
  # Get 96-nucleotide context of the mutation
  valid_mutations <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  if (mutation_type %in% valid_mutations){
    mut_context_96 <- paste0(substring(mut_context, first = 1, last = 1), 
                          "[", mutation_type, "]",
                          substring(mut_context, first = 3, last = 3))
  # Get reverse compliment if mutation type is on the negative strand
  } else {
    mut_context_96 <- paste0(reverse_complement_dna_sequence(substring(mut_context, first = 3, last = 3)),
                          "[", reverse_complement_dna_sequence(substring(mutation_type, first = 1, last = 1)), ">",
                          reverse_complement_dna_sequence(substring(mutation_type, first = 3, last = 3)), "]",
                          reverse_complement_dna_sequence(substring(mut_context, first = 1, last = 1)))
  }
  # Get mutation probabilities for the mutation context and add to data frame
  mut_context_probabilities <- sample_cluster_probabilities[mut_context_96,, drop = F]
  print(sum(mut_context_probabilities))
  snv_line_compl <- cbind(snv_line, 
                          data.frame(trinucleotide_context = mut_context,
                                     trinucleotide_mutation = mut_context_96),
                          mut_context_probabilities)
  probability_df <- rbind(probability_df, snv_line_compl)
  
}

probability_df["Unassigned"] <- 1 - round(rowSums(probability_df[,c(38:47)]), digits = 6)
openxlsx::write.xlsx(probability_df, 
                     file = "20240423_Overview_SNV_signature_probabilities_bootstrapCorrectedSignaturesMovedToUnassigned.xlsx")



#-------------------------- Ref cohort contributions --------------------------#

# Prepare data
selected_signatures <- fit_res_strict_final$signatures
ref_cohort_mut_mat <- MRsamples_ref_cohort_mut_mat[, !colnames(MRsamples_ref_cohort_mut_mat) %in% 
                                                      colnames(MR_mut_mat)]

# Strict refit
strict_refit_ref_cohort <- fit_to_signatures_strict(ref_cohort_mut_mat, 
                                                    selected_signatures, 
                                                    max_delta = 0.033)
fit_res_strict_ref_cohort <- strict_refit_ref_cohort$fit_res
fit_res_strict_ref_cohort$signatures <- selected_signatures

## Calculate relative and absolute contributions
# Relative
contribution_with_totals <- rbind(fit_res_strict_ref_cohort$contribution, 
                                  apply(fit_res_strict_ref_cohort$contribution, 2, sum))
rownames(contribution_with_totals)[nrow(contribution_with_totals)] <- "Total"
relative_contributions <- apply(contribution_with_totals, 2, calculate_relative_contribution)

# Absolute
contribution_with_totals <- rbind(contribution_with_totals, 
                                  apply(ref_cohort_mut_mat, 2, sum))
rownames(contribution_with_totals)[nrow(contribution_with_totals)] <- "OriginalTotal"
absolute_contributions <- round(apply(contribution_with_totals, 2, calculate_absolute_contribution))

# Add reconstructed cosine similarity to data frame
Reconstruction <- diag(cos_sim_matrix(ref_cohort_mut_mat, 
                                      fit_res_strict_ref_cohort$reconstructed[]))
absolute_contributions_withreconstr <- rbind(absolute_contributions, Reconstruction)

# Calculate bootstrap values
contri_bootstrap <- fit_to_signatures_bootstrapped(ref_cohort_mut_mat,
                                                   fit_res_strict_ref_cohort$signatures,
                                                   n_boots = 100,
                                                   method = "strict",
                                                   max_delta = 0.033)

# Calculate bootstrap percentage for each signature
contributions_withbootstrap <- rbind(absolute_contributions_withreconstr, 
                                     matrix(ncol = ncol(absolute_contributions_withreconstr), 
                                            nrow = ncol(contri_bootstrap), 
                                            dimnames = list(paste0(colnames(contri_bootstrap), "_bootstrapPerc"), 
                                                            colnames(absolute_contributions_withreconstr))))
for (patient in colnames(fit_res_strict_ref_cohort$contribution)){
  patient_bootstrap <- contri_bootstrap[str_remove(rownames(contri_bootstrap), "_\\d+$") == patient,]
  signature_percentage <- colSums(patient_bootstrap > 0)
  contributions_withbootstrap[grep("bootstrap", rownames(contributions_withbootstrap)), patient] <- signature_percentage
}

# Write info to out file, incl sample and cluster info
contributions_withbootstrap <- as.data.frame(t(contributions_withbootstrap))
rownames(contributions_withbootstrap) <- gsub("P_", "", rownames(contributions_withbootstrap))
sample_info <- ref_random_id_info[rownames(contributions_withbootstrap),]
contributions_withbootstrap <- cbind(sample_info, contributions_withbootstrap)
write.table(contributions_withbootstrap, quote = FALSE, sep = "\t",
            paste0(output_dir, "/Refit_refCohort_includingbootstrap_25032024.tsv"))




