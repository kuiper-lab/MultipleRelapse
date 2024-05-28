#!/usr/bin/env Rscript 

# Clear data 
rm(list = ls())

# Load packages
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)
library(MutationalPatterns)
library(stringr)
library(tidyr)
library(tidyverse)
library(readxl)
library(plyr)
library(lsa)
library(gtools)
library(cowplot)
library(grid)
library(gridExtra)
library(ggpubr)

# Set global parameters
# All input files should be located somewhere in PROJECT_DIR
PROJECT_DIR <- "/path/to/project/dir/"
RESULTS_DIR <- paste0(PROJECT_DIR, "/path/to/result/dir/")
RANDOM_ID_INFO <- "MR_cohortDescription.xlsx"
CENTROMERES_FILE <- "cytoBand.hg38.centromeresOnly.txt"
DEEPSEQ_INFO <- "deepsequencing_32039variants_depthAndAFinfo.txt"
SNV_EXT <- "snp-WGS.vepAnnotated.vcf.gz"

setwd(RESULTS_DIR)

# Load custom functions
FUNCTIONS_SCRIPT <- list.files(path = "/path/to/code/dir/", 
                               pattern = "createdFunctions_HM-ALL.R",
                               full.names = TRUE, recursive = TRUE)
source(FUNCTIONS_SCRIPT)

# Define filter parameters
REMOVE_GERMLINE_VARIANTS <- FALSE
DEEPSEQUENCING <- TRUE
split <- TRUE
deleteSmallClusters <- TRUE
MIN_AF <- 0.25
MinClustSize <- 75
noise_threshold <- 0  
noise_threshold_deep <- 0.01

# Obtain paths to input files 
centromeres <- list.files(path = PROJECT_DIR,
                          pattern = CENTROMERES_FILE, 
                          full.names = TRUE, 
                          recursive = TRUE)
vcf_vector <- list.files(path = PROJECT_DIR,
                         pattern = SNV_EXT, 
                         full.names = TRUE, 
                         recursive = TRUE)

# Read deep sequencing and random ID info
ds_data <- read.table(file = list.files(pattern = DEEPSEQ_INFO, 
                                        full.names = TRUE, 
                                        recursive = TRUE),
                      header = TRUE, sep = "\t")
random_id_info <- as.data.frame(read_excel(list.files(path = PROJECT_DIR,
                                                      pattern = RANDOM_ID_INFO,
                                                      recursive = TRUE,
                                                      full.names = TRUE))[1:29, c(1,2)])
colnames(random_id_info) <- c("skion_nr", "random_id")

# Only select files of samples for which both deep seq and WGS is performed
random_ids <- c("P0077", "P0122", "P0144", "P0164", "P0608", "P0180")
target_samples <- random_id_info$skion_nr[match(random_ids, random_id_info$random_id)]
names(target_samples) <- random_ids
vcf_vector <- vcf_vector[sapply(target_samples, function(x) {grep(x, vcf_vector)})]

# Define directory where to store output files
output_dir <- paste0(RESULTS_DIR, "figures")
if (!dir.exists(output_dir)){dir.create(output_dir)}



#---------------------------- Function definitions ---------------------------#

# Define cluster order
tumor_types <- c("ClDx", do.call(paste0, expand.grid("ClR", 1:8)))
cluster_types <- as.character(sapply(tumor_types, function(x) 
  paste0(x, c("", "_A", "_B", "_sub", "_sub_A", "_sub_B", "_sub_C", "_fall", "_fall_A", 
              "_fall_B", "_fall_sub", "_fall_sub_A", "_fall_sub_B"))))
na_clusters <- paste0("NA_", c(paste0(0, 1:9), 10:11))
cluster_order <- c(cluster_types, na_clusters)


splitVcfBasedOnCluster <- function (vcfObject, clusterDataframe) 
{
  #Adapted so "cluster" does not get added to the cluster name
  library(MutationalPatterns)
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
  require(logging)
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
                                            ggplot(data = cluster_specific_data_slice, aes(x = sample, y = AF, group = var)) + 
                                              geom_point(alpha = 0.2) + geom_line(alpha = 0.2) + facet_wrap(~cluster) + 
                                              facet_grid(cols = vars(cluster)) + ylim(0, 1) + 
                                              theme(panel.background = element_rect(fill = "white", colour = "grey50"), 
                                                    text = element_text(size = 16), 
                                                    legend.position = "none", panel.grid.major.y = element_line(colour = "black", size = 0.5), 
                                                    axis.text.x = axis_text_x, axis.title.x = axis_title_x, 
                                                    axis.text.y = axis_text_y, axis.title.y = axis_title_y)})
  names(cluster_allel_frequency_plots) <- unique(alleleFrequencyDataFrame$cluster)
  return(cluster_allel_frequency_plots)
}



#------------------------------ Function calling -----------------------------#

#load("wgs_variant_clusters_210324.Rda")
#load("deepseq_variant_clusters_210324.Rda")
#load("wgs_filtered_vcf_list_210324.Rda")

# Initiate objects to store results
sample_filtered_vcf_list <- rep(list(NA), length(random_ids))
names(sample_filtered_vcf_list) <- random_ids
sample_filtered_variants_clusters_df_list <- rep(list(NA), length(random_ids))
names(sample_filtered_variants_clusters_df_list) <- random_ids
deepseq_muts_df_list <- rep(list(NA), length(random_ids))
names(deepseq_muts_df_list) <- random_ids
qual_stats_df <- data.frame()


## Filter the vcf files and make mutation matrix for each timepoint in each patient
# Obtain the filtered mutations from WGS
# Take mutations per timepoint, get reads/af of these mutation in the deepseq data
# And plot values in histogram
for (vcf_path in vcf_vector) {
  skion_number <- str_extract(vcf_path, "\\d\\d\\d\\d+")
  random_id <- random_id_info[random_id_info$skion_nr == skion_number, 2]
  
  # Read VCF and rename clusters
  sample_vcf <- loadVcfFile(vcf_path,
                            ref_genome,
                            "chromosomes")
  GenomeInfoDb::genome(sample_vcf) <- "hg38"
  
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
  if (random_id == "P0164") {
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
  # Save filtered VCF to object
  sample_hg38_filtered <- sample_vcf_gonl_filtered[rownames(sample_filtered_df)]
  sample_filtered_vcf_list[[random_id]] <- sample_hg38_filtered
  
  # Extract allele frequencies
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
  clusters <- c()
  cluster_names <- c()
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

  # Split clusters where necessary
  if (split == TRUE) {
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
    if (random_id == "P0180") {
      sample_filtered_variants_clusters_df <- splitClusterOnTimePointSpecificAlleleFrequency(cluster_df = sample_filtered_variants_clusters_df,
                                                                                             time_point = "Dx",
                                                                                             cluster = "ClDx",
                                                                                             AF_threshold = 0.25)
    } 
    if (random_id == "P0608") {
      sample_filtered_variants_clusters_df <- splitClusterOnTimePointSpecificAlleleFrequency(cluster_df = sample_filtered_variants_clusters_df,
                                                                                             time_point = "Dx",
                                                                                             cluster = "ClDx",
                                                                                             AF_threshold = 0.15)
    } 
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
  }
  
  # Delete small clusters before renumbering of NA clusters
  cluster_count <- plyr::count(sample_filtered_variants_clusters_df$cluster)
  
  # Correct NA count after deleting the small ones
  # Only relevant if small NA clusters were deleted
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
  
  # Rename filtered WGS output
  sample_filtered_variants_clusters_df <- sample_filtered_variants_clusters_df_ordered

  # Parse deep sequencing data
  deepseq_muts_df <- data.frame()
  ds_muts_VAF10 <- 0
  for (rowname in rownames(sample_filtered_variants_clusters_df)) {
    clustername <- sample_filtered_variants_clusters_df[rowname,]$cluster
    
    # Obtain deep sequencing info of each variant found in WGS data
    coordinates <- strsplit(rowname, split="/")[[1]][1]
    alt <- strsplit(rowname, split = "/")[[1]][2]
    ds_variant <- ds_data[coordinates,]
    
    # Get total reads of the variant
    total_reads_values <- ds_variant[, grepl(paste0("Total_reads_", skion_number), colnames(ds_data)) &
                                       !grepl("F1|F2|F3|X-SK", colnames(ds_data))]

    # Get AF information
    alt_AF_TP_values <- ds_variant[, grepl(paste0(alt, "_AF_", skion_number), colnames(ds_data)) &
                                     !grepl("F1|F2|F3|X-SK", colnames(ds_data))]
    
    # Define new column names
    sample_time_points <- sapply(strsplit(colnames(alt_AF_TP_values), split = "_"), `[`, 3)
    colnames(alt_AF_TP_values) <- paste0(sample_time_points, "_read_based_AF")
    rownames(alt_AF_TP_values) <- paste0(coordinates, "/", alt)
    
    # Order columns by time points and calculate mean VAF
    sample_time_points_ordered <- orderSampleIds(colnames(alt_AF_TP_values))
    alt_AF_TP_values_df <- alt_AF_TP_values[, sample_time_points_ordered]
    alt_AF_TP_values <- as.numeric(alt_AF_TP_values_df)
    alt_AF_mean <- mean(alt_AF_TP_values) 
    
    # Add cluster info
    alt_AF_TP_values_df$cluster <- clustername
    
    # Only use working probes
    if (is.nan(alt_AF_mean) || is.na(alt_AF_mean) || sum(total_reads_values) < 100){
      next
    
      # Add deep sequencing variant to data frame
    } else {
      deepseq_muts_df <- rbind(deepseq_muts_df, 
                               alt_AF_TP_values)
      
      # Calculate quality metrics
      if (max(alt_AF_TP_values) >= 0.1) {
        ds_muts_VAF10 <- ds_muts_VAF10 + 1
      } else {
        next
      } 
    }
  }

  # Save quality stats to object
  qual_stats_df <- rbind(qual_stats_df,
                         data.frame(Sample = random_id,
                                    WGS_variants = nrow(sample_filtered_variants_clusters_df),
                                    DS_total_muts_100DP = nrow(deepseq_muts_df),
                                    DS_muts_100DP_aboveVAF10 = ds_muts_VAF10))

  # Save results to object
  sample_filtered_variants_clusters_df_list[[random_id]] <- sample_filtered_variants_clusters_df
  deepseq_muts_df_list[[random_id]] <- deepseq_muts_df
  
  # Calculate cluster sizes
  cluster_sizes_ds <- table(deepseq_muts_df$cluster)
  cluster_sizes_ds <- paste0("n=", cluster_sizes_ds[unique(deepseq_muts_df$cluster)])
  cluster_sizes_wgs <- table(sample_filtered_variants_clusters_df$cluster)
  cluster_sizes_wgs <- paste0("n=", cluster_sizes_wgs[unique(sample_filtered_variants_clusters_df$cluster)])

  # Generate plots of variants found in deep sequencing and WGS data
  sample_clusterSpecificAlleleFrequencyPlots <- plotClusterSpecificAlleleFrequency(alleleFrequencyDataFrame = deepseq_muts_df,
                                                                                   minimal = FALSE)
  sample_af_merged_plots_ds <- plot_grid(plotlist = sample_clusterSpecificAlleleFrequencyPlots, 
                                         labels = cluster_sizes_ds, ncol = 1, label_size = 10)
  sample_clusterSpecificAlleleFrequencyPlots <- plotClusterSpecificAlleleFrequency(alleleFrequencyDataFrame = sample_filtered_variants_clusters_df,
                                                                                   minimal = FALSE)
  sample_af_merged_plots_wgs <- plot_grid(plotlist = sample_clusterSpecificAlleleFrequencyPlots, 
                                          labels = cluster_sizes_wgs, ncol = 1, label_size = 10)
  
  # Save to PDF
  fn <- paste0(output_dir, "/AFoverview_", random_id, ".pdf")
  pdf(file = fn, width=10, height=20)
  gridExtra::grid.arrange(sample_af_merged_plots_ds, sample_af_merged_plots_wgs, 
                          ncol = 2, top = paste(random_id, "left: DS, right: WGS"))
  dev.off()
}

qual_stats_df



#---------------------------- Overview visualizations -------------------------#

## Generate plots for sample P0077
wgs_variant_clusters_df <- sample_filtered_variants_clusters_df_list[["P0077"]]
deepseq_variant_clusters_df <- deepseq_muts_df_list[["P0077"]]
identical(colnames(wgs_variant_clusters_df), colnames(deepseq_variant_clusters_df))
time_points <- c("Dx", "R1", "R2")
plots <- rep(list(NA), 3)

# Generate plot for each time point
for (i in 1:length(time_points)){
  time_point <- grep(time_points[i], colnames(deepseq_variant_clusters_df))
  plot_title <- c("Diagnosis", "Relapse 1", "Relapse 2")[i]
  ifelse(i == 1,  xcol <- rep(c("black", "red", "black"), each = 2),
                  xcol <- "black")
  
  # WGS plot
  wgs_plot <- ggplot(wgs_variant_clusters_df, 
                     aes(x = cluster, y = wgs_variant_clusters_df[,time_point])) + 
    geom_jitter(size = 0.7) + theme_linedraw() +
    ylim(-0.01, 1.01) +
    labs(title = ifelse(i == 1, "WGS\n", ""), subtitle = plot_title, y = "MAF", x = "") +
    geom_hline(yintercept = 0.01, linetype = "dashed", color = "red") +
    theme(axis.text.x = element_text(colour = xcol),
          panel.grid.major.x = element_blank(),
          panel.grid = element_line(color = "grey80"))
  
  # Deepseq plot
  ds_plot <- ggplot(deepseq_variant_clusters_df, 
                    aes(x = cluster, y = deepseq_variant_clusters_df[,time_point])) + 
    geom_jitter(size = 0.7) + theme_linedraw() +
    ylim(-0.01, 1.01) +
    labs(title = ifelse(i == 1, "DeepSeq\n", ""), subtitle = plot_title, y = "MAF", x = "") +
    geom_hline(yintercept = 0.01, linetype = "dashed", color = "red") +
    theme(axis.text.x = element_text(colour = xcol),
          panel.grid.major.x = element_blank(),
          panel.grid = element_line(color = "grey80"))
  
  combined_plot <- grid.arrange(wgs_plot, ds_plot, ncol = 2)
  plots[[i]] <- combined_plot
}

# Save to PDF
pdf(file = paste0(output_dir, "/dotplots_P", target_samples["P0077"], ".pdf"), 
    width = 11.5, height = 14)
grid.arrange(plots[[1]], plots[[2]], plots[[3]],
             top = textGrob("P0077", gp = gpar(fontsize = 14, font = 2)))
dev.off()




## Generate plots for sample P0122
wgs_variant_clusters_df <- sample_filtered_variants_clusters_df_list[["P0122"]]
deepseq_variant_clusters_df <- deepseq_muts_df_list[["P0122"]]
identical(colnames(wgs_variant_clusters_df), colnames(deepseq_variant_clusters_df))
time_points <- c("Dx", "R1", "R2")
plots <- rep(list(NA), 3)

# Generate plot for each time point
for (i in 1:length(time_points)){
  time_point <- grep(time_points[i], colnames(deepseq_variant_clusters_df))
  plot_title <- c("Diagnosis", "Relapse 1", "Relapse 2")[i]
  ifelse(i == 2,  xcol <- rep(c("black", "black", "red"), each = 2),
                  xcol <- "black")
  
  # WGS plot
  wgs_plot <- ggplot(wgs_variant_clusters_df, 
                     aes(x = cluster, y = wgs_variant_clusters_df[,time_point])) + 
    geom_jitter(size = 0.7) + theme_linedraw() +
    ylim(-0.01, 1.01) +
    labs(title = ifelse(i == 1, "WGS\n", ""), subtitle = plot_title, y = "MAF", x = "") +
    geom_hline(yintercept = 0.01, linetype = "dashed", color = "red") +
    theme(axis.text.x = element_text(colour = xcol),
          panel.grid.major.x = element_blank(),
          panel.grid = element_line(color = "grey80"))

  # Deepseq plot
  ds_plot <- ggplot(deepseq_variant_clusters_df, 
                    aes(x = cluster, y = deepseq_variant_clusters_df[,time_point])) + 
    geom_jitter(size = 0.7) + theme_linedraw() +
    ylim(-0.01, 1.01) +
    labs(title = ifelse(i == 1, "DeepSeq\n", ""), subtitle = plot_title, y = "MAF", x = "") +
    geom_hline(yintercept = 0.01, linetype = "dashed", color = "red") +
    theme(axis.text.x = element_text(colour = xcol),
          panel.grid.major.x = element_blank(),
          panel.grid = element_line(color = "grey80"))
  
  combined_plot <- grid.arrange(wgs_plot, ds_plot, ncol = 2)
  plots[[i]] <- combined_plot
}

# Save to PDF
pdf(file = paste0(output_dir, "/dotplots_P", target_samples["P0122"], ".pdf"), 
    width = 11.5, height = 14)
grid.arrange(plots[[1]], plots[[2]], plots[[3]],
             top = textGrob("P0122", gp = gpar(fontsize = 14, font = 2)))
dev.off()




## Generate plots for sample P0144
wgs_variant_clusters_df <- sample_filtered_variants_clusters_df_list[["P0144"]]
deepseq_variant_clusters_df <- deepseq_muts_df_list[["P0144"]]
time_points <- c("Dx", "R1", "R2", "R3", "R4", "R5")
plots <- rep(list(NA), 6)

# Ensure the same ordering of data frames
wgs_variant_clusters_df <- wgs_variant_clusters_df[,colnames(deepseq_variant_clusters_df)]
identical(colnames(wgs_variant_clusters_df), colnames(deepseq_variant_clusters_df))

# Generate plot for each time point
for (i in 1:length(time_points)){
  time_point <- grep(time_points[i], colnames(deepseq_variant_clusters_df))
  plot_title <- c("Diagnosis", "Relapse 1", "Relapse 2", "Relapse 3", 
                  "Relapse 4", "Relapse 5")[i]
  if (i == 4){
    xcol <- c(rep("black", 4), "red", "red", "black")
  } else {
    xcol <- "black"
  }
  
  # WGS plot
  wgs_plot <- ggplot(wgs_variant_clusters_df, 
                     aes(x = cluster, y = wgs_variant_clusters_df[,time_point])) + 
    geom_jitter(size = 0.7) + theme_linedraw() +
    ylim(-0.01, 1.01) +
    labs(title = ifelse(i == 1, "WGS\n", ""), subtitle = plot_title, y = "MAF", x = "") +
    geom_hline(yintercept = 0.01, linetype = "dashed", color = "red") +
    theme(axis.text.x = element_text(colour = xcol),
          panel.grid.major.x = element_blank(),
          panel.grid = element_line(color = "grey80"))
  
  # Deepseq plot
  ds_plot <- ggplot(deepseq_variant_clusters_df, 
                    aes(x = cluster, y = deepseq_variant_clusters_df[,time_point])) + 
    geom_jitter(size = 0.7) + theme_linedraw() +
    ylim(-0.01, 1.01) +
    labs(title = ifelse(i == 1, "DeepSeq\n", ""), subtitle = plot_title, y = "MAF", x = "") +
    geom_hline(yintercept = 0.01, linetype = "dashed", color = "red") +
    theme(axis.text.x = element_text(colour = xcol),
          panel.grid.major.x = element_blank(),
          panel.grid = element_line(color = "grey80"))
  
  combined_plot <- grid.arrange(wgs_plot, ds_plot, ncol = 2)
  plots[[i]] <- combined_plot
}

# Save to PDF
pdf(file = paste0(output_dir, "/dotplots_P", target_samples["P0144"], ".pdf"), 
    width = 11.5, height = 21)
grid.arrange(plots[[1]], plots[[2]], plots[[3]],
             plots[[4]], plots[[5]], plots[[6]],
             top = textGrob("P0144", gp = gpar(fontsize = 14, font = 2)),
             ncol = 1)
dev.off()




## Generate plots for sample P0164
wgs_variant_clusters_df <- sample_filtered_variants_clusters_df_list[["P0164"]]
deepseq_variant_clusters_df <- deepseq_muts_df_list[["P0164"]]
identical(colnames(wgs_variant_clusters_df), colnames(deepseq_variant_clusters_df))
time_points <- c("Dx", "R1", "R2")
plots <- rep(list(NA), 3)

# Generate plot for each time point
for (i in 1:length(time_points)){
  time_point <- grep(time_points[i], colnames(deepseq_variant_clusters_df))
  plot_title <- c("Diagnosis", "Relapse 1", "Relapse 2")[i]
  y_lim <- c(0.85, 0.67, 0.87)[i]
  
  # WGS plot
  wgs_plot <- ggplot(wgs_variant_clusters_df, 
                     aes(x = cluster, y = wgs_variant_clusters_df[,time_point])) + 
    geom_jitter(size = 0.7) + theme_linedraw() +
    ylim(-0.01, y_lim) + 
    labs(title = ifelse(i == 1, "WGS\n", ""), subtitle = plot_title, y = "MAF", x = "") +
    geom_hline(yintercept = 0.01, linetype = "dashed", color = "red") +
    theme(panel.grid.major.x = element_blank(),
          panel.grid = element_line(color = "grey80"))
  
  # Deepseq plot
  ds_plot <- ggplot(deepseq_variant_clusters_df, 
                    aes(x = cluster, y = deepseq_variant_clusters_df[,time_point])) + 
    geom_jitter(size = 0.7) + theme_linedraw() +
    ylim(-0.01, y_lim) + 
    labs(title = ifelse(i == 1, "DeepSeq\n", ""), subtitle = plot_title, y = "MAF", x = "") +
    geom_hline(yintercept = 0.01, linetype = "dashed", color = "red") +
    theme(panel.grid.major.x = element_blank(),
          panel.grid = element_line(color = "grey80"))
  
  combined_plot <- grid.arrange(wgs_plot, ds_plot, ncol = 2)
  plots[[i]] <- combined_plot
}

# Save to PDF
pdf(file = paste0(output_dir, "/dotplots_P", target_samples["P0164"], ".pdf"), 
    width = 11.5, height = 14)
grid.arrange(plots[[1]], plots[[2]], plots[[3]],
             top = textGrob("P0164", gp = gpar(fontsize = 14, font = 2)))
dev.off()




## Generate plots for sample P0608
wgs_variant_clusters_df <- sample_filtered_variants_clusters_df_list[["P0608"]]
deepseq_variant_clusters_df <- deepseq_muts_df_list[["P0608"]]
identical(colnames(wgs_variant_clusters_df), colnames(deepseq_variant_clusters_df))
time_points <- c("Dx", "R1", "R2", "R3")
plots <- rep(list(NA), 4)

# Generate plot for each time point
for (i in 1:length(time_points)){
  time_point <- grep(time_points[i], colnames(deepseq_variant_clusters_df))
  plot_title <- c("Diagnosis", "Relapse 1", "Relapse 2", "Relapse 3")[i]
  
  # WGS plot
  wgs_plot <- ggplot(wgs_variant_clusters_df, 
                     aes(x = cluster, y = wgs_variant_clusters_df[,time_point])) + 
    geom_jitter(size = 0.7) + theme_linedraw() +
    ylim(-0.01, 1.01) + 
    labs(title = ifelse(i == 1, "WGS\n", ""), subtitle = plot_title, y = "MAF", x = "") +
    geom_hline(yintercept = 0.01, linetype = "dashed", color = "red") +
    theme(panel.grid.major.x = element_blank(),
          panel.grid = element_line(color = "grey80"))
  
  # Deepseq plot
  ds_plot <- ggplot(deepseq_variant_clusters_df, 
                    aes(x = cluster, y = deepseq_variant_clusters_df[,time_point])) + 
    geom_jitter(size = 0.7) + theme_linedraw() +
    ylim(-0.01, 1.01) + 
    labs(title = ifelse(i == 1, "DeepSeq\n", ""), subtitle = plot_title, y = "MAF", x = "") +
    geom_hline(yintercept = 0.01, linetype = "dashed", color = "red") +
    theme(panel.grid.major.x = element_blank(),
          panel.grid = element_line(color = "grey80"))
  
  combined_plot <- grid.arrange(wgs_plot, ds_plot, ncol = 2)
  plots[[i]] <- combined_plot
}

# Save to PDF
pdf(file = paste0(output_dir, "/dotplots_P", target_samples["P0608"], ".pdf"), 
    width = 11.5, height = 16)
grid.arrange(plots[[1]], plots[[2]], 
             plots[[3]], plots[[4]],
             top = textGrob("P0608", gp = gpar(fontsize = 14, font = 2)),
             ncol = 1)
dev.off()




## Generate plots for sample P0180
wgs_variant_clusters_df <- sample_filtered_variants_clusters_df_list[["P0180"]]
deepseq_variant_clusters_df <- deepseq_muts_df_list[["P0180"]]
time_points <- c("Dx", "R1", "R2", "R3")
plots <- rep(list(NA), 4)

# Ensure the same ordering of data frames
wgs_variant_clusters_df <- wgs_variant_clusters_df[,colnames(deepseq_variant_clusters_df)]
identical(colnames(wgs_variant_clusters_df), colnames(deepseq_variant_clusters_df))

# Generate plot for each time point
for (i in 1:length(time_points)){
  time_point <- grep(time_points[i], colnames(deepseq_variant_clusters_df))
  plot_title <- c("Diagnosis", "Relapse 1", "Relapse 2", "Relapse 3")[i]
  if (i == 2){
    xcol <- c("black", "black", rep(c("black", "red"), 2), rep("black", 2))
  } else if (i == 3){
    xcol <- c(rep("black", 5), "red", "black", "red")
  } else {
    xcol <- "black"
  }
  
  # WGS plot
  wgs_plot <- ggplot(wgs_variant_clusters_df, 
                     aes(x = cluster, y = wgs_variant_clusters_df[,time_point])) + 
    geom_jitter(size = 0.7) + theme_linedraw() +
    ylim(-0.01, 1.01) +
    labs(title = ifelse(i == 1, "WGS\n", ""), subtitle = plot_title, y = "MAF", x = "") +
    geom_hline(yintercept = 0.01, linetype = "dashed", color = "red") +
    theme(axis.text.x = element_text(colour = xcol),
          panel.grid.major.x = element_blank(),
          panel.grid = element_line(color = "grey80"))
  
  # Deepseq plot
  ds_plot <- ggplot(deepseq_variant_clusters_df, 
                    aes(x = cluster, y = deepseq_variant_clusters_df[,time_point])) + 
    geom_jitter(size = 0.7) + theme_linedraw() +
    ylim(-0.01, 1.01) +
    labs(title = ifelse(i == 1, "DeepSeq\n", ""), subtitle = plot_title, y = "MAF", x = "") +
    geom_hline(yintercept = 0.01, linetype = "dashed", color = "red") +
    theme(axis.text.x = element_text(colour = xcol),
          panel.grid.major.x = element_blank(),
          panel.grid = element_line(color = "grey80"))
  
  combined_plot <- grid.arrange(wgs_plot, ds_plot, ncol = 2)
  plots[[i]] <- combined_plot
}

# Save to PDF
pdf(file = paste0(output_dir, "/dotplots_P", target_samples["P0180"], ".pdf"), 
    width = 11.5, height = 16)
grid.arrange(plots[[1]], plots[[2]], 
             plots[[3]], plots[[4]],
             top = textGrob("P0180", gp = gpar(fontsize = 14, font = 2)),
             ncol = 1)
dev.off()



#---------------------- Acquired and subclonal mutations ----------------------#

all_timepoints <- c("Dx", paste0("R", 1:5))
names(all_timepoints) <- c("Diagnosis", paste("Relapse", 1:5))

# Prepare data frame with target clusters
target_clusters_df <- rbind(data.frame(Sample = "P0077",
                                       timepoint = "Dx",
                                       clusters = c("ClDx_sub", "ClR1")),
                            data.frame(Sample = "P0122",
                                       timepoint = "R1",
                                       clusters = c("ClR1_sub", "ClR2")),
                            data.frame(Sample = "P0144",
                                       timepoint = "R3",
                                       clusters = c("ClR3_sub", "ClR4")),
                            data.frame(Sample = "P0180",
                                       timepoint = rep(c("R1", "R2"), each = 2),
                                       clusters = c("ClR1", "ClR2", "ClR2", "ClR3")))


# Generate jitter plot of target clusters
sample_plots <- rep(list(NA), 5)
i <- 1
for (random_id in unique(target_clusters_df$Sample)){
  sample_id <- target_samples[random_id]
  
  # Get WGS and DS variants for each sample
  deepseq_variant_clusters_df <- deepseq_muts_df_list[[random_id]]
  wgs_variant_clusters_df <- sample_filtered_variants_clusters_df_list[[random_id]]
  
  # Ensure the same ordering of data frames
  wgs_variant_clusters_df <- wgs_variant_clusters_df[, colnames(deepseq_variant_clusters_df)]
  
  # Get info of each time point where target clusters are present
  for (time_point in unique(target_clusters_df[target_clusters_df$Sample == random_id, 2])){
    next_timepoint <- all_timepoints[grep(time_point, all_timepoints) + 1]
    both_timepoints <- names(all_timepoints)[match(c(time_point, next_timepoint), all_timepoints)]
    
    # Obtain column and cluster information
    target_time_point <- colnames(deepseq_variant_clusters_df)[grep(time_point, colnames(deepseq_variant_clusters_df))]
    target_clusters <- target_clusters_df[target_clusters_df$Sample == random_id & 
                                            target_clusters_df$timepoint == time_point, 3]
    next_timepoint <- gsub(time_point, next_timepoint, target_time_point)
    names(both_timepoints) <- c(target_time_point, next_timepoint)
    
    # Subset WGS and DS variants for target columns and clusters
    target_cols <- c(target_time_point, next_timepoint, "cluster")
    subset_ds_variants <- deepseq_variant_clusters_df[deepseq_variant_clusters_df$cluster %in% target_clusters, 
                                                      colnames(deepseq_variant_clusters_df) %in% target_cols]
    subset_wgs_variants <- wgs_variant_clusters_df[wgs_variant_clusters_df$cluster %in% target_clusters, 
                                                   colnames(wgs_variant_clusters_df) %in% target_cols]
    
    # Add sequencing info to data frames and combine info
    subset_ds_variants["Type"] <- "DeepSeq"
    subset_wgs_variants["Type"] <- "WGS"
    combined_variant_info <- rbind(subset_ds_variants, subset_wgs_variants)
    combined_variant_info$Type <- factor(combined_variant_info$Type, levels = c("WGS", "DeepSeq"))
    
    # Define plotting function 
    #generate_ggplot <- function(df, target_time_point, plot_title, legend_pos = "none"){
    #  ggplot_object <- ggplot(df, 
    #         aes(x = cluster, y = df[,target_time_point], color = Type)) + 
    #    geom_point(size = 0.7, position = position_jitterdodge(dodge.width = 0.85, jitter.width = 0.7)) + 
    #    guides(colour = guide_legend(override.aes = list(size = 2))) +
    #    theme_linedraw() + ylim(-0.01, 1.01) +
    #    scale_color_manual(values = c(DeepSeq = "#0079ce", WGS = "#228b57")) +
    #    labs(subtitle = plot_title, y = "MAF", x = "") +
    #    geom_hline(yintercept = 0.01, linetype = "dashed", color = "red") +
    #    theme(panel.grid.major.x = element_blank(),
    #          panel.grid = element_line(color = "grey80"),
    #          legend.title = element_blank(),
    #          legend.text = element_text(size = 11),
    #          legend.position = legend_pos)
    #  return(ggplot_object)
    #}
    
    # Get legend 
    #sample_ggplot <- generate_ggplot(df = combined_variant_info, target_time_point, plot_title = "",
    #                                legend_pos = "right")
    #legend <- as_ggplot(get_legend(sample_ggplot))
    
    # Generate plot of first timepoint
    #first_ggplot <- generate_ggplot(df = combined_variant_info,
    #                                target_time_point,
    #                                plot_title = as.character(both_timepoints[target_time_point]))
    #second_ggplot <- generate_ggplot(df = combined_variant_info,
    #                                 next_timepoint,
    #                                 plot_title = as.character(both_timepoints[next_timepoint]))
    #sample_plots[[i]] <- grid.arrange(first_ggplot, second_ggplot, ncol = 1, 
    #                                  top = textGrob(random_id, gp = gpar(fontsize = 13, font = 2)))
    
    # Melt and parse data frame
    combined_variant_info_melt <- reshape2::melt(combined_variant_info)
    colnames(combined_variant_info_melt) <- c("cluster", "Type", "timepoint", "AF")
    combined_variant_info_melt$timepoint <- both_timepoints[match(combined_variant_info_melt$timepoint, names(both_timepoints))]
    combined_variant_info_melt$Type <- factor(combined_variant_info_melt$Type, levels = c("WGS", "DeepSeq"))
    
    # Define plotting function 
    generate_ggplot <- function(df, legend_pos = "none"){
      ggplot_object <- ggplot(df, aes(x = cluster, y = AF, color = Type)) + 
        geom_point(size = 0.3, position = position_jitterdodge(dodge.width = 0.85, jitter.width = 0.7)) + 
        guides(colour = guide_legend(fill = FALSE, override.aes = list(size = 2))) +
        theme_linedraw() + ylim(-0.01, 1.01) +
        facet_wrap(. ~ timepoint, ncol = 1) +
        scale_y_continuous(expand = c(0.008, 0.008), limits = c(0, 1)) +
        scale_color_manual(values = c(DeepSeq = "#0079ce", WGS = "#228b57")) +
        labs(title = random_id, y = "MAF", x = "") +
        geom_hline(yintercept = 0.01, linetype = "dashed", color = "red") +
        theme(panel.grid.major.x = element_blank(),
              panel.grid = element_line(color = "grey80"),
              plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
              strip.background = element_rect(fill = "grey90"),
              strip.text = element_text(color = "black", size = 9.5),
              legend.title = element_blank(),
              legend.text = element_text(size = 10),
              legend.position = legend_pos)
      return(ggplot_object)
    }
    # Get legend 
    sample_ggplot <- generate_ggplot(df = combined_variant_info_melt, legend_pos = "right")
    legend <- as_ggplot(get_legend(sample_ggplot)) 
    
    # Generate plot
    sample_ggplot <- generate_ggplot(df = combined_variant_info_melt)
    sample_ggplot <- sample_ggplot + theme(axis.title.y = element_blank())
    
    # Save to object
    sample_plots[[i]] <- sample_ggplot
    i <- i + 1
   }
}
pdf("figures/wrongClustersOverview_newFormat.pdf", width = 13, height = 6)
grid.arrange(sample_plots[[1]], sample_plots[[2]], sample_plots[[3]], 
             sample_plots[[4]], sample_plots[[5]], legend, 
             layout_matrix = matrix(c(rep(c(1,2,3,4,5), each = 5), 6, 6), 
                                    byrow = T, nrow = 1))
dev.off()

#grid.arrange(sample_plots[[1]], sample_plots[[2]], 
#             sample_plots[[3]], sample_plots[[4]], 
#             sample_plots[[5]], legend, 
#             layout_matrix = matrix(c(rep(c(rep(1, 11), NA, 
#                                            rep(2, 11), NA, 
#                                            rep(3, 11), NA, 
#                                            rep(4, 11), NA, 
#                                            rep(5, 11)), 20),
#                                      c(rep(NA, 29), 6, rep(NA, 29))), 
#                                    byrow = T, nrow = 21))
#grid.arrange(sample_plots[[1]], sample_plots[[2]], sample_plots[[3]], 
#             sample_plots[[4]], sample_plots[[5]], legend, 
#             layout_matrix = matrix(c(rep(c(1,2,3,4,5), each = 3), 6), 
#                                    byrow = T, nrow = 1))



## Determine acquired and subclonal variants per target cluster in each sample
variant_info_df <- data.frame()
wgs_variant_status_info <- list()
for (random_id in unique(target_clusters_df$Sample)){
  sample_id <- target_samples[random_id]
  wgs_variant_info <- c()
  
  # Get WGS and DS variants for each sample
  deepseq_variant_clusters_df <- deepseq_muts_df_list[[random_id]]
  wgs_variant_clusters_df <- sample_filtered_variants_clusters_df_list[[random_id]]
  
  # Ensure the same ordering of data frames
  wgs_variant_clusters_df <- wgs_variant_clusters_df[, colnames(deepseq_variant_clusters_df)]
  
  # Get info of each time point where target clusters are present
  for (time_point in unique(target_clusters_df[target_clusters_df$Sample == random_id, 2])){
    target_time_point <- colnames(deepseq_variant_clusters_df)[grep(time_point, colnames(deepseq_variant_clusters_df))]
    target_clusters <- target_clusters_df[target_clusters_df$Sample == random_id & 
                                            target_clusters_df$timepoint == time_point, 3]
    
    # Define which cluster is acquired and which is subclonal
    names(target_clusters) <- ifelse(grepl("sub", target_clusters) | 
                                       random_id == "P0180" & target_clusters == "ClR1" |
                                       random_id == "P0180" & target_clusters == "ClR2" & time_point == "R2",
                                     "Subclonal", "Acquired")
    
    # Loop over each target cluster within the target time point
    for (target_cluster in target_clusters){
      status <- names(target_clusters)[target_clusters == target_cluster]
      acquired <- 0
      subclonal <- 0
      
      # Get WGS variants belonging to the target cluster and time point 
      subset_wgs_muts_df <- wgs_variant_clusters_df[wgs_variant_clusters_df$cluster == target_cluster, 
                                                    target_time_point, drop = FALSE]
      wgs_variants <- rownames(subset_wgs_muts_df)
      
      # Get the DS allele frequency for each variant found in the WGS data
      for (variant in wgs_variants){
        deepseq_AF <- deepseq_variant_clusters_df[variant, target_time_point]
        
        # Check if the variant is acquired or subclonal 
        if ( !is.nan(deepseq_AF) && !is.na(deepseq_AF)) {
          if (deepseq_AF <= noise_threshold_deep) {
            acquired <- acquired + 1
            wgs_variant_info <- c(wgs_variant_info, paste(variant, target_cluster, target_clusters["Acquired"], sep = "-"))
            
          } else if (deepseq_AF > noise_threshold_deep) {
            subclonal <- subclonal + 1
            wgs_variant_info <- c(wgs_variant_info, paste(variant, target_cluster, target_clusters["Subclonal"], sep = "-"))
          }
        } else {
          next
        }
      }
      
      # Calculate the validation_rate
      if (status == "Subclonal"){
        validation_rate <- subclonal / (acquired + subclonal)
      } else {
        validation_rate <- acquired / (acquired + subclonal)
      }
      
      # Add all info to data frame
      cluster_info_df <- data.frame(Sample = sample_id,
                                    Random_id = random_id,
                                    Time_point = time_point,
                                    Target_cluster = target_cluster,
                                    Status = status,
                                    Total_variants = (acquired + subclonal),
                                    Acquired = acquired,
                                    Subclonal = subclonal,
                                    Validation_rate = validation_rate)
      variant_info_df <- rbind(variant_info_df, cluster_info_df)
    }
  }
  # Add WGS variant status to data frame
  wgs_variant_status_info[[random_id]] <- wgs_variant_info
}
variant_info_df <- variant_info_df[order(variant_info_df$Status), ]
variant_info_df



## Correct and visualize acquired and subclonal variants per target cluster 
# Generate mutation matrices of original and corrected clusters
overview_mut_mat <- matrix(nrow = 96)
for (random_id in unique(target_clusters_df$Sample)){
  
  # Get VCF info
  sample_vcf <- sample_filtered_vcf_list[[random_id]]

  # Extract variant info
  variant_info <- wgs_variant_status_info[[random_id]]
  original_clusters <- sapply(strsplit(variant_info, "-"), `[`, 2)
  corrected_clusters <- sapply(strsplit(variant_info, "-"), `[`, 3)
  variant_info <- sapply(strsplit(variant_info, "-"), `[`, 1)
  
  # Create mutation matrix of original and corrected clusters
  for (cluster in unique(original_clusters)){
    
    # Get VCF entries of original clusters
    orig_vars <- sample_vcf[unique(variant_info[original_clusters == cluster]), ]
    
    # Get VCF entries of corrected clusters
    if (cluster == "ClR2" & random_id == "P0180"){
      corr_1 <- unique(variant_info[corrected_clusters == "ClR1"])
      corr_2 <- unique(variant_info[corrected_clusters == "ClR2"])
      corr_3 <- unique(variant_info[corrected_clusters == "ClR3"])
      correct_vars_2 <- corr_2[!corr_2 %in% c(intersect(corr_1, corr_2),
                                              intersect(corr_3, corr_2))]
      correct_vars <- sample_vcf[correct_vars_2, ]
    } else {
      correct_vars <- sample_vcf[unique(variant_info[corrected_clusters == cluster]), ]
    }
    
    # Convert to GRanges and create mutation matrix
    cluster_mut_mat <- cbind(mut_matrix(granges(orig_vars), ref_genome = ref_genome),
                             mut_matrix(granges(correct_vars), ref_genome = ref_genome))
    colnames(cluster_mut_mat) <- c(paste0(random_id, "-Original ", cluster),
                                   paste0(random_id, "-Corrected ", cluster))
    
    # Save to object
    overview_mut_mat <- cbind(overview_mut_mat, cluster_mut_mat)
  }
}
overview_mut_mat <- overview_mut_mat[, colnames(overview_mut_mat) != ""]

# Generate mutation profiles
i <- 1
overview_plot <- list()
for (random_id in unique(target_clusters_df$Sample)){
  sample_mut_mat <- overview_mut_mat[, grepl(random_id, colnames(overview_mut_mat))]
  colnames(sample_mut_mat) <- paste0(gsub(paste0(random_id, "-"), "", colnames(sample_mut_mat)),
                                     "\n(n = ", colSums(sample_mut_mat), ")")

  # Function to generate profile
  generate_profile <- function(mut_mat, ymax, title, show_x) {
    status <- strsplit(colnames(mut_mat), " ")[[1]][1]
    colnames(mut_mat) <- gsub(paste0(status, " "), "", colnames(mut_mat))
    profile <- plot_96_profile(mut_mat,
                               ymax = ymax,
                               condensed = F) + 
      labs(title = title, x = "Sequence context") +
      theme(axis.text = element_text(size = 10, color = "black"),
            axis.text.x = element_text(color = c("black", "transparent")),
            plot.title = element_text(hjust = 0.5, size = 12),
            axis.title = element_text(size = 10),
            axis.title.y = element_blank(),
            strip.text = element_text(size = 10.5))
    if (isFALSE(show_x)){
      profile <- profile + theme(axis.text.x = element_blank(),
                                 axis.title.x = element_blank(),
                                 axis.ticks.x = element_blank())
    }
    return(profile)
  }

  # Generate mutation_profile of original clusters
  ymax <- c(0.1, 0.28, 0.1, 0.21)[i]
  original_profile <- generate_profile(mut_mat = sample_mut_mat[,grepl("Original", colnames(sample_mut_mat))],
                                       title = ifelse(i == 1, "Original profile", ""), 
                                       show_x = ifelse(i == 4, TRUE, FALSE), ymax = ymax)
  corrected_profile <- generate_profile(mut_mat = sample_mut_mat[,grepl("Corrected", colnames(sample_mut_mat))],
                                        title = ifelse(i == 1, "Corrected profile", ""), 
                                        show_x = ifelse(i == 4, TRUE, FALSE), ymax = ymax)
  if (i == 1){
    overview_plot[[random_id]] <- grid.arrange(original_profile, corrected_profile, ncol = 2,
                                  top = textGrob(random_id, gp = gpar(fontsize = 12, font = 2)))
  } else {
    overview_plot[[random_id]] <- grid.arrange(original_profile + theme(strip.text.x.top = element_blank(), plot.title = element_blank()), 
                                   corrected_profile + theme(strip.text.x.top = element_blank(), plot.title = element_blank()), 
                                   ncol = 2,
                                   top = textGrob(random_id, gp = gpar(fontsize = 12, font = 2)))
  }
  i <- i + 1
}

# Combine figures and save to PDF
pdf("figures/OriginalAndCorrectedMutProfiles.pdf", width = 9, height = 10)
cowplot::plot_grid(overview_plot[[1]], overview_plot[[2]], overview_plot[[3]], overview_plot[[4]],
             ncol = 1, rel_heights = c(1.05,0.8,0.8,1.3))
dev.off()

# Visualise original vs reconstructed cosine similarities
orig_mut_mat <- overview_mut_mat[,grepl("Original", colnames(overview_mut_mat))]
corr_mut_mat <-  overview_mut_mat[,grepl("Corrected", colnames(overview_mut_mat))]
colnames(orig_mut_mat) <- gsub("-Original", "", colnames(orig_mut_mat))
pdf("figures/origVsReconstructedCosSim_correctedProfiles.pdf", width = 6, height = 4)
plot_original_vs_reconstructed(orig_mut_mat[,rev(colnames(orig_mut_mat))],
                               corr_mut_mat[,rev(colnames(corr_mut_mat))], y_intercept = 0.85) + 
  coord_flip() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title = element_text(size = 10))
dev.off()
pdf("figures/origVsReconstructedCosSim_correctedProfiles_rev.pdf", width = 7, height = 5)
plot_original_vs_reconstructed(orig_mut_mat,
                               corr_mut_mat,
                               y_intercept = 0.85,
                               ylims = c(0,1))
dev.off()

# Print cosine similarities
data.frame(Sample = colnames(orig_mut_mat), 
           Cossim = diag(cos_sim_matrix(orig_mut_mat, corr_mut_mat)))



#---------------------------------- Statistics --------------------------------#

## Calculate validation rate for each time point in each cluster
stats_df <- data.frame()
for (i in (1:length(sample_filtered_variants_clusters_df_list))){
  
  # Get WGS and DS variants for each sample
  deepseq_muts_df <- deepseq_muts_df_list[[i]]
  sample_filtered_variants_clusters_df <- sample_filtered_variants_clusters_df_list[[i]]

  # Determine correctly and wrongly assigned variants in each cluster at each time point
  for (clustername in unique(sample_filtered_variants_clusters_df$cluster)) {
    for (timepoint in colnames(sample_filtered_variants_clusters_df)[1:(length(colnames(sample_filtered_variants_clusters_df))-1)]) {
      correct <- 0
      wrong <- 0
      
      for (rowname in rownames(sample_filtered_variants_clusters_df[sample_filtered_variants_clusters_df$cluster == clustername,])) {
        WGS_AF <- sample_filtered_variants_clusters_df[rowname,timepoint]
        Deepseq_AF <- deepseq_muts_df[rowname,timepoint]

        # Determine if WGS variants are assigned correctly based on DS info
        if (is.nan(Deepseq_AF) || is.na(Deepseq_AF)) {
          next
        } else {
          if (WGS_AF == 0 && Deepseq_AF <= noise_threshold_deep) {
            correct <- correct + 1
          } else if (WGS_AF == 0 && Deepseq_AF > noise_threshold_deep) {
            wrong <- wrong + 1
          } else if (WGS_AF > 0 && Deepseq_AF > noise_threshold_deep) {
            correct <- correct + 1
          } else if (WGS_AF > 0 && Deepseq_AF <= noise_threshold_deep) {
            wrong <- wrong + 1
          } 
        }
      }

      # Calculate stats and save to object
      total <- wrong + correct
      validationrate <- correct / total
      sample_stats_df <- data.frame(RandomID = names(sample_filtered_variants_clusters_df_list)[i],
                                    Cluster = clustername,
                                    TimePoint = timepoint,
                                    Total = total,
                                    Correct = correct,
                                    ValidationRate = validationrate)
      stats_df <- rbind(stats_df, sample_stats_df)
    }
  }
}
stats_df

# Calculate total validation rates per sample
detach(package:plyr)
library(dplyr)
compl_stats_df <- stats_df %>% 
                    group_by(RandomID) %>% 
                    summarise(Total_vars = sum(Total),
                              Correct_vars = sum(Correct))
compl_stats_df <- as.data.frame(compl_stats_df)
compl_stats_df["ValidationRate"] <- compl_stats_df$Correct_vars / compl_stats_df$Total_vars
compl_stats_df["Sample"] <- target_samples[match(compl_stats_df$RandomID, names(target_samples))]
compl_stats_df


## Calculate mean standard deviation of variability in DS and WGS data
stdev_df <- data.frame()
for (sample in names(deepseq_muts_df_list)){

  # Get WGS and DS variants for each sample
  deepseq_df <- deepseq_muts_df_list[[sample]]
  wgs_df <- sample_filtered_variants_clusters_df_list[[sample]]
  n_tp <- ncol(wgs_df) - 1

  # Get WGS and DS variants for each cluster
  for (clustername in unique(wgs_df$cluster)) {
    cluster_deepseq_df <- deepseq_df[deepseq_df$cluster == clustername, 1:n_tp]
    cluster_wgs_df <- wgs_df[wgs_df$cluster == clustername, 1:n_tp]
    
    # Reorder WGS info to match order of DS info
    cluster_wgs_df <- cluster_wgs_df[, names(cluster_deepseq_df)]
    
    # Calculate stdev for each time point
    deepseq_stdev <- apply(cluster_deepseq_df, 2, sd)
    wgs_stdev <- apply(cluster_wgs_df, 2, sd)
    
    # Create data frame of stdev info
    cluster_stdev_info <- data.frame(Sample = sample,
                                     Timepoint = names(deepseq_stdev),
                                     Cluster = clustername,
                                     DeepSeq_stdev = deepseq_stdev,
                                     WGS_stdev = wgs_stdev, row.names = NULL)
    
    # Add info to object
    stdev_df <- rbind(stdev_df, cluster_stdev_info)
  }
}
mean(stdev_df$DeepSeq_stdev)
mean(stdev_df$WGS_stdev)
#openxlsx::write.xlsx(stdev_df, file = "DeepSeq_WGS_StdevInfo.xlsx")

# Check mean stdev for all clusters that contain a mutation
stdev_df_subset <- stdev_df[stdev_df$WGS_stdev > 0, ]
mean(stdev_df_subset$DeepSeq_stdev)
mean(stdev_df_subset$WGS_stdev)



#------------------------------- Miscellaneous ------------------------------#

library("reshape2")
library("plyr")

## Visualizing the frequency of VAFs found in the DS data
freq_df <- data.frame()
for (sample in names(deepseq_muts_df_list)){
  deepseq_df <- deepseq_muts_df_list[[sample]]
  sample_VAFs <- melt(deepseq_df)[,3]
  sample_VAFs <- round_any(sample_VAFs, 
                           accuracy = 0.01, 
                           f = ceiling)
  sample_VAFs_freq <- data.frame(VAF = as.numeric(names(table(sample_VAFs))),
                                 Frequency = as.numeric(table(sample_VAFs)),
                                 Sample = sample)
  freq_df <- rbind(freq_df,
                   sample_VAFs_freq[sample_VAFs_freq$VAF != 0, ])
}

# Generate plot and save to PDF
pdf("deepSequencing_VAFoverview.pdf", width = 15, height = 10)
ggplot(freq_df, aes(x = VAF, y = Frequency)) +
  facet_grid(Sample ~ ., scales = "free") +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_continuous(breaks = seq(0.01,1,0.01), expand = c(0.001, 0.001)) +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_line(color = "grey95"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = -90))
dev.off()


## Visualizing VAFs of variants not corresponding to each sample
all_targets <- unlist(sapply(1:length(target_samples), function(x) {
                      rownames(deepseq_muts_df_list[[x]])}))
freq_df_2 <- data.frame()
for (sample in names(deepseq_muts_df_list)){
  sample_id <- target_samples[sample]
  sample_variants <- rownames(deepseq_muts_df_list[[sample]])
  nonsample_variants <- sapply(strsplit(setdiff(all_targets, sample_variants), "/"), `[`, 1)
  nonsample_variants_ds_info <- ds_data[nonsample_variants,
                                        grepl(paste0("[ACGT]_AF_", sample_id), colnames(ds_data)) &
                                        !grepl("F1|F2|F3|X-SK", colnames(ds_data))]

  # Obtain VAFs of non-reference nucleotides
  sample_VAFs <- c()
  for (i in (1:nrow(nonsample_variants_ds_info))){
    row <- nonsample_variants_ds_info[i, ]
    ref_allele <- strsplit(rownames(row), "_")[[1]][2]
    row_VAFs <- row[, !grepl(paste0(ref_allele, "_"), colnames(row))]
    sample_VAFs <- c(sample_VAFs, as.numeric(row_VAFs))
  }
  
  # Determine frequencies
  sample_VAFs <- round_any(sample_VAFs, 
                           accuracy = 0.01, 
                           f = ceiling)
  sample_VAFs_freq <- data.frame(VAF = as.numeric(names(table(sample_VAFs))),
                                 Frequency = as.numeric(table(sample_VAFs)),
                                 Sample = sample)
  freq_df_2 <- rbind(freq_df_2, sample_VAFs_freq)
}

# Generate plot and save to PDF
pdf("deepSequencing_VAFoverview_non-assignedVariantsPerSample_0removed.pdf", width = 8, height = 7)
#test <- freq_df_2[freq_df_2$VAF != 0,]
ggplot(freq_df_2, aes(x = VAF, y = Frequency)) +
  facet_grid(Sample ~ ., scales = "free") +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_continuous(breaks = seq(0, 1, 0.01), expand = c(0.001, 0.001),
                     limits = c(-0.01, 0.4)) +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_line(color = "grey95"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = -90))
dev.off()


## Determine frequency of variants with a VAF <= 0.01
vaf_stats_df <- data.frame()
for (sample in random_ids){
  sample_freq <- freq_df_2[freq_df_2$Sample == sample, ]
  total_freq <- sum(sample_freq$Frequency)
  
  below_cutoff_freq <- sum(sample_freq[sample_freq$VAF <= 0.01, 2])
  around_cutoff_freq <- as.numeric(sample_freq[sample_freq$VAF == 0.02, 2])
  above_cutoff_freq <- sum(sample_freq[sample_freq$VAF > 0.02, 2])
  
  validation_rate <- round((below_cutoff_freq * 100 / total_freq), digits = 3)
  around_cutoff_rate <- round((around_cutoff_freq * 100 / total_freq), digits = 3)
  error_rate <- round((above_cutoff_freq * 100 / total_freq), digits = 3)
  cat(paste0("Sample ", sample, ".\n", validation_rate, "% variants found <= 0.01.\n", 
             around_cutoff_rate, "% variants found between 0.01 and 0.02.\n",
             error_rate, "% variants found >0.02.\n\n"))
  
  vaf_stats_df <- rbind(vaf_stats_df,
                        data.frame(Below = below_cutoff_freq * 100 / total_freq,
                                   Around = around_cutoff_freq * 100 / total_freq,
                                   Above = above_cutoff_freq * 100 / total_freq))
}
colMeans(vaf_stats_df)

