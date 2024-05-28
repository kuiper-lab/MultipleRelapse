#!/usr/bin/env Rscript 

## Load VCF file
loadVcfFile <- function (vcfFilePath, bsgGenomeName, chromosomeSet) 
{
  require(VariantAnnotation)
  require(logging)
  CHROMOSOMES_TO_LOAD <- loadPredefinedChromosomeSet(bsgGenomeName, 
                                                     chromosomeSet)
  referenceBsGenome <- loadReferenceBsGenome(bsgGenomeName)
  ref_style <- seqlevelsStyle(referenceBsGenome)
  ref_organism <- GenomeInfoDb::organism(referenceBsGenome)
  vcf <- VariantAnnotation::readVcf(vcfFilePath, bsgGenomeName)
  loginfo(paste0("loaded ", dim(vcf)[[1]], " variants from file: ", 
                 vcfFilePath))
  loginfo(paste0("changing style from: ", paste(seqlevelsStyle(vcf), 
                                                collapse = ", "), " to: ", ref_style, " for vcf-file: ", 
                 vcfFilePath))
  seqlevelsStyle(vcf) <- ref_style
  loginfo(paste0("keep seq levels ", paste(CHROMOSOMES_TO_LOAD, 
                                           collapse = ", "), " for vcf-file: ", vcfFilePath))
  chromosomes_present_and_to_keep <- CHROMOSOMES_TO_LOAD[CHROMOSOMES_TO_LOAD %in% 
                                                           levels(unique(seqnames(rowRanges(vcf))))]
  vcf <- keepSeqlevels(x = vcf, value = chromosomes_present_and_to_keep, 
                       pruning.mode = "coarse")
  return(vcf)
}


## Load chromosome set
loadPredefinedChromosomeSet <- function (referenceGenomeBuild = c("BSgenome.Hsapiens.UCSC.hg19", 
                                                                  "BSgenome.Hsapiens.NCBI.GRCh38", "BSgenome.Hsapiens.UCSC.hg38"), 
                                         chromosomeSet = c("chromosomes", "autosomes", "autosomes_and_x")) 
{
  require(logging)
  referenceGenomeBuild <- match.arg(referenceGenomeBuild)
  chromosomeSet <- match.arg(chromosomeSet)
  HG_ALL_CHROMOSOMES <- c("chr1", "chr2", "chr3", "chr4", "chr5", 
                          "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", 
                          "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", 
                          "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
  HG_AUTOSOMES <- c("chr1", "chr2", "chr3", "chr4", "chr5", 
                    "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", 
                    "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", 
                    "chr19", "chr20", "chr21", "chr22")
  HG_AUTOSOMES_AND_X <- c("chr1", "chr2", "chr3", "chr4", "chr5", 
                          "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", 
                          "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", 
                          "chr19", "chr20", "chr21", "chr22", "chrX")
  B_ALL_CHROMOSOMES <- c("1", "2", "3", "4", "5", "6", "7", 
                         "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", 
                         "18", "19", "20", "21", "22", "X", "Y")
  B_AUTOSOMES <- c("1", "2", "3", "4", "5", "6", "7", "8", 
                   "9", "10", "11", "12", "13", "14", "15", "16", "17", 
                   "18", "19", "20", "21", "22")
  B_AUTOSOMES_AND_X <- c("1", "2", "3", "4", "5", "6", "7", 
                         "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", 
                         "18", "19", "20", "21", "22", "X")
  predefinedChromosomeSet <- NULL
  if ((referenceGenomeBuild == "BSgenome.Hsapiens.NCBI.GRCh38") && 
      (chromosomeSet == "chromosomes")) {
    predefinedChromosomeSet <- B_ALL_CHROMOSOMES
  }
  if ((referenceGenomeBuild == "BSgenome.Hsapiens.NCBI.GRCh38") && 
      (chromosomeSet == "autosomes")) {
    predefinedChromosomeSet <- B_AUTOSOMES
  }
  if ((referenceGenomeBuild == "BSgenome.Hsapiens.NCBI.GRCh38") && 
      (chromosomeSet == "autosomes_and_x")) {
    predefinedChromosomeSet <- B_AUTOSOMES_AND_X
  }
  if ((referenceGenomeBuild == "BSgenome.Hsapiens.UCSC.hg19") && 
      (chromosomeSet == "chromosomes")) {
    predefinedChromosomeSet <- HG_ALL_CHROMOSOMES
  }
  if ((referenceGenomeBuild == "BSgenome.Hsapiens.UCSC.hg19") && 
      (chromosomeSet == "autosomes")) {
    predefinedChromosomeSet <- HG_AUTOSOMES
  }
  if ((referenceGenomeBuild == "BSgenome.Hsapiens.UCSC.hg19") && 
      (chromosomeSet == "autosomes_and_x")) {
    predefinedChromosomeSet <- HG_AUTOSOMES_AND_X
  }
  if ((referenceGenomeBuild == "BSgenome.Hsapiens.UCSC.hg38") && 
      (chromosomeSet == "chromosomes")) {
    predefinedChromosomeSet <- HG_ALL_CHROMOSOMES
  }
  if ((referenceGenomeBuild == "BSgenome.Hsapiens.UCSC.hg38") && 
      (chromosomeSet == "autosomes")) {
    predefinedChromosomeSet <- HG_AUTOSOMES
  }
  if ((referenceGenomeBuild == "BSgenome.Hsapiens.UCSC.hg38") && 
      (chromosomeSet == "autosomes_and_x")) {
    predefinedChromosomeSet <- HG_AUTOSOMES_AND_X
  }
  loginfo(paste0("loaded: ", chromosomeSet, " for genome build: ", 
                 referenceGenomeBuild))
  return(predefinedChromosomeSet)
}


## Load reference genome
loadReferenceBsGenome <- function (bsGenomeReferenceGenomeName) 
{
  require(BSgenome)
  require(logging)
  library(bsGenomeReferenceGenomeName, character.only = TRUE)
  ref_genome <- base::get(bsGenomeReferenceGenomeName)
  ref_organism <- GenomeInfoDb::organism(ref_genome)
  ref_style <- seqlevelsStyle(ref_genome)
  genome_name <- genome(ref_genome)[[1]]
  loginfo(paste0("Loading genome: ", bsGenomeReferenceGenomeName, 
                 " (organism: ", ref_organism, " with reference style: ", 
                 ref_style, ")"))
  return(ref_genome)
}

## Convert VCF to data frame
convertVcfObjectToDataFrame <- function (vcfObject, addReadDepthInfo = TRUE, removeMultiAllelicSites = FALSE) 
{
  require(ensemblVEP)
  if (removeMultiAllelicSites) {
    vcfObject <- removeMultiAllelicSitesFromVcf(vcfObject)
  }
  rowRange <- extractRowRangesAsDataFrame(vcfObject)
  vep <- as.data.frame(parseCSQToGRanges(vcfObject))
  if (addReadDepthInfo) {
    read_depth <- as.data.frame(extractAllelFrequenciesFromVcf(vcfObject))
    merged_df <- cbind(rowRange, vep[, 7:ncol(vep)], read_depth)
  }
  else {
    merged_df <- cbind(rowRange, vep[, 7:ncol(vep)])
  }
  row.names(merged_df) <- row.names(vcfObject)
  return(merged_df)
}


## Filter VCF on centromeric variants
excludeVariantsInCentromericRegions <- function (vcfObject, bed_file_containing_centromeric_regions, 
                                                 centromeric_variant_mode = FALSE) 
{
  require(logging)
  require(bedr)
  centromeres_bed_file <- read.table(file = bed_file_containing_centromeric_regions, 
                                     header = FALSE, sep = "\t")
  centromeres_bed_file <- centromeres_bed_file[, 1:3]
  colnames(centromeres_bed_file) <- c("chr", "start", "end")
  centromeres_bed_file$chr <- as.character(centromeres_bed_file$chr)
  vcfObject_as_bed <- extractRowRangesAsDataFrame(vcfObject)
  vcfObject_as_bed$chr <- as.character(vcfObject_as_bed$chr)
  vcfObject_as_bed$end <- vcfObject_as_bed$end + 1
  variants_outside_centromere_idx <- !in.region(vcfObject_as_bed, 
                                                centromeres_bed_file, check.sort = FALSE)
  variants_inside_centromere_idx <- !variants_outside_centromere_idx
  loginfo(paste0("Filtering on centromeric regions of file: ", 
                 bed_file_containing_centromeric_regions))
  loginfo(paste0("Started with: ", dim(vcfObject)[1], " variants"))
  loginfo(paste0("Kept: ", sum(variants_outside_centromere_idx), 
                 " variants"))
  loginfo(paste0("Removed: ", sum(!variants_outside_centromere_idx), 
                 " variants"))
  if (centromeric_variant_mode) {
    return(vcfObject[variants_inside_centromere_idx])
  }
  else {
    return(vcfObject[variants_outside_centromere_idx])
  }
}


## Filter VCF for somatic variants
filterOnReadsInNormalSamples <- function (vcfObject, allowedNumberOfReads = 0) 
{
  require(logging)
  af_df <- extractAllelFrequenciesFromVcf(vcfObject)
  read_count_alt_idx <- endsWith(colnames(af_df), "_read_count_alt")
  read_count_alt_df <- af_df[, read_count_alt_idx]
  read_count_alt_remission_idx <- getRemissionIds(colnames(read_count_alt_df))
  read_count_alt_remission_df <- read_count_alt_df[, read_count_alt_remission_idx]
  read_count_alt_remission_df_boolean <- read_count_alt_remission_df >= 
    (allowedNumberOfReads + 1)
  read_count_alt_remission_df_boolean <- as.data.frame(read_count_alt_remission_df_boolean)
  filter <- as.vector(apply(read_count_alt_remission_df_boolean, 
                            1, any))
  read_count_alt_remission_df_boolean$filter <- filter
  variants_to_remove <- rownames(read_count_alt_remission_df_boolean[read_count_alt_remission_df_boolean$filter, 
  ])
  filtered_vcfObject <- vcfObject[!(rownames(vcfObject) %in% 
                                      variants_to_remove), ]
  loginfo(paste0("Filtering on alternative read counts in the normal samples"))
  loginfo(paste0("Used sample ids: ", paste(read_count_alt_remission_idx, 
                                            collapse = ", ")))
  loginfo(paste0("Started with: ", nrow(vcfObject), " variants"))
  loginfo(paste0("Kept: ", nrow(filtered_vcfObject), " variants"))
  loginfo(paste0("Removed: ", (nrow(vcfObject) - nrow(filtered_vcfObject)), 
                 " variants"))
  return(filtered_vcfObject)
}


## Filter VCF on GnomAD allele frequencies
filterOnGnomadAlleleFrequency <- function (vcfObject, gnomadGenomesAllelFrequencyThreshold, useNA = TRUE) 
{
  require(logging)
  vcfObject_as_df <- convertVcfObjectToDataFrame(vcfObject, 
                                                 addReadDepthInfo = FALSE)
  gnomADg_AF_idx <- which(colnames(vcfObject_as_df) == "gnomADg_AF")
  values_na <- is.na(vcfObject_as_df[[gnomADg_AF_idx]])
  values_smaller_than_threshold <- as.numeric(vcfObject_as_df[[gnomADg_AF_idx]]) <= 
    gnomadGenomesAllelFrequencyThreshold
  values_smaller_than_threshold <- replace(values_smaller_than_threshold, 
                                           is.na(values_smaller_than_threshold), FALSE)
  merged_indeces <- Reduce("|", list(c(values_na), c(values_smaller_than_threshold)))
  rows_to_keep <- rownames(vcfObject_as_df)[merged_indeces]
  filtered_vcf <- vcfObject[rownames(vcfObject) %in% rows_to_keep, 
  ]
  loginfo(paste0("Filtering on GnomAD genomes allele frequencies"))
  loginfo(paste0("Started with: ", dim(vcfObject)[1], " variants"))
  loginfo(paste0("Kept: ", nrow(filtered_vcf), " variants"))
  loginfo(paste0("Removed: ", (dim(vcfObject)[1] - nrow(filtered_vcf)), 
                 " variants"))
  return(filtered_vcf)
}


## Filter VCF on GoNL allele frequencies
filterOnGonlAlleleFrequency <- function (vcfObject, gonlAllelFrequencyThreshold) 
{
  library(TVTB)
  require(logging)
  vcfObject_as_df <- convertVcfObjectToDataFrame(vcfObject, 
                                                 addReadDepthInfo = FALSE)
  goNL_AF_idx <- which(colnames(vcfObject_as_df) == "goNL_AF")
  values_na <- is.na(vcfObject_as_df[[goNL_AF_idx]])
  values_smaller_than_threshold <- as.numeric(vcfObject_as_df[[goNL_AF_idx]]) <= 
    gonlAllelFrequencyThreshold
  values_smaller_than_threshold <- replace(values_smaller_than_threshold, 
                                           is.na(values_smaller_than_threshold), FALSE)
  merged_indeces <- Reduce("|", list(c(values_na), c(values_smaller_than_threshold)))
  rows_to_keep <- rownames(vcfObject_as_df)[merged_indeces]
  filtered_vcf <- vcfObject[rownames(vcfObject) %in% rows_to_keep, 
  ]
  loginfo(paste0("Filtering on GoNL genomes allele frequencies"))
  loginfo(paste0("Started with: ", dim(vcfObject)[1], " variants"))
  loginfo(paste0("Kept: ", nrow(filtered_vcf), " variants"))
  loginfo(paste0("Removed: ", (dim(vcfObject)[1] - nrow(filtered_vcf)), 
                 " variants"))
  return(filtered_vcf)
}


## Order sample based on time point
orderSampleIds <- function (sampleIDS) 
{
  require(logging)
  d_sample <- c()
  f_sample <- c()
  relapse_fractions <- c()
  for (sample in sampleIDS) {
    if (grepl("D", sample, )) {
      d_sample <- sample
      loginfo(paste0("found the diagnosis: ", sample, "\n"))
    }
    else if (grepl("_F", sample)) {
      f_sample <- c(f_sample, sample)
      loginfo(paste0("found the full remission: ", sample, 
                     "\n"))
    }
    else if (grepl("CR[[:digit:]]", sample)) {
      f_sample <- c(f_sample, sample)
      loginfo(paste0("found the full remission: ", sample, 
                     "\n"))
    }
    else if (grepl("R[[:digit:]]", sample)) {
      relapse_fractions <- c(relapse_fractions, sample)
    }
  }
  if (length(relapse_fractions) > 0) {
    ordered_relapse_fractions <- relapse_fractions[order(as.numeric(str_sub(relapse_fractions, 
                                                                            -1)), relapse_fractions)]
    for (i in 1:length(ordered_relapse_fractions)) {
      loginfo(paste0("found replase fraction(s): ", ordered_relapse_fractions[i], 
                     "\n"))
    }
  }
  else {
    loginfo(paste0("found no replase fraction(s)!\n"))
    ordered_relapse_fractions <- c()
  }
  return(c(f_sample, d_sample, ordered_relapse_fractions))
}


## Extract alelle frequencies from VCF
extractAllelFrequenciesFromVcf <- function (vcfOject) 
{
  vcf_based_AD <- geno(vcfOject)$AD
  vcf_based_AF <- geno(vcfOject)$AF
  vcf_based_sample_names <- colnames(vcf_based_AD)
  column_names <- rep(NA, times = length(vcf_based_sample_names) * 
                        4)
  for (i in 1:length(colnames(geno(vcfOject)$AD))) {
    column_names[((i * 4) - 3)] <- paste(colnames(vcf_based_AD)[i], 
                                         "_read_count_ref", sep = "")
    column_names[((i * 4) - 2)] <- paste(colnames(vcf_based_AD)[i], 
                                         "_read_count_alt", sep = "")
    column_names[((i * 4) - 1)] <- paste(colnames(vcf_based_AD)[i], 
                                         "_read_based_AF", sep = "")
    column_names[((i * 4))] <- paste(colnames(vcf_based_AD)[i], 
                                     "_mutect_2_based_AF", sep = "")
  }
  af_matrix <- matrix(nrow = nrow(vcf_based_AD), ncol = length(column_names))
  colnames(af_matrix) <- column_names
  rownames(af_matrix) <- rownames(vcf_based_AD)
  for (i in 1:length(vcf_based_sample_names)) {
    loginfo(paste0("processing ", nrow(af_matrix), " variants from sample ", 
                   vcf_based_sample_names[i]))
    for (j in 1:nrow(af_matrix)) {
      af_matrix[j, ((i * 4) - 3)] <- vcf_based_AD[[j, i]][1]
      af_matrix[j, ((i * 4) - 2)] <- vcf_based_AD[[j, i]][2]
      af_matrix[j, ((i * 4) - 1)] <- computeAllelelFrequency(af_matrix[j, 
                                                                       ((i * 4) - 3)], af_matrix[j, ((i * 4) - 2)])
      af_matrix[j, ((i * 4))] <- vcf_based_AF[[j, i]]
    }
  }
  return(af_matrix)
}


## Remove variants of specific clusters from data frame
removeClustersFromDataFrame <- function (cluster_data_frame, clusters_to_remove) 
{
  variants_to_remove_idx <- cluster_data_frame$cluster %in% 
    clusters_to_remove
  variants_to_keep_idx <- !(variants_to_remove_idx)
  filtered_cluster_data_frame <- cluster_data_frame[variants_to_keep_idx, 
  ]
  return(filtered_cluster_data_frame)
}


## Extract row ranges
extractRowRangesAsDataFrame <- function (vcfObject) 
{
  chromosome <- as.vector(seqnames(vcfObject))
  start <- start(ranges(rowRanges(vcfObject)))
  end <- end(ranges(rowRanges(vcfObject)))
  quality <- qual(vcfObject)
  reference_nucleotide <- as.character(ref(vcfObject))
  alternative_nucleotide <- as.character(unlist(alt(vcfObject)))
  filter <- filt(vcfObject)
  df <- data.frame(chr = chromosome, start = start, end = end, 
                   ref = reference_nucleotide, alt = alternative_nucleotide, 
                   quality = quality, filter = filter)
  return(df)
}


## Compute allele frequency
computeAllelelFrequency <- function (ref_read_count, alt_read_count) 
{
  AF_vector <- NA
  if (alt_read_count == 0) {
    AF_vector <- 0
  }
  else {
    AF_vector <- alt_read_count/(ref_read_count + alt_read_count)
  }
  return(AF_vector)
}


## Get remission IDs
getRemissionIds <- function (idVector) 
{
  matched_sample <- c()
  if (length(grep("[0-9]{1,}_F", idVector)) >= 1) {
    f_match <- (idVector[(grep("[0-9]{1,}_F", idVector))])
    matched_sample <- c(matched_sample, f_match)
    loginfo(paste0("found the full remission: ", f_match, 
                   "\n"))
  }
  if (length(grep("[0-9]{1,}CR[0-9]{1,}", idVector)) >= 1) {
    cr_match <- (idVector[(grep("[0-9]{1,}CR[0-9]{1,}", idVector))])
    loginfo(paste0("found the full remission: ", cr_match, 
                   "\n"))
    matched_sample <- c(matched_sample, cr_match)
  }
  if (length(grep("[0-9]{1,}FU-R[0-9]{1,}", idVector)) >= 1) {
    cr_match <- (idVector[(grep("[0-9]{1,}FU-R[0-9]{1,}", 
                                idVector))])
    loginfo(paste0("found the full remission: ", cr_match, 
                   "\n"))
    matched_sample <- c(matched_sample, cr_match)
  }
  if (length(grep("[0-9]{1,}F[0-9]{1,}", idVector)) >= 1) {
    cr_match <- (idVector[(grep("[0-9]{1,}F[0-9]{1,}", idVector))])
    loginfo(paste0("found the full remission: ", cr_match, 
                   "\n"))
    matched_sample <- c(matched_sample, cr_match)
  }
  if (is.null(matched_sample)) {
    loginfo(paste0("No matching full remission sample in: ", 
                   idVector, "\n"))
  }
  return(matched_sample)
}


## Filter variants on a timepoint's specific allele frequency
filterVariantsOnTimePointSpecificAlleleFrequency <- function (cluster_df, time_point, cluster, AF_threshold, variants_to_keep = FALSE) 
{
  cluster_specific_cluster_df_idx <- cluster_df$cluster == 
    cluster
  cluster_specific_cluster_df <- cluster_df[cluster_specific_cluster_df_idx, 
  ]
  cluster_specific_filtered_variants_clusters_cluster_AF <- cluster_specific_cluster_df[, 
                                                                                        time_point]
  cluster_specific_filtered_variants_clusters_filtered_idx <- cluster_specific_filtered_variants_clusters_cluster_AF <= 
    AF_threshold
  filtered_variant_identifiers <- cluster_specific_cluster_df[!(cluster_specific_filtered_variants_clusters_filtered_idx), 
  ]
  if (variants_to_keep) {
    filtered_variant_identifiers <- cluster_specific_cluster_df[cluster_specific_filtered_variants_clusters_filtered_idx, 
    ]
  }
  return(filtered_variant_identifiers)
}




