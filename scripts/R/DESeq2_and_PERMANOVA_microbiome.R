
############################################################################
#         Author: Juan Jovel (juan@dayhoff.ai)                            #
###########################################################################

# To install DESeq2, run the following commands:
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")

# All other required libraries:
# "xlsx", "ecodist", "ggplot2"... etc
# can be installed with command:
# install.packages("library") [for individual libraries]
# or: install.packages(c("library1","library2","libraryn")) [for multiple libraries]

library(DESeq2)
library(ecodist)
library(RColorBrewer)
library(pheatmap)
library(ape)
library(vegan)
library(ggplot2)
library(EnhancedVolcano)

rm(list = ls())
setwd("/Users/juanjovel/OneDrive/jj/UofC/data_analysis/heshanthi/DESeq2_analysis")

prefix = 'Avian_ASV_DMV'
# Create names for output files
heatmap_file <- paste0(prefix, '_HCheatmap.pdf')
PCoA_file    <- paste0(prefix, '_PCoA.pdf')
volcano_file <- paste0(prefix, '_volcanoPlot.pdf')
res_file     <- paste0(prefix, '_q0.05_FC2.tsv')
allRes_file  <- paste0(prefix, '_AllRes.tsv')

metadata_file <- 'metadata.tsv'
taxa_file     <- 'GTDB_RDP.tsv'
otus_file     <- 'seqtab_nochimeras.csv'

metadata  <- read.table(metadata_file, header = T, row.names = 1, "\t")
taxonomy  <- read.table(taxa_file, header = T, row.names = 1)
otu_table <- read.csv(otus_file, header = T, row.names = 1)

# Ensure alignment between taxonomy and OTU table
if (!all(rownames(taxonomy) %in% colnames(otu_table))) {
  stop("Mismatch between taxonomy and OTU table")
}

# Reorder taxonomy to match the order of colnames in otu_table
taxonomy <- taxonomy[colnames(otu_table), ]

# List of taxonomic levels and their corresponding abbreviations
taxa_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
abbreviations <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")
names(abbreviations) <- taxa_levels

# Function to extract the last non-NA taxon for a given row
extract_last_taxon <- function(row) {
  last_non_na <- tail(which(!is.na(row)), 1)
  if (length(last_non_na) == 0) {
    return("Unknown")
  } else {
    taxon_level <- taxa_levels[last_non_na]
    taxon_abbr <- abbreviations[taxon_level]
    taxon_name <- row[last_non_na]
    
    if (taxon_level == "Species") {
      previous_taxon_name <- row[last_non_na - 1]
      return(paste0(taxon_abbr, previous_taxon_name, "_", taxon_name))
    } else {
      return(paste0(taxon_abbr, taxon_name))
    }
  }
}

# Apply the function to each row of the taxonomy data frame
last_taxon_names <- apply(taxonomy, 1, extract_last_taxon)

# Replace NAs with "Unknown" or any placeholder if needed
last_taxon_names[is.na(last_taxon_names)] <- "Unknown"

# Replace the colnames of otu_table with the extracted taxonomy names
colnames(otu_table) <- last_taxon_names

# Filter metadata
without_controls <- metadata$Infection == "Control" | metadata$Infection == "IBV"
clean_metadata   <- metadata[without_controls,]
keep_rows        <- !(row.names(clean_metadata) == "s_325")
clean_metadata   <-  clean_metadata[keep_rows,]
clean_metadata$InfectionTime <- paste0(clean_metadata$Infection, "_", clean_metadata$Time)

otu_table_t <- t(otu_table)
current_sample_names <- colnames(otu_table_t)
new_sample_names <- paste0('s_', current_sample_names)
colnames(otu_table_t) <- new_sample_names

# Remove low abundance ASVs
# Logical vector to retain rows not matching the criteria
data_5up          <- subset(otu_table_t, rowMeans(otu_table_t) >= 5)
filtered_data_5up <- data_5up[keep_rows, ]

metadata_trachea <- clean_metadata[clean_metadata$Tissue == "Trachea", ]
data_trachea     <- filtered_data_5up[, row.names(metadata_trachea)]

metadata_cecum   <- clean_metadata[clean_metadata$Tissue == "Cecum", ]
data_cecum       <- filtered_data_5up[, row.names(metadata_cecum)]

# Remove decimals because DESeq2 only handles integer
data_trachea <- round(as.matrix(data_trachea, digits = 0))
data_cecum   <- round(as.matrix(data_cecum, digits = 0))

# Create cecum control only
metadata_cecum_control   <- metadata_cecum[metadata_cecum$Infection == "Control", ]
data_cecum_control <- data_cecum[, row.names(metadata_cecum_control)]    

                                           # Import data into DESEq2 object and define your design
dds_trachea <- DESeqDataSetFromMatrix(countData = data_trachea, colData = metadata_trachea, 
                                 design =~ InfectionTime)

dds_cecum   <- DESeqDataSetFromMatrix(countData = data_cecum, colData = metadata_cecum, 
                                      design =~ InfectionTime)

dds_cecum_control <- DESeqDataSetFromMatrix(countData = data_cecum_control, colData = metadata_cecum_control, 
                                      design =~ InfectionTime)

# Apply logarithmic transformation to create exploratory plots 
# (Hierarchical clustering and PCoA)
dds_trachea <- estimateSizeFactors(dds_trachea)
dds_cecum  <- estimateSizeFactors(dds_cecum)
dds_cecum_control  <- estimateSizeFactors(dds_cecum_control)

vst_trachea <- varianceStabilizingTransformation(dds_trachea, fitType = 'local')
vst_cecum   <- varianceStabilizingTransformation(dds_cecum, fitType = 'local')
vst_cecum_control   <- varianceStabilizingTransformation(dds_cecum_control, fitType = 'local')

rld_trachea <- rlog(dds_trachea, fitType = 'local')
rld_cecum   <- rlog(dds_cecum, fitType = 'local')
rld_cecum_control   <- rlog(dds_cecum_control, fitType = 'local')

# Hierarchical clustering
bc_dist_trachea <- bcdist(t(assay(rld_trachea)))
bc_dist_cecum   <- bcdist(t(assay(rld_cecum)))
bc_dist_cecum_control   <- bcdist(t(assay(rld_cecum_control)))

sampleDistMatrix_trachea <- as.matrix(bc_dist_trachea)
sampleDistMatrix_cecum   <- as.matrix(bc_dist_cecum)
sampleDistMatrix_cecum_control   <- as.matrix(bc_dist_cecum_control)

generate_pcoa <- function(bc_dist, metadata, PCoA_file)
{
  # Generate Principal Coordinates (PCoA)
  PCo <- pco(bc_dist)
  
  # Put principal coordinates into a dataframe
  
  PCo_df <- data.frame(PC1 = PCo$vectors[,1], PC2 = PCo$vectors[,2])
  col_names <- c("PCoA1", "PCoA2")
  colnames(PCo_df) <- col_names
  
  # Plot PCoA
  pdf(file=PCoA_file, width = 10, height = 8)
  p <- ggplot(PCo_df, aes(PCoA1, PCoA2, color=metadata$Time)) +
    geom_point(shape=19, size=7) +
    scale_color_manual(values=c("firebrick1", "dodgerblue", "limegreen")) + 
    geom_text(aes(label=metadata_cecum_control$Bird, color='black')) +
    theme_bw()
  print(p)
  dev.off()
}

file_name <- paste0(prefix, "_trachea_PCoA.pdf")
generate_pcoa(bc_dist_trachea, metadata_trachea, file_name)

file_name <- paste0(prefix, "_cecum_PCoA.pdf")
generate_pcoa(bc_dist_cecum, metadata_cecum, file_name)

file_name <- paste0(prefix, "_cecum_control_PCoA.pdf")
generate_pcoa(bc_dist_cecum_control, metadata_cecum_control, file_name)



generate_pcoa_infTime <- function(bc_dist, metadata, PCoA_file) {
  # Ensure the Time variable is a factor
  metadata$Time <- as.factor(metadata$Time)
  
  # Generate Principal Coordinates (PCoA)
  PCo <- pco(bc_dist)
  
  # Put principal coordinates into a dataframe
  PCo_df <- data.frame(PC1 = PCo$vectors[, 1], PC2 = PCo$vectors[, 2])
  col_names <- c("PCoA1", "PCoA2")
  colnames(PCo_df) <- col_names

  # Combine PCoA data with metadata
  PCo_df <- cbind(PCo_df, metadata)
  
  # Plot PCoA
  pdf(file = PCoA_file, width = 10, height = 8)
  p <- ggplot(PCo_df, aes(x = PCoA1, y = PCoA2, color = Infection, shape = Time)) +
    geom_point(size = 7) +
    scale_color_manual(values = c("firebrick1", "dodgerblue")) +
    scale_shape_manual(values = c(16, 17, 18)) +  # Shapes for different levels of Time
    geom_text(aes(label = row.names(PCo_df)), color = 'black', hjust = -0.1, vjust = -0.1) +
    theme_bw()
  print(p)
  dev.off()
}

file_name <- paste0(prefix, "_cecum_PCoA_perInfectionTime.pdf")
generate_pcoa_infTime(bc_dist_cecum, metadata_cecum, file_name)

file_name <- paste0(prefix, "_trachea_PCoA_perInfectionTime.pdf")
generate_pcoa_infTime(bc_dist_trachea, metadata_trachea, file_name)

# Generate PCoA plots per time point function with BH correction for p-values
generate_pcoa_per_time_point <- function(bc_dist, metadata, PCoA_file_prefix) {
  # Ensure the Time variable is a factor
  metadata$Time <- as.factor(metadata$Time)
  
  # Generate Principal Coordinates (PCoA)
  PCo <- pco(bc_dist)
  
  # Put principal coordinates into a dataframe
  PCo_df <- data.frame(PC1 = PCo$vectors[, 1], PC2 = PCo$vectors[, 2])
  col_names <- c("PCoA1", "PCoA2")
  colnames(PCo_df) <- col_names
  
  # Combine PCoA data with metadata
  PCo_df <- cbind(PCo_df, metadata)
  
  # Loop through each time point and collect p-values
  time_points <- levels(metadata$Time)
  p_values <- c()
  
  for (time_point in time_points) {
    # Filter data for the current time point
    PCo_df_time <- subset(PCo_df, Time == time_point)
    
    # Subset the distance matrix
    sample_names <- rownames(PCo_df_time)
    bc_dist_time <- as.dist(as.matrix(bc_dist)[sample_names, sample_names])
    
    # Replace negative dissimilarities with zero
    bc_dist_time <- as.matrix(bc_dist_time)
    bc_dist_time[bc_dist_time < 0] <- 0
    bc_dist_time <- as.dist(bc_dist_time)
    
    # Perform PERMANOVA analysis using adonis2
    if (length(unique(PCo_df_time$Infection)) > 1) {
      perm <- adonis2(bc_dist_time ~ Infection, data = PCo_df_time)
      p_value <- perm$`Pr(>F)`[1]
    } else {
      p_value <- NA
    }
    
    # Collect p-values
    p_values <- c(p_values, p_value)
  }
  
  # Apply BH correction to p-values
  p_values_adj <- p.adjust(p_values, method = "BH")
  
  # Generate plots with adjusted p-values
  for (i in seq_along(time_points)) {
    time_point <- time_points[i]
    p_value_adj <- p_values_adj[i]
    
    # Filter data for the current time point
    PCo_df_time <- subset(PCo_df, Time == time_point)
    
    # Subset the distance matrix
    sample_names <- rownames(PCo_df_time)
    bc_dist_time <- as.dist(as.matrix(bc_dist)[sample_names, sample_names])
    
    # Replace negative dissimilarities with zero
    bc_dist_time <- as.matrix(bc_dist_time)
    bc_dist_time[bc_dist_time < 0] <- 0
    bc_dist_time <- as.dist(bc_dist_time)
    
    # Define the output file name
    PCoA_file <- paste0(PCoA_file_prefix, "_Time", time_point, ".pdf")
    
    # Plot PCoA
    pdf(file = PCoA_file, width = 10, height = 8)
    p <- ggplot(PCo_df_time, aes(x = PCoA1, y = PCoA2, color = Infection)) +
      geom_point(size = 7) +
      scale_color_manual(values = c("firebrick1", "dodgerblue")) +
      geom_text(aes(label = row.names(PCo_df_time)), color = 'black', hjust = -0.1, vjust = -0.1) +
      ggtitle(paste("PCoA Plot for Time Point", time_point, "\nPERMANOVA p=", format(p_value_adj, digits = 4))) +
      theme_bw()
    print(p)
    dev.off()
  }
}

file_name_prefix <- paste0(prefix, "_cecum_PCoA")
generate_pcoa_per_time_point(bc_dist_cecum, metadata_cecum, file_name_prefix)

file_name_prefix <- paste0(prefix, "_trachea_PCoA")
generate_pcoa_per_time_point(bc_dist_trachea, metadata_trachea, file_name_prefix)


# Transform data using the harmonic mean
gm_mean <- function(x, na.rm=T) {
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

geoMean   <- apply(counts(dds_cecum), 1, gm_mean)
dds_cecum <- estimateSizeFactors(dds_cecum, geoMean = geoMean)

geoMean     <- apply(counts(dds_trachea), 1, gm_mean)
dds_trachea <- estimateSizeFactors(dds_trachea, geoMean = geoMean)

list_dds <- list(dds_cecum, dds_trachea)
names(list_dds) <- c("dds_cecum", "dds_trachea")

contrasts <- list(
  ctrl_vs_BVI_1dpi  = c(group = "InfectionTime", "IBV_1", "Control_1"),
  ctrl_vs_BVI_4dpi  = c(group = "InfectionTime", "IBV_4", "Control_4"),
  ctrl_vs_BVI_10dpi = c(group = "InfectionTime", "IBV_10", "Control_10")
)

makeVolcanoPlot <- function(df, vp_file){
  keyvals <- ifelse(
    df$log2FoldChange < -1 & df$padj < 0.05, 'forestgreen',
    ifelse(df$log2FoldChange > 1 & df$padj < 0.05, 'firebrick1',
           'dodgerblue1'))
  
  keyvals[is.na(keyvals)] <- 'black'
  names(keyvals)[keyvals == 'firebrick1'] <- 'High'
  names(keyvals)[keyvals == 'dodgerblue1'] <- 'small FC'
  names(keyvals)[keyvals == 'forestgreen'] <- 'Low'
  
  evp <- EnhancedVolcano(df,
                         lab = NA,
                         x = 'log2FoldChange',
                         y = 'pvalue',
                         pCutoff = 0.01,
                         FCcutoff = 1,
                         colCustom = keyvals,
                         title = NULL,
                         subtitle = NULL,
                         colAlpha = 0.6,
                         shape = 19,
                         pointSize = 5,
                         labSize = 3,
                         gridlines.major = TRUE,
                         gridlines.minor = TRUE,
                         border = 'full',
                         borderWidth = 0.25,
                         borderColour = 'black')
  
  png(vp_file)
  print(evp)
  dev.off()
}

exportTable <- function(df, filename){
  write.table(df, filename, quote = F, sep='\t', row.names = F)
}

for (dds_name in names(list_dds)) {
  tissue <- gsub("dds_", "", dds_name)
  
  # Extract the object using its name from the list
  object <- list_dds[[dds_name]]
  object <- DESeq(object, fitType = 'local')
  
  for (contrast_name in names(contrasts)) {
    my_results        <- paste(tissue, contrast_name, "res", sep = "_")
    cat("Contrast being processed: ", my_results, "\n")
    contrast_value    <- contrasts[[contrast_name]]
    my_result         <- results(object, contrast = contrast_value)
    resdata           <- merge(as.data.frame(my_result), as.data.frame(counts(object, normalized=TRUE)), by="row.names", sort=FALSE)
    names(resdata)[1] <- "taxa"
    outfile_norm_name <- paste(tissue, contrast_name, "allData.tsv", sep = "_")
    exportTable(resdata, outfile_norm_name)
    sign_results      <- subset(resdata, padj < 0.05 & abs(log2FoldChange) > 1)
    outfile_name      <- paste(tissue, contrast_name, "q0.05_FC2.tsv", sep = "_")
    exportTable(sign_results, outfile_name)
    vp_file           <- paste(tissue, contrast_name, "volcano_plot.png", sep = "_")
    makeVolcanoPlot(resdata, vp_file)
  }
}






