library(DESeq2)

# Assuming you have a DESeqDataSet `dds` already prepared and analyzed
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = col_data, design = ~ condition)
dds <- DESeq(dds)

# Obtain results for two contrasts
res1 <- results(dds, contrast = c("condition", "TreatmentA", "Control"))
res2 <- results(dds, contrast = c("condition", "TreatmentB", "Control"))

# Extract log2 fold changes and standard errors
log2FC1 <- res1$log2FoldChange
log2FC2 <- res2$log2FoldChange
SE1 <- res1$lfcSE
SE2 <- res2$lfcSE

# Compare fold changes using Z-test
z_stat <- (log2FC1 - log2FC2) / sqrt(SE1^2 + SE2^2)
p_values <- 2 * pnorm(-abs(z_stat))

# Multiple testing correction
adjusted_p_values <- p.adjust(p_values, method = "BH")

# Create summary table
comparison_results <- data.frame(
  Gene = rownames(res1),
  log2FC1 = log2FC1,
  SE1 = SE1,
  log2FC2 = log2FC2,
  SE2 = SE2,
  Z = z_stat,
  p_value = p_values,
  adjusted_p_value = adjusted_p_values
)

# View the top results
head(comparison_results)