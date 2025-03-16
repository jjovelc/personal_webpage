# Load required libraries
library(ggplot2)
library(reshape2)
library(dplyr)
library(RColorBrewer)

# Set seed for reproducibility
set.seed(123)

# Define sample names
samples <- c(paste0("Control_", 1:5), paste0("IBD_", 1:5))

# Define 15 common bacterial families
bacterial_families <- c(
  "Bacteroidaceae", "Lachnospiraceae", "Ruminococcaceae", 
  "Prevotellaceae", "Enterobacteriaceae", "Bifidobacteriaceae", 
  "Akkermansiaceae", "Rikenellaceae", "Christensenellaceae", 
  "Clostridiaceae", "Lactobacillaceae", "Veillonellaceae",
  "Erysipelotrichaceae", "Desulfovibrionaceae", "Oscillospiraceae"
)

# Create a matrix to hold abundances
abundance_matrix <- matrix(0, nrow = length(samples), ncol = length(bacterial_families))
rownames(abundance_matrix) <- samples
colnames(abundance_matrix) <- bacterial_families

# Fill with simulated data with differences between Control and IBD
# Control samples will have higher beneficial bacteria
for (i in 1:5) {
  # Controls have higher Lachnospiraceae, Ruminococcaceae, Bacteroidaceae
  abundance_matrix[paste0("Control_", i), "Lachnospiraceae"] <- runif(1, 20, 30)
  abundance_matrix[paste0("Control_", i), "Ruminococcaceae"] <- runif(1, 15, 25)
  abundance_matrix[paste0("Control_", i), "Bacteroidaceae"] <- runif(1, 15, 25)
  abundance_matrix[paste0("Control_", i), "Bifidobacteriaceae"] <- runif(1, 8, 15)
  abundance_matrix[paste0("Control_", i), "Akkermansiaceae"] <- runif(1, 5, 10)
  
  # Fill the rest with lower abundances
  for (family in setdiff(bacterial_families, c("Lachnospiraceae", "Ruminococcaceae", "Bacteroidaceae", "Bifidobacteriaceae", "Akkermansiaceae"))) {
    abundance_matrix[paste0("Control_", i), family] <- runif(1, 0.5, 5)
  }
  
  # IBD samples have higher Enterobacteriaceae, Clostridiaceae, and lower diversity
  abundance_matrix[paste0("IBD_", i), "Enterobacteriaceae"] <- runif(1, 20, 35)
  abundance_matrix[paste0("IBD_", i), "Clostridiaceae"] <- runif(1, 10, 20)
  abundance_matrix[paste0("IBD_", i), "Bacteroidaceae"] <- runif(1, 10, 20)
  abundance_matrix[paste0("IBD_", i), "Lachnospiraceae"] <- runif(1, 5, 15)
  abundance_matrix[paste0("IBD_", i), "Prevotellaceae"] <- runif(1, 8, 12)
  
  # Fill the rest with lower abundances
  for (family in setdiff(bacterial_families, c("Enterobacteriaceae", "Clostridiaceae", "Bacteroidaceae", "Lachnospiraceae", "Prevotellaceae"))) {
    abundance_matrix[paste0("IBD_", i), family] <- runif(1, 0.1, 3)
  }
}

# Normalize to ensure each sample sums to 100%
for (i in 1:nrow(abundance_matrix)) {
  abundance_matrix[i, ] <- abundance_matrix[i, ] / sum(abundance_matrix[i, ]) * 100
}

# Convert to data frame for ggplot
abundance_df <- as.data.frame(abundance_matrix)
abundance_df$Sample <- rownames(abundance_df)

# Add group information
abundance_df$Group <- ifelse(grepl("Control", abundance_df$Sample), "Control", "IBD")

# Reshape to long format for ggplot
abundance_long <- reshape2::melt(abundance_df, id.vars = c("Sample", "Group"),
                                 variable.name = "Family", value.name = "RelativeAbundance")

# Order samples by group
abundance_long$Sample <- factor(abundance_long$Sample, 
                                levels = c(paste0("Control_", 1:5), paste0("IBD_", 1:5)))

# Generate a colorful palette for the families
n_colors <- length(bacterial_families)
family_colors <- colorRampPalette(brewer.pal(8, "Set3"))(n_colors)

# Create stacked barplot
ggplot(abundance_long, aes(x = Sample, y = RelativeAbundance, fill = Family)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = family_colors) +
  labs(
    title = "Bacterial Family Composition in Control vs IBD Samples",
    x = "",
    y = "Relative Abundance (%)",
    fill = "Bacterial Family"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5)
  ) +
  facet_grid(~ Group, scales = "free_x", space = "free_x")

# Save the plot
ggsave("family_level_barplot.png", width = 12, height = 8, dpi = 300)

# Print the first few rows of the data for verification
head(abundance_long)