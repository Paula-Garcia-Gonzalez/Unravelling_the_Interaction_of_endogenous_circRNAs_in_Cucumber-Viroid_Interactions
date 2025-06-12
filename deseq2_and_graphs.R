# Install packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", dependencies = TRUE)
}
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2", dependencies = TRUE)
}
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap", dependencies = TRUE)
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", dependencies = TRUE)
}
if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  install.packages("RColorBrewer", dependencies = TRUE)
}
if (!requireNamespace("plotly", quietly = TRUE)) {
  install.packages("plotly", dependencies = TRUE)
}
if (!requireNamespace("gridExtra", quietly = TRUE)) {
  install.packages("gridExtra", dependencies = TRUE)
}

# Load packages
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(plotly)
library(gridExtra)

# Read count table
counts <- read.csv("final_counts_zeros.csv", row.names = 1)
head(counts)

# Split table into 3 by time points
table_1 <- counts[, grepl("^R00_1", colnames(counts))]  # Select columns starting with "R00_1"
table_2 <- counts[, grepl("^R00_2", colnames(counts))]  # Select columns starting with "R00_2"
table_3 <- counts[, grepl("^R00_3", colnames(counts))]  # Select columns starting with "R00_3"

# Save the tables as separate files
write.table(table_1, "table_1.csv", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(table_2, "table_2.csv", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(table_3, "table_3.csv", sep = "\t", quote = FALSE, row.names = TRUE)

# Sample information
coldata_1 <- read.csv("coldata_1.csv", row.names = 1)
coldata_2 <- read.csv("coldata_2.csv", row.names = 1)
coldata_3 <- read.csv("coldata_3.csv", row.names = 1)

# Convert Condition column to factor
coldata_1$Condition <- as.factor(coldata_1$Condition)
coldata_2$Condition <- as.factor(coldata_2$Condition)
coldata_3$Condition <- as.factor(coldata_3$Condition)

# Filtering function
filter_table <- function(table) {
  table[rowSums(table >= 5) >= 3, ]  # Keep rows with at least 3 columns having counts >= 5
}

# Filter the tables
table_1_filtered <- filter_table(table_1)
table_2_filtered <- filter_table(table_2)
table_3_filtered <- filter_table(table_3)

# Save filtered tables if needed
write.table(table_1_filtered, "table_1_filtered.csv", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(table_2_filtered, "table_2_filtered.csv", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(table_3_filtered, "table_3_filtered.csv", sep = "\t", quote = FALSE, row.names = TRUE)

########################################################################################################################

# TIME POINT 1:

# Create DESeqDataSet
dds_1 <- DESeqDataSetFromMatrix(countData = table_1_filtered,
                                colData = coldata_1,
                                design = ~ Condition)

# Specify control group as reference
dds_1$Condition <- relevel(dds_1$Condition, ref = "control")

# Run DESeq
dds_out_1 <- DESeq(dds_1)

# View results names
results_names_1 <- resultsNames(dds_out_1)
head(results_names_1)

# Loop through results
for (comp in results_names_1) {
  if (comp != 'Intercept') {
    # Get results
    res_deseq_1 <- results(dds_out_1, name = comp, alpha = 0.05)
    resLFC_1 <- lfcShrink(dds_out_1, coef = comp, res = res_deseq_1)
    
    # Optional: save or analyze results
    print(paste("Results for comparison:", comp))
    print(summary(resLFC_1))
  }
}

# Export results to CSV
write.csv(as.data.frame(resLFC_1), file = paste0("deseq2_results_", comp, ".csv"))

# Consistent color scheme for conditions
condition_colors <- c("control" = "#00FF00", "infected" = "#FF69B4")

######################################
# PCA (Principal Component Analysis) #
######################################

# Obtain PCA data
vsd_1 <- varianceStabilizingTransformation(dds_out_1, blind = FALSE) # Normalization preserving experimental design
pca_data_1 <- plotPCA(vsd_1, intgroup = "Condition", returnData = TRUE) # Get data for ggplot
percentVar_1 <- round(100 * attr(pca_data_1, "percentVar")) # Explained variance percentages

# Plot using ggplot2
library(ggplot2)
pca_plot <- ggplot(pca_data_1, aes(PC1, PC2, color = Condition)) + 
  geom_point(size = 3) +
  geom_text(aes(label = name), vjust = -0.5, size = 3) + # Add labels
  scale_color_manual(values = condition_colors) + # Apply color scheme manually
  xlab(paste0("PC1: ", percentVar_1[1], "% variance")) + 
  ylab(paste0("PC2: ", percentVar_1[2], "% variance")) + 
  theme_minimal() +
  ggtitle("PCA for Time Point 1") # Add title

# Display the plot
print(pca_plot)

# Save the plot with white background
ggsave("PCA_plot_1_labels.png", plot = pca_plot, width = 8, height = 6, dpi = 300, bg = "white")

# Create 3D PCA plot (first attempt - note: this part has a small error you fix below)
ppca_3d_plot_1 <- plot_ly(pca_data_3d_1, 
                          x = ~PC1, 
                          y = ~PC2, 
                          z = ~PC3, 
                          color = ~Condition, 
                          colors = c("control" = "#00FF00", "infected" = "#FF69B4"),
                          type = "scatter3d", 
                          mode = "markers+text", 
                          marker = list(size = 5),
                          text = rownames(pca_data_3d_1), 
                          textposition = "top center", 
                          hoverinfo = "text+x+y+z", 
                          hovertext = ~paste("Name:", rownames(pca_data_3d_1),
                                             "<br>Condition:", Condition,
                                             "<br>PC1:", round(PC1, 2), 
                                             "<br>PC2:", round(PC2, 2),
                                             "<br>PC3:", round(PC3, 2))) %>%
  layout(scene = list(
    xaxis = list(title = "PC1"),
    yaxis = list(title = "PC2"),
    zaxis = list(title = "PC3")),
    title = "3D PCA with Labels")

pca_3d_plot

htmlwidgets::saveWidget(as_widget(pca_3d_plot), "PCA_3D_Time_1_labels.html")

##########
# 3D PCA #
##########

# Obtain PCA data
vsd_1 <- varianceStabilizingTransformation(dds_out_1, blind = FALSE) # Normalization
pca_res_1 <- prcomp(t(assay(vsd_1))) # Compute PCA using transposed data

# Create data frame with first three principal components
pca_data_3d_1 <- as.data.frame(pca_res_1$x[, 1:3]) # Select PC1, PC2 and PC3
pca_data_3d_1$Condition <- colData(vsd_1)$Condition # Add conditions column
colnames(pca_data_3d_1)[1:3] <- c("PC1", "PC2", "PC3") # Rename columns

# Explained variance percentages
percentVar_3d_1 <- round(100 * pca_res_1$sdev^2 / sum(pca_res_1$sdev^2), 1)

# Create 3D PCA plot (corrected version)
pca_3d_plot_1 <- plot_ly(pca_data_3d_1, 
                         x = ~PC1, 
                         y = ~PC2, 
                         z = ~PC3, 
                         color = ~Condition, 
                         colors = c("control" = "#00FF00", "infected" = "#FF69B4"),
                         type = "scatter3d", 
                         mode = "markers+text", 
                         marker = list(size = 5),
                         text = rownames(pca_data_3d_1), 
                         textposition = "top center", 
                         hoverinfo = "text+x+y+z", 
                         hovertext = ~paste("Name:", rownames(pca_data_3d_1),
                                            "<br>Condition:", Condition,
                                            "<br>PC1:", round(PC1, 2), 
                                            "<br>PC2:", round(PC2, 2),
                                            "<br>PC3:", round(PC3, 2))) %>%
  layout(scene = list(
    xaxis = list(title = "PC1"),
    yaxis = list(title = "PC2"),
    zaxis = list(title = "PC3")),
    title = "3D PCA with Labels")

pca_3d_plot_1

htmlwidgets::saveWidget(as_widget(pca_3d_plot_1), "PCA_3D_Time_1_labels.html")

################
# VOLCANO PLOT #
################

# Define thresholds
log2fc_threshold <- 1  # Change this as needed
padj_threshold <- 0.05

# Classify genes
resLFC_1$diffexpressed <- "Not differentially expressed"
resLFC_1$diffexpressed[resLFC_1$log2FoldChange > log2fc_threshold & resLFC_1$padj < padj_threshold] <- "Up-regulated"
resLFC_1$diffexpressed[resLFC_1$log2FoldChange < -log2fc_threshold & resLFC_1$padj < padj_threshold] <- "Down-regulated"

# INDEPENDENT

# Calculate the maximum absolute log2FoldChange value
max_abs_value_1 <- max(abs(resLFC_1$log2FoldChange), na.rm = TRUE)

# Create the volcano plot
volcano_1 <- ggplot(data = as.data.frame(resLFC_1), aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed)) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0, col = "grey") +  # Vertical line at x = 0
  geom_hline(yintercept = -log10(padj_threshold), col = "grey") +  # Horizontal line for adjusted p-value threshold
  scale_color_manual(values = c("Not differentially expressed" = "snow3", 
                                "Up-regulated" = "#FF6F61", 
                                "Down-regulated" = "#6EC5E9"),
                     name = "Differential Expression",
                     labels = c("Not differentially expressed", "Up-regulated", "Down-regulated")) +
  labs(x = "Log2 Fold Change", 
       y = "-Log10 (Adjusted p-value)", 
       title = "Volcano Plot for Time Point 1") +
  scale_x_continuous(limits = c(-max_abs_value_1 - 0.2, max_abs_value_1 + 0.2)) +  # Limit X axis
  guides(color = guide_legend(title = NULL)) +  # Remove legend title
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate X axis text
  theme_bw()  # Use white background theme

# Display the plot
print(volcano_1)

# Save the plot as PNG
ggsave("volcano_plot_1.png", plot = volcano_1, width = 8, height = 6, dpi = 300, bg = "white")

# COMPARISON BETWEEN THE 3 TIME POINTS (SAME AXES)

# Calculate the maximum absolute log2FoldChange value across all time points
max_log2FoldChange_all <- max(c(abs(resLFC_1$log2FoldChange), abs(resLFC_2$log2FoldChange), abs(resLFC_3$log2FoldChange)))
max_neg_log10_padj_all <- max(c(-log10(resLFC_1$padj), -log10(resLFC_2$padj), -log10(resLFC_3$padj)), na.rm = TRUE)  # Ignore NA values

# Create the volcano plot with adjusted axes
volcano_1_adjusted <- ggplot(data = as.data.frame(resLFC_1), aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed)) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0, col = "grey") +
  geom_hline(yintercept = -log10(padj_threshold), col = "grey") +
  scale_color_manual(values = c("Not differentially expressed" = "snow3", 
                                "Up-regulated" = "#FF6F61", 
                                "Down-regulated" = "#6EC5E9"),
                     name = "Differential Expression",
                     labels = c("Not differentially expressed", "Up-regulated", "Down-regulated")) +
  labs(x = "Log2 Fold Change", 
       y = "-Log10 (Adjusted p-value)", 
       title = "Volcano Plot for Time Point 1") +
  scale_x_continuous(limits = c(-max_log2FoldChange_all - 0.2, max_log2FoldChange_all + 0.2)) +  
  scale_y_continuous(limits = c(0, max_neg_log10_padj_all + 1)) +  
  guides(color = guide_legend(title = NULL)) +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
  theme_bw()

# Display the adjusted plot
print(volcano_1_adjusted)

# Save the adjusted plot as PNG
ggsave("volcano_plot_1_adjusted.png", plot = volcano_1_adjusted, width = 8, height = 6, dpi = 300, bg = "white")

########################################################################################################################

# TIME POINT 2:

# Create DESeqDataSet
dds_2 <- DESeqDataSetFromMatrix(countData = table_2_filtered,
                                colData = coldata_2,
                                design = ~ Condition)

# Set control group as reference
dds_2$Condition <- relevel(dds_2$Condition, ref = "control")

# Run DESeq
dds_out_2 <- DESeq(dds_2)

# View results names
results_names_2 <- resultsNames(dds_out_2)
head(results_names_2)

# Iterate through results
for (comp in results_names_2) {
  if (comp != 'Intercept') {
    # Get results
    res_deseq_2 <- results(dds_out_2, name = comp, alpha = 0.05)
    resLFC_2 <- lfcShrink(dds_out_2, coef = comp, res = res_deseq_2)
    
    # Optional: save or analyze results
    print(paste("Results for comparison:", comp))
    print(summary(resLFC_2))
  }
}

# Export results to CSV
write.csv(as.data.frame(resLFC_2), file = paste0("deseq2_results_", comp, "_T2.csv"))

# Consistent color scheme for conditions
condition_colors <- c("control" = "#00FF00", "infected" = "#FF69B4")


# Get PCA data
vsd_2 <- varianceStabilizingTransformation(dds_out_2, blind = FALSE) # Normalization preserving experimental design
pca_data_2 <- plotPCA(vsd_2, intgroup = "Condition", returnData = TRUE) # Extract PCA data for ggplot
percentVar_2 <- round(100 * attr(pca_data_2, "percentVar")) # Percent variance explained

# Plot PCA with ggplot2
library(ggplot2)
pca_plot <- ggplot(pca_data_2, aes(PC1, PC2, color = Condition)) + 
  geom_point(size = 3) +
  geom_text(aes(label = name), vjust = -0.5, size = 3) + # Add labels
  scale_color_manual(values = condition_colors) + # Apply color scheme manually
  xlab(paste0("PC1: ", percentVar_2[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_2[2], "% variance")) +
  theme_minimal() +
  ggtitle("PCA for Time Point 2")

# Display plot
print(pca_plot)

# Save plot with white background
ggsave("PCA_plot_2.png", plot = pca_plot, width = 8, height = 6, dpi = 300, bg = "white")


##########
# 3D PCA #
##########

# Get PCA data
vsd_2 <- varianceStabilizingTransformation(dds_out_2, blind = FALSE) # Normalization
pca_res_2 <- prcomp(t(assay(vsd_2))) # Compute PCA on transposed data

# Create dataframe with first three principal components
pca_data_3d_2 <- as.data.frame(pca_res_2$x[, 1:3]) # Select PC1, PC2, PC3
pca_data_3d_2$Condition <- colData(vsd_2)$Condition # Add condition as column
colnames(pca_data_3d_2)[1:3] <- c("PC1", "PC2", "PC3") # Rename columns

# Variance explained percentages
percentVar_3d_2 <- round(100 * pca_res_2$sdev^2 / sum(pca_res_2$sdev^2), 1)

# Create 3D PCA plot
pca_3d_plot_2 <- plot_ly(pca_data_3d_2, 
                         x = ~PC1, 
                         y = ~PC2, 
                         z = ~PC3, 
                         color = ~Condition, 
                         colors = c("control" = "#00FF00", "infected" = "#FF69B4"),
                         type = "scatter3d", 
                         mode = "markers+text",
                         marker = list(size = 5),
                         text = rownames(pca_data_3d_2),
                         textposition = "top center",
                         hoverinfo = "text+x+y+z",
                         hovertext = ~paste("Name:", rownames(pca_data_3d_2),
                                            "<br>Condition:", Condition,
                                            "<br>PC1:", round(PC1, 2), 
                                            "<br>PC2:", round(PC2, 2),
                                            "<br>PC3:", round(PC3, 2))) %>%
  layout(scene = list(
    xaxis = list(title = "PC1"),
    yaxis = list(title = "PC2"),
    zaxis = list(title = "PC3")),
    title = "3D PCA with Labels")

# Display 3D PCA plot
pca_3d_plot

# Save interactive 3D PCA as HTML file
htmlwidgets::saveWidget(as_widget(pca_3d_plot_2), "PCA_3D_Time_2_labels.html")


################
# VOLCANO PLOT #
################

# Define thresholds
log2fc_threshold <- 1  # Change this as needed
padj_threshold <- 0.05

# Classify genes
resLFC_2$diffexpressed <- "Not differentially expressed"
resLFC_2$diffexpressed[resLFC_2$log2FoldChange > log2fc_threshold & resLFC_2$padj < padj_threshold] <- "Up-regulated"
resLFC_2$diffexpressed[resLFC_2$log2FoldChange < -log2fc_threshold & resLFC_2$padj < padj_threshold] <- "Down-regulated"

# INDEPENDENT PLOT

# Calculate the maximum absolute value of log2FoldChange
max_abs_value_2 <- max(abs(resLFC_2$log2FoldChange), na.rm = TRUE)

# Create volcano plot
volcano_2 <- ggplot(data = as.data.frame(resLFC_2), aes(x = 
                                                          log2FoldChange, y = -log10(padj), col = diffexpressed)) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0, col = "grey") +  # Vertical line at x = 0
  geom_hline(yintercept = -log10(padj_threshold), col = "grey") +  # Horizontal line for adjusted p-value threshold
  scale_color_manual(values = c("Not differentially expressed" = 
                                  "snow3", "Up-regulated" = "#FF6F61", "Down-regulated" = "#6EC5E9"),
                     name = "Differential Expression",
                     labels = c("Not differentially expressed", "Up-regulated", "Down-regulated")) +
  labs(x = "Log2 Fold Change", 
       y = "-Log10 (Adjusted p-value)", 
       title = "Volcano Plot for Time Point 2") +
  scale_x_continuous(limits = c(-max_abs_value_2 - 0.2, max_abs_value_2 + 0.2)) +  # Limit X axis
  guides(color = guide_legend(title = NULL)) +  # Remove legend title
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate X axis labels
  theme_bw()  # Use white background theme

# Display the plot
print(volcano_2)

# Save the plot as PNG file
ggsave("volcano_plot_2.png", plot = volcano_2, width = 8, height = 6, dpi = 300, bg = "white")


# COMPARISON ACROSS THE 3 TIME POINTS (SAME AXES)

# Calculate the maximum absolute value of log2FoldChange across all time points
max_log2FoldChange_all <- max(c(abs(resLFC_1$log2FoldChange), abs(resLFC_2$log2FoldChange), abs(resLFC_3$log2FoldChange)))
max_neg_log10_padj_all <- max(c(-log10(resLFC_1$padj), -log10(resLFC_2$padj), -log10(resLFC_3$padj)), na.rm = TRUE)  # Ignore NA values

# Create adjusted volcano plot for comparison
volcano_2_adjusted <- ggplot(data = as.data.frame(resLFC_2), aes(x = 
                                                                   log2FoldChange, y = -log10(padj), col = diffexpressed)) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0, col = "grey") +  # Vertical line at x = 0
  geom_hline(yintercept = -log10(padj_threshold), col = "grey") +  # Horizontal line for adjusted p-value threshold
  scale_color_manual(values = c("Not differentially expressed" = 
                                  "snow3", "Up-regulated" = "#FF6F61", "Down-regulated" = "#6EC5E9"),
                     name = "Differential Expression",
                     labels = c("Not differentially expressed", "Up-regulated", "Down-regulated")) +
  labs(x = "Log2 Fold Change", 
       y = "-Log10 (Adjusted p-value)", 
       title = "Volcano Plot for Time Point 2") +
  scale_x_continuous(limits = c(-max_log2FoldChange_all - 0.2, max_log2FoldChange_all + 0.2)) +  # Limit X axis
  scale_y_continuous(limits = c(0, max_neg_log10_padj_all + 1)) +  # Limit Y axis
  guides(color = guide_legend(title = NULL)) +  # Remove legend title
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate X axis labels
  theme_bw()  # Use white background theme

# Display the plot
print(volcano_2_adjusted)

# Save the plot as PNG file
ggsave("volcano_plot_2_adjusted.png", plot = volcano_2_adjusted, width = 8, height = 6, dpi = 300, bg = "white")

########################################################################################################################

# TIME POINT 3:

# Create DESeqDataSet
dds_3 <- DESeqDataSetFromMatrix(countData = table_3_filtered,
                                colData = coldata_3,
                                design = ~ Condition)

# Set control group as reference
dds_3$Condition <- relevel(dds_3$Condition, ref = "control")

# Run DESeq
dds_out_3 <- DESeq(dds_3)

# View results
results_names_3 <- resultsNames(dds_out_3)
head(results_names_3)

# Iterate through results
for (comp in results_names_3){
  if (comp != 'Intercept') {
    # Get results
    res_deseq_3 <- results(dds_out_3, name = comp, alpha = 0.05)
    resLFC_3 <- lfcShrink(dds_out_3, coef = comp, res = res_deseq_3)
    
    # Optional: save or analyze results
    print(paste("Results for comparison:", comp))
    print(summary(resLFC_3))
  }
}

# Export results to CSV
write.csv(as.data.frame(resLFC_3), file = paste0("deseq3_results_", comp, "_T3.csv"))

# Consistent color scheme for conditions
condition_colors <- c("control" = "#00FF00", "infected" = "#FF69B4") 


#############################################
# PCA (Principal Component Analysis)        #
#############################################

# Obtain PCA data
vsd_3 <- varianceStabilizingTransformation(dds_out_3, blind = FALSE)  # Normalization while preserving experimental design
pca_data_3 <- plotPCA(vsd_3, intgroup = "Condition", returnData = TRUE)  # Get data for ggplot
percentVar_3 <- round(100 * attr(pca_data_3, "percentVar"))  # Percentage of explained variance

# Plot using ggplot2
library(ggplot2)
pca_plot <- ggplot(pca_data_3, aes(PC1, PC2, color = Condition)) + 
  geom_point(size = 3) +
  geom_text(aes(label = name), vjust = -0.5, size = 3) +  # Add labels
  scale_color_manual(values = condition_colors) +  # Apply manual color scheme
  xlab(paste0("PC1: ", percentVar_3[1], "% variance")) + 
  ylab(paste0("PC2: ", percentVar_3[2], "% variance")) + 
  theme_minimal() +
  ggtitle("PCA for Time Point 3")

# Display the plot
print(pca_plot)

# Save the plot with white background
ggsave("PCA_plot_3.png", plot = pca_plot, width = 8, height = 6, dpi = 300, bg = "white")



##########
# 3D PCA #
##########

# Obtain PCA data
vsd_3 <- varianceStabilizingTransformation(dds_out_3, blind = FALSE)  # Normalization
pca_res_3 <- prcomp(t(assay(vsd_3)))  # Compute PCA on transposed data

# Create a data frame with the first three principal components
pca_data_3d_3 <- as.data.frame(pca_res_3$x[, 1:3])  # Select PC1, PC2, and PC3
pca_data_3d_3$Condition <- colData(vsd_3)$Condition  # Add conditions as a column
colnames(pca_data_3d_3)[1:3] <- c("PC1", "PC2", "PC3")  # Rename columns

# Explained variance percentages
percentVar_3d_3 <- round(100 * pca_res_3$sdev^2 / sum(pca_res_3$sdev^2), 1)

# Create 3D PCA plot
pca_3d_plot_3 <- plot_ly(pca_data_3d_3, 
                         x = ~PC1, 
                         y = ~PC2, 
                         z = ~PC3, 
                         color = ~Condition, 
                         colors = c("control" = "#00FF00", "infected" = "#FF69B4"),
                         type = "scatter3d", 
                         mode = "markers+text",  # Show points and labels
                         marker = list(size = 5),
                         text = rownames(pca_data_3d_3),  # Use row names as labels
                         textposition = "top center",  # Text position
                         hoverinfo = "text+x+y+z",  # Hover information
                         hovertext = ~paste("Name:", rownames(pca_data_3d_3),
                                            "<br>Condition:", Condition,
                                            "<br>PC1:", round(PC1, 2), 
                                            "<br>PC2:", round(PC2, 2),
                                            "<br>PC3:", round(PC3, 2))) %>%
  layout(scene = list(
    xaxis = list(title = "PC1"),
    yaxis = list(title = "PC2"),
    zaxis = list(title = "PC3")),
    title = "3D PCA with Labels")

pca_3d_plot_3

htmlwidgets::saveWidget(as_widget(pca_3d_plot_3), "PCA_3D_TimePoint_3_labels.html")

################
# VOLCANO PLOT #
################

# Define thresholds
log2fc_threshold <- 1  # Change this if necessary
padj_threshold <- 0.05

# Classify genes
resLFC_3$diffexpressed <- "Not differentially expressed"
resLFC_3$diffexpressed[resLFC_3$log2FoldChange > log2fc_threshold & resLFC_3$padj < padj_threshold] <- "Up-regulated"
resLFC_3$diffexpressed[resLFC_3$log2FoldChange < -log2fc_threshold & resLFC_3$padj < padj_threshold] <- "Down-regulated"


# INDEPENDENT PLOT

# Calculate the maximum absolute value of log2FoldChange
max_abs_value_3 <- max(abs(resLFC_3$log2FoldChange), na.rm = TRUE)

# Check if data exists for each category
include_not_differentially_expressed <- any(resLFC_3$diffexpressed == "Not differentially expressed")
include_up_regulated <- any(resLFC_3$diffexpressed == "Up-regulated")
include_down_regulated <- any(resLFC_3$diffexpressed == "Down-regulated")

# Dynamically create breaks and labels
breaks_list <- c()
labels_list <- c()

# Add "Not differentially expressed" if present
if (include_not_differentially_expressed) {
  breaks_list <- c(breaks_list, "Not differentially expressed")
  labels_list <- c(labels_list, "Not differentially expressed")
}

# Add "Up-regulated" if present
if (include_up_regulated) {
  breaks_list <- c(breaks_list, "Up-regulated")
  labels_list <- c(labels_list, "Up-regulated")
}

# Add "Down-regulated" if present
if (include_down_regulated) {
  breaks_list <- c(breaks_list, "Down-regulated")
  labels_list <- c(labels_list, "Down-regulated")
}

# Create the plot
volcano_3 <- ggplot(data = as.data.frame(resLFC_3), aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed)) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0, col = "grey") +
  geom_hline(yintercept = -log10(padj_threshold), col = "grey") +
  scale_color_manual(
    values = c("Not differentially expressed" = "snow3",
               "Up-regulated" = "#FF6F61",
               "Down-regulated" = "#6EC5E9"),
    name = "Differential Expression",
    breaks = breaks_list,
    labels = labels_list
  ) +
  labs(x = "Log2 Fold Change", 
       y = "-Log10 (Adjusted p-value)", 
       title = "Volcano Plot for Time Point 3") +
  scale_x_continuous(limits = c(-max_abs_value_3 - 0.2, max_abs_value_3 + 0.2)) +
  guides(color = guide_legend(title = NULL)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_bw()

# Display the plot
print(volcano_3)

# Save the plot as PNG file
ggsave("volcano_plot_3.png", plot = volcano_3, width = 8, height = 6, dpi = 300, bg = "white")


# COMPARISON BETWEEN ALL 3 TIME POINTS (SAME AXES)

# Calculate the maximum absolute value of log2FoldChange
max_log2FoldChange_all <- max(c(abs(resLFC_1$log2FoldChange), abs(resLFC_2$log2FoldChange), abs(resLFC_3$log2FoldChange)))
max_neg_log10_padj_all <- max(c(-log10(resLFC_1$padj), -log10(resLFC_2$padj), -log10(resLFC_3$padj)), na.rm = TRUE)  # Ignore NA values

# Create the adjusted volcano plot
volcano_3_adjusted <- ggplot(data = as.data.frame(resLFC_3), aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed)) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0, col = "grey") +
  geom_hline(yintercept = -log10(padj_threshold), col = "grey") +
  scale_color_manual(
    values = c("Not differentially expressed" = "snow3",
               "Up-regulated" = "#FF6F61",
               "Down-regulated" = "#6EC5E9"),
    name = "Differential Expression",
    breaks = breaks_list,
    labels = labels_list
  ) +
  labs(x = "Log2 Fold Change", 
       y = "-Log10 (Adjusted p-value)", 
       title = "Adjusted Volcano Plot for Time Point 3") +
  scale_x_continuous(limits = c(-max_log2FoldChange_all - 0.2, max_log2FoldChange_all + 0.2)) +
  scale_y_continuous(limits = c(0, max_neg_log10_padj_all + 1)) +
  guides(color = guide_legend(title = NULL)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_bw()

# Display the adjusted plot
print(volcano_3_adjusted)

# Save the adjusted plot as PNG file
ggsave("volcano_plot_3_adjusted.png", plot = volcano_3_adjusted, width = 8, height = 6, dpi = 300, bg = "white")


###########################################################################################

# COMBINE THE 3 TIME POINTS INTO A VOLCANO PLOT

############################################
# 3 plots side by side with same legend #
############################################

# Make sure to add a 'Time' column to distinguish time points
resLFC_1$Time <- "Time 1"
resLFC_2$Time <- "Time 2"
resLFC_3$Time <- "Time 3"

# Define thresholds
log2fc_threshold <- 1  # Change if necessary
padj_threshold <- 0.05

# Classify genes for each dataset
for (resLFC in list(resLFC_1, resLFC_2, resLFC_3)) {
  resLFC$diffexpressed <- "Not differentially expressed"
  resLFC$diffexpressed[resLFC$log2FoldChange > log2fc_threshold & resLFC$padj < padj_threshold] <- "Up-regulated"
  resLFC$diffexpressed[resLFC$log2FoldChange < -log2fc_threshold & resLFC$padj < padj_threshold] <- "Down-regulated"
}

# Combine all datasets into a single one
combined_data <- rbind(resLFC_1, resLFC_2, resLFC_3)

# Calculate the maximum values for X and Y axes
max_log2FoldChange_all <- max(c(abs(resLFC_1$log2FoldChange), abs(resLFC_2$log2FoldChange), abs(resLFC_3$log2FoldChange)))
max_neg_log10_padj_all <- max(c(-log10(resLFC_1$padj), -log10(resLFC_2$padj), -log10(resLFC_3$padj)), na.rm = TRUE)

# Create the combined plot (side by side)
volcano_combined <- ggplot(combined_data, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed)) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0, col = "grey") +  # Vertical line at x = 0
  geom_hline(yintercept = -log10(padj_threshold), col = "grey") +  # Horizontal line for adjusted p-value
  scale_color_manual(values = c("Not differentially expressed" = "snow3", 
                                "Up-regulated" = "#FF6F61", 
                                "Down-regulated" = "#6EC5E9"),
                     name = "Differential Expression",
                     labels = c("Not differentially expressed", "Up-regulated", "Down-regulated")) +
  facet_wrap(~ Time) +  # Split by time point
  labs(x = "Log2 Fold Change", 
       y = "-Log10 (Adjusted p-value)", 
       title = "Combined Volcano Plot") +
  scale_x_continuous(limits = c(-max_log2FoldChange_all - 0.2, max_log2FoldChange_all + 0.2)) +  # Limit X axis
  scale_y_continuous(limits = c(0, max_neg_log10_padj_all + 1)) +  # Limit Y axis
  guides(color = guide_legend(title = NULL)) +  # Remove legend title
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate X-axis labels
  theme_bw()  # Use white theme

# Show the plot
print(volcano_combined)

# Save the combined plot as PNG file
ggsave("volcano_plot_combined.png", plot = volcano_combined, width = 12, height = 8, dpi = 300, bg = "white")



###############
# Superposed Plot #
###############

library(ggplot2)

# Make sure to add a 'Time' column to distinguish time points
resLFC_1$Time <- "Time 1"
resLFC_2$Time <- "Time 2"
resLFC_3$Time <- "Time 3"

# Define thresholds
log2fc_threshold <- 1  # Change if necessary
padj_threshold <- 0.05

# Classify genes for each dataset
for (resLFC in list(resLFC_1, resLFC_2, resLFC_3)) {
  resLFC$diffexpressed <- "Not differentially expressed"
  resLFC$diffexpressed[resLFC$log2FoldChange > log2fc_threshold & resLFC$padj < padj_threshold] <- "Up-regulated"
  resLFC$diffexpressed[resLFC$log2FoldChange < -log2fc_threshold & resLFC$padj < padj_threshold] <- "Down-regulated"
}

# Combine all datasets into a single one
combined_data <- rbind(resLFC_1, resLFC_2, resLFC_3)

# Calculate the maximum values for X and Y axes
max_log2FoldChange_all <- max(c(abs(resLFC_1$log2FoldChange), abs(resLFC_2$log2FoldChange), abs(resLFC_3$log2FoldChange)))
max_neg_log10_padj_all <- max(c(-log10(resLFC_1$padj), -log10(resLFC_2$padj), -log10(resLFC_3$padj)), na.rm = TRUE)

# Create the combined superposed plot
volcano_combined <- ggplot(combined_data, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed, shape = Time)) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0, col = "grey") +  # Vertical line at x = 0
  geom_hline(yintercept = -log10(padj_threshold), col = "grey") +  # Horizontal line for adjusted p-value
  scale_color_manual(values = c("Not differentially expressed" = "snow3", 
                                "Up-regulated" = "#FF6F61", 
                                "Down-regulated" = "#6EC5E9"),
                     name = "Differential Expression",
                     labels = c("Not differentially expressed", "Up-regulated", "Down-regulated")) +
  scale_shape_manual(values = c(16, 17, 18)) +  # Different shapes for time points
  labs(x = "Log2 Fold Change", 
       y = "-Log10 (Adjusted p-value)", 
       title = "Combined Volcano Plot (Superposed)") +
  scale_x_continuous(limits = c(-max_log2FoldChange_all - 0.2, max_log2FoldChange_all + 0.2)) +
  scale_y_continuous(limits = c(0, max_neg_log10_padj_all + 1)) +
  guides(color = guide_legend(title = "Expression Status"), 
         shape = guide_legend(title = "Time")) +  # Legends for expression and time
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme_bw()

# Show the plot
print(volcano_combined)

# Save the superposed plot as PNG file
ggsave("volcano_plot_combined_superposed.png", plot = volcano_combined, width = 12, height = 8, dpi = 300, bg = "white")
