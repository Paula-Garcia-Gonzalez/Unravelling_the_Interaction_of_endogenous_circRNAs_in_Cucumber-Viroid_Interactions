# Load required packages
library(edgeR)

# Load count table
counts <- read.csv("final_counts_zeros.csv", row.names = 1)
head(counts)

# Split count table into 3 time points
table_1 <- counts[, grepl("^R00_1", colnames(counts))]  # Select columns starting with "R00_1"
table_2 <- counts[, grepl("^R00_2", colnames(counts))]  # Select columns starting with "R00_2"
table_3 <- counts[, grepl("^R00_3", colnames(counts))]  # Select columns starting with "R00_3"

# Save each table into separate files
write.table(table_1, "table_1.csv", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(table_2, "table_2.csv", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(table_3, "table_3.csv", sep = "\t", quote = FALSE, row.names = TRUE)

# Load sample metadata
coldata_1 <- read.csv("coldata_1.csv", row.names = 1)
coldata_2 <- read.csv("coldata_2.csv", row.names = 1)
coldata_3 <- read.csv("coldata_3.csv", row.names = 1)

# Make sure Condition is a factor
coldata_1$Condition <- as.factor(coldata_1$Condition)
coldata_2$Condition <- as.factor(coldata_2$Condition)
coldata_3$Condition <- as.factor(coldata_3$Condition)

# Filtering function
filter_table <- function(table) {
  table[rowSums(table >= 5) >= 3, ]  # Keep genes with at least 3 samples with count >= 5
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

# Create DGEList object for filtered data
dge_1 <- DGEList(counts = table_1_filtered, group = coldata_1$Condition)

# Normalize the data
dge_1 <- calcNormFactors(dge_1)

# Create design matrix
design_1 <- model.matrix(~ coldata_1$Condition)
colnames(design_1) <- c("Intercept", "ConditionInfected")

# Estimate dispersion
dge_1 <- estimateDisp(dge_1, design_1)

# Fit GLM model
fit_1 <- glmFit(dge_1, design_1)

# Perform likelihood ratio test for condition effect
lrt_1 <- glmLRT(fit_1, coef = 2)

# Display top differentially expressed genes
topTags(lrt_1)

# Export full results
results_1 <- as.data.frame(topTags(lrt_1, n = Inf))
write.csv(results_1, "results_time_1.csv")

########################################################################################################################

# TIME POINT 2:

# Create DGEList object
dge_2 <- DGEList(counts = table_2_filtered, group = coldata_2$Condition)

# Normalize
dge_2 <- calcNormFactors(dge_2)

# Design matrix
design_2 <- model.matrix(~ coldata_2$Condition)
colnames(design_2) <- c("Intercept", "ConditionInfected")

# Dispersion estimation
dge_2 <- estimateDisp(dge_2, design_2)

# GLM fit
fit_2 <- glmFit(dge_2, design_2)

# Likelihood ratio test
lrt_2 <- glmLRT(fit_2, coef = 2)

# Top DE genes
topTags(lrt_2)

# Export results
results_2 <- as.data.frame(topTags(lrt_2, n = Inf))
write.csv(results_2, "results_time_2.csv")

########################################################################################################################

# TIME POINT 3:

# Create DGEList object
dge_3 <- DGEList(counts = table_3_filtered, group = coldata_3$Condition)

# Normalize
dge_3 <- calcNormFactors(dge_3)

# Design matrix
design_3 <- model.matrix(~ coldata_3$Condition)
colnames(design_3) <- c("Intercept", "ConditionInfected")

# Dispersion estimation
dge_3 <- estimateDisp(dge_3, design_3)

# GLM fit
fit_3 <- glmFit(dge_3, design_3)

# Likelihood ratio test
lrt_3 <- glmLRT(fit_3, coef = 2)

# Top DE genes
topTags(lrt_3)

# Export results
results_3 <- as.data.frame(topTags(lrt_3, n = Inf))
write.csv(results_3, "results_time_3.csv")
