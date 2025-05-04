########-------------------DOWNLOAD DATASET FROM ENCODE PORTAL --------------------------------########
## This should be done in the RStduio Terminal/ Linux Terminal 

# HSC
# curl -O -L https://www.encodeproject.org/files/ENCFF247FEJ/@@download/ENCFF247FEJ.tsv
# curl -O -L https://www.encodeproject.org/files/ENCFF064MKY/@@download/ENCFF064MKY.tsv


# CMP 
# curl -O -L https://www.encodeproject.org/files/ENCFF623OLU/@@download/ENCFF623OLU.tsv
# curl -O -L https://www.encodeproject.org/files/ENCFF691MHW/@@download/ENCFF691MHW.tsv

 
# CFU-E (ScriptSeq)
# curl -O -L https://www.encodeproject.org/files/ENCFF667IDY/@@download/ENCFF667IDY.tsv
# curl -O -L https://www.encodeproject.org/files/ENCFF655LMK/@@download/ENCFF655LMK.tsv


# Erythroblast (ScriptSeq)
# curl -O -L https://www.encodeproject.org/files/ENCFF342WUL/@@download/ENCFF342WUL.tsv
# curl -O -L https://www.encodeproject.org/files/ENCFF858JHF/@@download/ENCFF858JHF.tsv



#########-------------------INSTALL REQUIRED LIBRARIES ----------------------------------------########

library(ggrepel)
library(ggplot2)
library(apeglm)
library(dplyr)
library(readr)
library(tximport)
library(DESeq2)
library(SummarizedExperiment)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(tibble)
library(pheatmap)
library(RColorBrewer)


# Create a vector of file paths

file_All_Cells <- c("./Dataste_RNA_ScriptSeq/HSC/ENCFF247FEJ.tsv", # HSC
                    "./Dataste_RNA_ScriptSeq/HSC/ENCFF064MKY.tsv", # HSC
                    "./Dataste_RNA_ScriptSeq/CMP/ENCFF623OLU.tsv", # CMP
                    "./Dataste_RNA_ScriptSeq/CMP/ENCFF691MHW.tsv", # CMP
                    "./Dataste_RNA_ScriptSeq/CFU-E/ENCFF667IDY.tsv", # CFU-E
                    "./Dataste_RNA_ScriptSeq/CFU-E/ENCFF655LMK.tsv", # CFU- E
                    "./Dataste_RNA_ScriptSeq/ERY/ENCFF342WUL.tsv", # ERY
                    "./Dataste_RNA_ScriptSeq/ERY/ENCFF858JHF.tsv") # ERY

names(file_All_Cells) <- c("ENCFF247FEJ", "ENCFF064MKY", "ENCFF623OLU", "ENCFF691MHW",
                           "ENCFF667IDY", "ENCFF655LMK", "ENCFF342WUL", "ENCFF858JHF")


# Create a sample metadata data frame.
sample_info <- data.frame(
  sample = names(file_All_Cells),
  cell_type = c("HSC", "HSC",
                "CMP", "CMP",
                "CFU_E", "CFU_E",
                "ERY", "ERY")
)
rownames(sample_info) <- sample_info$sample

# Set HSC as the reference level
sample_info$cell_type <- relevel(factor(sample_info$cell_type), ref = "HSC")

# Import RSEM gene-level files using tximport
txi.rsem <- tximport(file_All_Cells, type = "rsem", txIn = FALSE, txOut = FALSE)
txi.rsem$length[txi.rsem$length == 0] <- 1

# Now create the DESeq2 dataset using the filtered object
ddsTxi <- DESeqDataSetFromTximport(txi.rsem, colData = sample_info, design = ~ cell_type)

# Run the DESeq2 analysis
dds <- DESeq(ddsTxi)

# view the counts
View(counts(dds))

# check normalized data
normalizationFactors(dds)

# Plot Dispersion Estimates
plotDispEsts(dds)

# check resultNames for dds
resultsNames(dds)

# Variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)

# Plot PCA plot
plotPCA(vsd, intgroup = "cell_type")

############################------------------PAIRWISE COMPARISONS------------------------------------------#################

##################################################
#--------------Compare HSC vs. CMP--------------
###################################################

# Step 1: get statistical results
res_HSC_vs_CMP <- results(dds, contrast = c("cell_type", "CMP", "HSC"), alpha = 0.05)
summary(res_HSC_vs_CMP)

# View(res_HSC_vs_CMP) by order of pvalue
head(res_HSC_vs_CMP[order(res_HSC_vs_CMP$pvalue), ])

## Save the unshrunken results to compare
res_HSC_CMP_unshrunken <- res_HSC_vs_CMP

resultsNames(dds)

# Step 2: get shrunken log2FCs for visualization
res_HSC_vs_CMP <- lfcShrink(dds, coef="cell_type_CMP_vs_HSC", type="apeglm")

# Step 3: plot the MA.compare shrunken vs unshrunken
plotMA(res_HSC_CMP_unshrunken, ylim=c(-2,2))
plotMA(res_HSC_vs_CMP, ylim=c(-2,2))


################################################
# 1. Differential Expression Analysis & Annotation
################################################

# Set padj cutoff
padj.cutoff <- 0.05

# Convert the DESeq2 results table to a tibble and then to a data.frame,
# moving the rownames (Ensembl IDs) into a column called "gene"
res_HSC_vs_CMP_tb <- as_tibble(as.data.frame(res_HSC_vs_CMP), rownames = "gene") %>% 
  as.data.frame()

# Set the rownames from the "gene" column and then remove the "gene" column
rownames(res_HSC_vs_CMP_tb) <- res_HSC_vs_CMP_tb$gene
res_HSC_vs_CMP_tb$gene <- NULL

# Create a helper function to annotate genes based on rownames (cleaning the Ensembl IDs)
annotate_genes <- function(df) {
  # Remove version numbers from the rownames (e.g., ENSMUSG00000000078.6 -> ENSMUSG00000000078)
  gene_ids <- gsub("\\..*", "", rownames(df))
  
  # Map cleaned Ensembl IDs to gene symbols
  df$symbol <- mapIds(org.Mm.eg.db,
                      keys = gene_ids,
                      column = "SYMBOL",
                      keytype = "ENSEMBL",
                      multiVals = "first")
  
  # Also add the cleaned gene IDs as a new column for later filtering
  df$gene_id <- gene_ids
  return(df)
}

# Annotate the DE results table
res_HSC_vs_CMP_tb <- annotate_genes(res_HSC_vs_CMP_tb)

# Subset the DE results to keep only significant genes (padj < padj.cutoff)
sig_HSC_vs_CMP <- res_HSC_vs_CMP_tb %>% 
  filter(padj < padj.cutoff)
summary(sig_HSC_vs_CMP)

# Save all significant genes for HSC vs CMP
write.csv(sig_HSC_vs_CMP, './Results/HSC_CMP_All_SigGenes.all.csv')

# Extract the subset of DEGs with |log2FoldChange| > 2 for further analyses
High_sig_HSC_vs_CMP <- res_HSC_vs_CMP_tb %>% 
  filter(padj < padj.cutoff & abs(log2FoldChange) > 2)

# Save all significant genes for HSC vs CMP with log2foldchange> 2
write.csv(High_sig_HSC_vs_CMP, './Results/HSC_CMP_High_SigGeneswithlog2foldchangegreaterthan2.all.csv')

length(High_sig_HSC_vs_CMP)

##############################
# 2. Normalized Counts Extraction & Annotation for Plotting
##############################

# Extract normalized counts from dds, convert to a data frame, and move rownames (Ensembl IDs) into a column
normalized_counts <- counts(dds, normalized = TRUE) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "ensembl")

# Create a column "gene_id" by removing version numbers from the ensembl column
normalized_counts <- normalized_counts %>% 
  mutate(gene_id = gsub("\\..*", "", ensembl))

# Annotate normalized counts: map cleaned gene IDs to gene symbols
normalized_counts$symbol <- mapIds(org.Mm.eg.db,
                                   keys = normalized_counts$gene_id,
                                   column = "SYMBOL",
                                   keytype = "ENSEMBL",
                                   multiVals = "first")

##############################
# 3. Filter Normalized Counts for Significant Genes
##############################

# Check if sig_HSC_vs_CMP has a "symbol" column; if not, annotate it first.
if (!("symbol" %in% colnames(sig_HSC_vs_CMP))) {
  message("Column 'symbol' not found in sig_HSC_vs_CMP. Annotating the DE results...")
  sig_HSC_vs_CMP <- annotate_genes(sig_HSC_vs_CMP)
}

# Check if a 'gene_id' column exists in sig_HSC_vs_CMP; if not, create it from the rownames.
if (!("gene_id" %in% colnames(sig_HSC_vs_CMP))) {
  sig_HSC_vs_CMP$gene_id <- gsub("\\..*", "", rownames(sig_HSC_vs_CMP))
}

# Prepare the significant DE results subset.
# Rename the "symbol" column from DE results to distinguish it from the normalized counts annotation.
sig_subset <- sig_HSC_vs_CMP %>% 
  dplyr::select(gene_id, symbol) %>% 
  dplyr::rename(symbol_sig = symbol)

# Now, from your normalized counts, select the identifier columns and the sample columns.
# This assumes normalized_counts already contains 'gene_id', 'symbol', and your sample columns.
norm_sig_HSC_vs_CMP <- normalized_counts %>% 
  dplyr::select(gene_id, symbol, ENCFF247FEJ, ENCFF064MKY, ENCFF623OLU, ENCFF691MHW) %>% 
  # Join the DE results subset by gene_id so that only significant genes are retained.
  inner_join(sig_subset, by = "gene_id") %>% 
  # Replace the gene symbol in normalized_counts with the one from the DE results.
  mutate(symbol = symbol_sig) %>% 
  dplyr::select(-symbol_sig)


#--------------------PLOTS HEATMAP AND VOLCANO PLOTS-------------------------------------------------

#-------------HEATMAP -----------------------------------
### Set a color palette
heat_colors <- brewer.pal(6, "YlOrRd")

### Run pheatmap using the metadata data frame for the annotation
pheatmap(norm_sig_HSC_vs_CMP[3:6], 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = sample_info, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20)


#----------------------------Volcano Plot----------------------------------

# Prepare the DE results table with a new column for -log10(padj)
res_table <- res_HSC_vs_CMP_tb %>%
  dplyr::mutate(neg_log10padj = -log10(padj)) %>% 
  # Define a threshold flag for consistency (this column is optional here)
  dplyr::mutate(threshold_OE = (padj < 0.05 & abs(log2FoldChange) > 2))

# 3) Identify top genes among those passing the threshold
#    For example, take the top 10 "up" and top 10 "down" by padj.
top_genes <- res_table %>%
  dplyr::filter(threshold_OE, abs(log2FoldChange) > 0) %>%
  dplyr::arrange(padj) %>%
  dplyr::slice(1:50)

#top_down <- res_table %>%
  #dplyr::filter(threshold_OE, log2FoldChange < 0) %>%
  #dplyr::arrange(padj) %>%
  #dplyr::slice(1:15)

# Combine them; these will be the labeled points
#top_genes <- bind_rows(top_up, top_down)

# 4) Volcano plot
volcano_plot <- ggplot(res_table, aes(x = log2FoldChange, y = neg_log10padj)) +
  # Plot all genes
  geom_point(aes(color = threshold_OE), alpha = 0.6, size = 2) +
  # Color scale: threshold-passing genes in teal, others in gray
  scale_color_manual(values = c("FALSE" = "gray70", "TRUE" = "red")) +
  # Label only the top genes among those passing threshold
  geom_text_repel(
    data = top_genes,
    aes(label = symbol),
    box.padding = 0.5,
    point.padding = 0.3,
    max.overlaps = Inf,
    size = 4
  ) +
  # Titles and theme
  ggtitle("Vocalno Plot for Top 50 genes") +
  xlab("log2FoldChange") +
  ylab("-log10 adjusted p-value") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

# 5) Print the plot
print(volcano_plot)




##---------------------------- Compare HSC vs. ERY-------------------------------------------

##################################################
#--------------Compare HSC vs. CMP--------------
###################################################

# Step 1: get statistical results
res_HSC_vs_ERY <- results(dds, contrast = c("cell_type", "ERY", "HSC"), alpha = 0.05)
summary(res_HSC_vs_ERY)

# View(res_HSC_vs_CMP) by order of pvalue
head(res_HSC_vs_ERY[order(res_HSC_vs_ERY$pvalue), ])

## Save the unshrunken results to compare
res_HSC_ERY_unshrunken <- res_HSC_vs_ERY

# Step 2: get shrunken log2FCs for visualization
res_HSC_vs_ERY <- lfcShrink(dds, coef = "cell_type_ERY_vs_HSC", res = res_HSC_vs_ERY)

# Step 3: plot the MA.compare shrunken vs unshrunken
plotMA(res_HSC_ERY_unshrunken, ylim=c(-2,2))
plotMA(res_HSC_vs_ERY, ylim=c(-2,2))


################################################
# 1. Differential Expression Analysis & Annotation
################################################

# Set padj cutoff
padj.cutoff <- 0.05

# Convert the DESeq2 results table to a tibble and then to a data.frame,
# moving the rownames (Ensembl IDs) into a column called "gene"
res_HSC_vs_ERY_tb <- as_tibble(as.data.frame(res_HSC_vs_ERY), rownames = "gene") %>% 
  as.data.frame()

# Set the rownames from the "gene" column and then remove the "gene" column
rownames(res_HSC_vs_ERY_tb) <- res_HSC_vs_ERY_tb$gene
res_HSC_vs_ERY_tb$gene <- NULL

# Create a helper function to annotate genes based on rownames (cleaning the Ensembl IDs)
annotate_genes <- function(df) {
  # Remove version numbers from the rownames (e.g., ENSMUSG00000000078.6 -> ENSMUSG00000000078)
  gene_ids <- gsub("\\..*", "", rownames(df))
  
  # Map cleaned Ensembl IDs to gene symbols
  df$symbol <- mapIds(org.Mm.eg.db,
                      keys = gene_ids,
                      column = "SYMBOL",
                      keytype = "ENSEMBL",
                      multiVals = "first")
  
  # Also add the cleaned gene IDs as a new column for later filtering
  df$gene_id <- gene_ids
  return(df)
}

# Annotate the DE results table
res_HSC_vs_ERY_tb <- annotate_genes(res_HSC_vs_ERY_tb)

# Subset the DE results to keep only significant genes (padj < padj.cutoff)
sig_HSC_vs_ERY <- res_HSC_vs_ERY_tb %>% 
  filter(padj < padj.cutoff)

# Save all significant genes for HSC vs CMP
write.csv(sig_HSC_vs_ERY, './Results/HSC_ERY_All_SigGenes.all.csv')

# Extract the subset of DEGs with |log2FoldChange| > 8 for further analyses
High_sig_HSC_vs_ERY <- res_HSC_vs_ERY_tb %>% 
  filter(padj < padj.cutoff & abs(log2FoldChange) > 8)

# Save all significant genes for HSC vs CMP with log2foldchange> 2
write.csv(High_sig_HSC_vs_ERY, './Results/HSC_ERY_High_SigGeneswithlog2foldchangegreaterthan8.all.csv')

##############################
# 2. Normalized Counts Extraction & Annotation for Plotting
##############################

# Extract normalized counts from dds, convert to a data frame, and move rownames (Ensembl IDs) into a column
normalized_counts <- counts(dds, normalized = TRUE) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "ensembl")

# Create a column "gene_id" by removing version numbers from the ensembl column
normalized_counts <- normalized_counts %>% 
  mutate(gene_id = gsub("\\..*", "", ensembl))

# Annotate normalized counts: map cleaned gene IDs to gene symbols
normalized_counts$symbol <- mapIds(org.Mm.eg.db,
                                   keys = normalized_counts$gene_id,
                                   column = "SYMBOL",
                                   keytype = "ENSEMBL",
                                   multiVals = "first")

##############################
# 3. Filter Normalized Counts for Significant Genes
##############################

# Check if sig_HSC_vs_ERY has a "symbol" column; if not, annotate it first.
if (!("symbol" %in% colnames(sig_HSC_vs_ERY))) {
  message("Column 'symbol' not found in sig_HSC_vs_ERY. Annotating the DE results...")
  sig_HSC_vs_ERY <- annotate_genes(sig_HSC_vs_ERY)
}

# Check if a 'gene_id' column exists in sig_HSC_vs_CMP; if not, create it from the rownames.
if (!("gene_id" %in% colnames(sig_HSC_vs_ERY))) {
  sig_HSC_vs_ERY$gene_id <- gsub("\\..*", "", rownames(sig_HSC_vs_ERY))
}

# Prepare the significant DE results subset.
# Rename the "symbol" column from DE results to distinguish it from the normalized counts annotation.
sig_subset_ERY <- sig_HSC_vs_ERY %>% 
  dplyr::select(gene_id, symbol) %>% 
  dplyr::rename(symbol_sig = symbol)

# Now, from your normalized counts, select the identifier columns and the sample columns.
# This assumes normalized_counts already contains 'gene_id', 'symbol', and your sample columns.
norm_sig_HSC_vs_ERY <- normalized_counts %>% 
  dplyr::select(gene_id, symbol, ENCFF247FEJ, ENCFF064MKY, ENCFF342WUL, ENCFF858JHF) %>% 
  # Join the DE results subset by gene_id so that only significant genes are retained.
  inner_join(sig_subset_ERY, by = "gene_id") %>% 
  # Replace the gene symbol in normalized_counts with the one from the DE results.
  mutate(symbol = symbol_sig) %>% 
  dplyr::select(-symbol_sig)


#--------------------PLOTS HEATMAP AND VOLCANO PLOTS-------------------------------------------------

#-------------HEATMAP -----------------------------------
### Set a color palette
heat_colors <- brewer.pal(6, "YlOrRd")

### Run pheatmap using the metadata data frame for the annotation
pheatmap(norm_sig_HSC_vs_ERY[3:6], 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = sample_info, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20)


#----------------------------Volcano Plot----------------------------------

# Prepare the DE results table with a new column for -log10(padj)
res_table_ERY <- res_HSC_vs_ERY_tb %>%
  dplyr::mutate(neg_log10padj = -log10(padj)) %>% 
  # Define a threshold flag for consistency (this column is optional here)
  dplyr::mutate(threshold_OE = (padj < 0.05 & abs(log2FoldChange) > 8))

# 3) Identify top genes among those passing the threshold
#    For example, take the top 10 "up" and top 10 "down" by padj.
top_genes_ERY <- res_table_ERY %>%
  dplyr::filter(threshold_OE, abs(log2FoldChange) > 0) %>%
  dplyr::arrange(padj) %>%
  dplyr::slice(1:50)

#top_down_ERY <- res_table_ERY %>%
 # dplyr::filter(threshold_OE, log2FoldChange < 0) %>%
  #dplyr::arrange(padj) %>%
  #dplyr::slice(1:10)

# Combine them; these will be the labeled points
#top_genes_ERY <- bind_rows(top_up_ERY, top_down_ERY)

# 4) Volcano plot
volcano_plot_ERY <- ggplot(res_table_ERY, aes(x = log2FoldChange, y = neg_log10padj)) +
  # Plot all genes
  geom_point(aes(color = threshold_OE), alpha = 0.6, size = 2) +
  # Color scale: threshold-passing genes in teal, others in gray
  scale_color_manual(values = c("FALSE" = "gray70", "TRUE" = "red")) +
  # Label only the top genes among those passing threshold
  geom_text_repel(
    data = top_genes_ERY,
    aes(label = symbol),
    box.padding = 0.5,
    point.padding = 0.4,
    max.overlaps = 50,
    size = 4
  ) +
  # Titles and theme
  ggtitle(" Volcano Plot for Top 50 genes") +
  xlab("log2FoldChange") +
  ylab("-log10 adjusted p-value") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

# 5) Print the plot
print(volcano_plot_ERY)


# 1. Use variance-stabilizing transformation (recommended for clustering)
vsd <- vst(dds, blind = FALSE)

# 2. Extract the expression matrix
expr_mat <- assay(vsd)  # rows = genes, columns = samples

# 3. Optionally subset to highly variable genes (for cleaner clustering)
top_var_genes <- head(order(rowVars(expr_mat), decreasing = TRUE), 5000)
expr_top <- expr_mat[top_var_genes, ]

# 4. Plot hierarchical clustering heatmap **of samples**
pheatmap(expr_top,
         scale = "row",
         color = heat_colors,
         show_rownames = FALSE,
         show_colnames = TRUE,
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         annotation_col = as.data.frame(colData(dds)["cell_type"]))


















































