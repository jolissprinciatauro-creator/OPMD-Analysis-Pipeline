
# Author: Joliss Princia Tauro
# Description: End-to-end analysis of Oral Premalignant Disorders (OPMD) 
# including QC, DGE, and Functional Enrichment.

# 1. Load Libraries
library(GEOquery)
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(ComplexHeatmap)
library(enrichplot)

# 2. Data Retrieval
gse_id <- "GSE10174" 
gse <- getGEO(gse_id, destdir = ".", getGPL = TRUE)
expr_data <- exprs(gse[[1]])

# 3. Quality Control: Principal Component Analysis (PCA)
# Helps visualize if Control vs Disease samples cluster correctly
pca_res <- prcomp(t(expr_data), scale. = TRUE)
pca_df <- as.data.frame(pca_res$x)
pca_df$Group <- c(rep("Control", 3), rep("Disease", 3)) # Adjust based on metadata

ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4) +
  theme_minimal() +
  labs(title = "PCA: Sample Clustering for GSE10174", x = "PC1", y = "PC2")

# 4. Differential Expression Analysis (DEA)
design <- model.matrix(~0 + factor(c(rep("Control", 3), rep("Disease", 3))))
colnames(design) <- c("Control", "Disease")

fit <- lmFit(expr_data, design)
cont_matrix <- makeContrasts(Disease-Control, levels=design)
fit2 <- contrasts.fit(fit, cont_matrix)
fit2 <- eBayes(fit2)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)

# 5. Visualization: Volcano Plot
tT$diffexpressed <- "NO"
tT$diffexpressed[tT$logFC > 1 & tT$adj.P.Val < 0.05] <- "UP"
tT$diffexpressed[tT$logFC < -1 & tT$adj.P.Val < 0.05] <- "DOWN"

ggplot(tT, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed)) +
  geom_point(alpha=0.4) +
  theme_minimal() +
  scale_color_manual(values=c("blue", "grey", "red")) +
  labs(title="Differential Gene Expression in OPMD", x="log2 Fold Change", y="-log10 Adjusted P-Value")

# 6. Visualization: Heatmap of Top 50 Genes
top_genes <- head(rownames(tT), 50)
matrix_top <- expr_data[top_genes, ]
# Z-score normalization for the heatmap
matrix_scaled <- t(scale(t(matrix_top)))

Heatmap(matrix_scaled, name = "Expression", 
        column_title = "Samples", 
        row_title = "Top 50 DEGs",
        show_row_names = TRUE)

# 7. Pathway Enrichment (GO)
# Filter for significant genes
sig_genes <- rownames(tT[tT$diffexpressed != "NO", ])

go_results <- enrichGO(gene = sig_genes, 
                       OrgDb = org.Hs.eg.db, 
                       keyType = "SYMBOL", 
                       ont = "BP", 
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05)

# 8. Advanced Visualization: Gene-Concept Network (Cnetplot)
# This shows how genes are shared across the top 5 biological pathways
cnetplot(go_results, categorySize="pvalue", showCategory = 5)

# 9. Visualization: Dotplot of Pathways
dotplot(go_results, showCategory=15) + labs(title="Top Biological Processes in OPMD")