suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(GenomicRanges)
  library(ggplot2)
  library(NatParksPalettes)
  library(cowplot)
  library(knitr)
})


## read in data 
dars_all <- read.csv('source_data/TableS2_DAR_between_AD_HC.txt',sep = '\t')
degs_all <- read.csv('source_data/TableS3_DEG_between_AD_HC.txt', sep = '\t')

## filter data by p-value and log2fc cut offs for significant data 
dars_all <- dars_all %>% 
  filter(LR_padj < 0.05) %>% 
  filter(avg_log2FC > 0.125 | avg_log2FC < -0.125)


degs_all <- degs_all %>% 
  filter(MAST_padj < 0.05) %>% 
  filter(avg_log2FC > 0.125 | avg_log2FC < -0.125)


# Ensure the cell_type column in degs_gr matches the groupings
degs_all$cell_type <- sapply(degs_all$cell_type, function(x) {
  if (x %in% c("B_intermediate", "B_memory", "B_naive", "Plasmablast")) {
    return("B_Cells")
  } else if (x %in% c("CD14_Mono", "CD16_Mono")) {
    return("Monocytes")
  } else if (x %in% c("ASDC", "cDC1", "cDC2", "pDC")) {
    return("Dendritic_Cells")
  } else if (x %in% c("CD4_CTL", "CD4_Naive", "CD4_Proliferating", "CD4_TCM", "CD4_TEM", "Treg")) {
    return("CD4+_T_Cells")
  } else if (x %in% c("CD8_Naive", "CD8_Proliferating", "CD8_TCM", "CD8_TEM", "MAIT")) {
    return("CD8+_T_Cells")
  } else if (x %in% c("NK", "NK_Proliferating", "NK_CD56bright")) {
    return("NK_Cells")
  } else if (x %in% c("ILC", "dnT", "gdT")) {
    return("Other_T_Cells")
  } else if (x %in% c("Platelet", "Eryth", "HSPC", "Doublet")) {
    return("Other")
  } else {
    return(x) # Return the original value if no match is found
  }
})




## do a few quick summary visualizations
dars_plt <- ggplot(dars_all, aes(x = cell_type, y = avg_log2FC)) + 
  geom_violin(aes(fill = cell_type)) +
  geom_boxplot(width = 0.1) +
  ggtitle('Log2FC Chromatin Accessibility Based on Cell Type') +
  scale_fill_manual(values = natparks.pals('Yellowstone',6)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

peaktype_plt <- ggplot(dars_all %>% group_by(cell_type, peakType) %>% tally(), aes(x = cell_type,y = n)) + 
  geom_col(aes(fill=peakType)) +
  scale_fill_manual(values = natparks.pals('DeathValley',6)) +
  ggtitle('Peak Annotation Based on Cell Type') +
  ylab('Number of Peaks') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

degs_plt <- ggplot(degs_all, aes(x = cell_type, y = avg_log2FC)) + 
  geom_violin(aes(fill = cell_type)) +
  geom_boxplot(width = 0.1) +
  scale_fill_manual(values = natparks.pals('Yellowstone',n=17)) +
  ggtitle('Log2FC Gene Expression Based on Cell Type') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

summary_overlaps <- degs_all %>% group_by(cell_type,DAR_DEG_overlap) %>% tally()

summary_overlaps %>%
  kable(caption = "Counts of DAR-DEG Overlaps by Cell Type",
    col.names = c("Cell Type", "DAR-DEG Overlap", "Count"))

dars_summary <- plot_grid(dars_plt, peaktype_plt, ncol = 2)
plts_summary <- plot_grid(dars_plt, degs_plt, ncol = 1)

#pca to see any clustering with cell type or peak type

pca <- prcomp(dars_all[, c("avg_log2FC", "pct.1", "pct.2")], scale. = TRUE)
pca_df <- data.frame(pca$x)
pca_df$cell_type <- dars_all$cell_type
pca_df$peak_type <- dars_all$peakType

ggplot(pca_df, aes(x = PC1, y = PC2, color = cell_type)) +
  geom_jitter() +
  labs(title = "PCA Plot by Cell Type") #possibly some clustering
#can compare to the log2fc chromatin accessibility based on Cell type violin plot

ggplot(pca_df, aes(x = PC1, y = PC2, color = peak_type)) +
  geom_point() +
  labs(title = "PCA Plot by Peak Type") #no clustering




