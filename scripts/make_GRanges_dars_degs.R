suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(GenomicRanges)
  library(ggplot2)
  library(NatParksPalettes)
  library(cowplot)
  library(stringr)
  library(GenomicFeatures)
  library(ensembldb)
  library(EnsDb.Hsapiens.v86) 
})


## read in data 
dars_all <- read.csv('source_data/TableS2_DAR_between_AD_HC.txt',sep = '\t')
degs_all <- read.csv('source_data/TableS3_DEG_between_AD_HC.txt', sep = '\t')

# make DARs Granges object ----------------------------

# turn DARs into Granges object for subsequent annotation steps 

dars_sitename_splt <- str_split(dars_all$sitename, pattern = '-',simplify = T)
dars_chr <- dars_sitename_splt[,1]
dars_start <- as.numeric(dars_sitename_splt[,2])
dars_end <- as.numeric(dars_sitename_splt[,3])

dars_gr <- GRanges(seqnames = dars_chr, 
                   ranges = IRanges(start = dars_start, end = dars_end),
                   strand = "*")
dars_gr$cell_type <- dars_all$cell_type
dars_gr$logFC_dar <- dars_all$avg_log2FC

# make DEGs Granges object ------------------------------------------------

# in data as provided, DEGs are just gene_name, but no genomic coordinates
# annotate with genomic coordinates

features <- degs_all$feature

# Load EnsDb and filter gene coordinates
edb <- EnsDb.Hsapiens.v86
g <- genes(edb)
g <- keepStandardChromosomes(g, pruning.mode = 'coarse') # standard chromosomes

# Filter the genes that match with DEGs from downloaded data
degs_coords <- g[g$symbol %in% features] %>% as.data.frame()

# Merge DEGs data with gene coordinates so we can have coordinates for our genes
degs_coords_celltype <- merge(degs_all, degs_coords, by.x = 'feature', by.y = 'symbol')

# Filter mapped coordinates by checking that start and end values are not NA
degs_coords_mapped <- degs_coords_celltype[!is.na(degs_coords_celltype$start) & !is.na(degs_coords_celltype$end), ]

# Create GRanges object
degs_gr <- GRanges(
  seqnames = paste0('chr', degs_coords_mapped$seqnames),  # correct format for seqnames
  ranges = IRanges(start = degs_coords_mapped$start, end = degs_coords_mapped$end),
  strand = "*"  # placeholder bc strand not relevant 
)

# Add metadata from downloaded data 
degs_gr$cell_type <- degs_coords_celltype$cell_type
degs_gr$gene_name <- degs_coords_celltype$feature
degs_gr$gene_biotype <- degs_coords_celltype$gene_biotype
degs_gr$DAR_deg_overlap <- degs_coords_celltype$DAR_DEG_overlap
degs_gr$log2FC_deg <- degs_coords_celltype$avg_log2FC

degs_gr$cell_type <- sapply(degs_gr$cell_type, function(x) {
  if (x %in% c("B_intermediate", "B_memory", "B_naive", "Plasmablast")) {
    "B_Cells"
  } else if (x %in% c("CD14_Mono", "CD16_Mono")) {
    "Monocytes"
  } else if (x %in% c("ASDC", "cDC1", "cDC2", "pDC")) {
    "Dendritic_Cells"
  } else if (x %in% c("CD4_CTL", "CD4_Naive", "CD4_Proliferating", "CD4_TCM", "CD4_TEM", "Treg")) {
    "CD4+_T_Cells"
  } else if (x %in% c("CD8_Naive", "CD8_Proliferating", "CD8_TCM", "CD8_TEM", "MAIT")) {
    "CD8+_T_Cells"
  } else if (x %in% c("NK", "NK_Proliferating", "NK_CD56bright")) {
    "NK_Cells"
  } else if (x %in% c("ILC", "dnT", "gdT")) {
    "Other_T_Cells"
  } else if (x %in% c("Platelet", "Eryth", "HSPC", "Doublet")) {
    "Other"
  } else {
    x # Return original value if no match is found
  }
})

