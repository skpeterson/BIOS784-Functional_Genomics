suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(GenomicRanges)
  library(ggplot2)
  library(NatParksPalettes)
  library(cowplot)
  library(stringr)
  library(GenomicFeatures)
  library(AnnotationHub)
  library(ChIPseeker)
})


## read in data 
dars_all <- read.csv('source_data/TableS2_DAR_between_AD_HC.txt',sep = '\t')
degs_all <- read.csv('source_data/TableS3_DEG_between_AD_HC.txt', sep = '\t')

# explore annotating DARs a few different ways ----------------------------

# turn DARs into Granges object for subsequent annotation steps 

dars_sitename_splt <- str_split(dars_all$sitename, pattern = '-',simplify = T)
dars_chr <- dars_sitename_splt[,1]
dars_start <- as.numeric(dars_sitename_splt[,2])
dars_end <- as.numeric(dars_sitename_splt[,3])

dars_gr <- GRanges(seqnames = dars_chr, 
                   ranges = IRanges(start = dars_start, end = dars_end),
                   strand = "*")
dars_gr$cell_type <- dars_all$cell_type



# turn DEGs into Granges object for subsequent annotation steps 

# in data as provided, DEGs are just gene_name, but no genomic coordinates
# annotate with genomic coordinates

# Connect to the Ensembl database
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# Create a vector of feature (gene) names that you want to query
features <- degs_all$feature

gene_coordinates <- getBM(
  attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
  filters = "hgnc_symbol",
  values = features,
  mart = ensembl
)


ah <- AnnotationHub()
orgs <- subset(ah, ah$rdataclass == "OrgDb")
orgdb <- query(orgs, "Homo sapiens")[[1]]


gene_info <- as.data.frame(select(orgdb, keys = keys(orgdb), columns = c("GENENAME", "GENETYPE", "ENSEMBL")))

protein_coding_genes <- gene_info %>%
  filter(GENETYPE == "protein_coding")

peak_annotation_coding <- annotatePeak(dars_gr, 
                                       TxDb = protein_coding_genes,
                                       level = 'gene')

