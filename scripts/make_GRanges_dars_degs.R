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

# make Granges objects ----------------------------

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

# pull out hgnc_symbols we want to query
features <- degs_all$feature

# perform query
gene_coordinates <- getBM(
  attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
  filters = "hgnc_symbol",
  values = features,
  mart = ensembl
)

# many random chrs, let's keep only coordinates that are one standard chrs
standard_chromosomes <- as.character(c(1:22, "X", "Y"))
gene_coordinates_standard <- gene_coordinates[gene_coordinates$chromosome_name %in% standard_chromosomes, ] %>% distinct()

# combine as downloaded degs data with the genomic coordinates
degs_all_coords <- inner_join(degs_all, gene_coordinates_standard, by = c("feature" = "hgnc_symbol"))

# check for NAs, cannot have NAs in the df that is used to make GRanges
sum(is.na(degs_all_coords$start_position))

# make GRanges object with differentially expressed genes
degs_gr <- GRanges(seqnames = degs_all_coords$chromosome_name,
                   ranges = IRanges(start = degs_all_coords$start_position, end = degs_all_coords$end_position),
                   strand = "*")


