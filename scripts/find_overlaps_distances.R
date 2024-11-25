# Load necessary libraries
library(GenomicRanges)
library(plyranges)


# Define promoter regions for DEGs (e.g., 2 kb upstream and 500 bp downstream of the TSS)
degs_promoters <- promoters(degs_gr, upstream = 2000, downstream = 500)


# Step 1: Identify unique cell types
cell_types <- intersect(unique(degs_gr$cell_type), unique(dars_gr$cell_type))

# Step 2: Initialize a list to store results
overlap_hits <- list()

# Step 3: Loop over each cell type and compute overlaps
for (ct in cell_types) {
  # Subset DEGs promoters and DARs for the current cell type
  degs_subset <- degs_promoters[degs_promoters$cell_type == ct]
  dars_subset <- dars_gr[dars_gr$cell_type == ct]
  
  # Find overlaps with DARs as the first argument
  overlaps <- findOverlaps(dars_subset, degs_subset)
  
  # Store the Hits object
  overlap_hits[[ct]] <- overlaps
}

# Step 5: Calculate unique DARs per cell type
dars_near_promoters_by_cell_type <- sapply(overlap_hits, function(hit) {
  length(unique(queryHits(hit)))
})

dars_by_cell_type <- sapply(cell_types, function(ct) {
  dars_subset <- dars_gr[dars_gr$cell_type == ct]
  length(unique(dars_subset))
})






# old code below
# Calculate overlaps between DARs and DEG promoters
overlaps <- findOverlaps(dars_gr, degs_promoters)
dars_near_promoters <- length(unique(queryHits(overlaps)))
print(dars_near_promoters)  # Number of unique DARs overlapping DEG promoters

# Calculate distance to the nearest DEG promoter for each DAR
dist_to_nearest <- distanceToNearest(dars_gr, degs_promoters)

# Create a data frame summarizing DAR-DEG promoter relationships
results <- data.frame(
  query = queryHits(dist_to_nearest),    # DAR index
  subject = subjectHits(dist_to_nearest), # DEG promoter index
  distance = mcols(dist_to_nearest)$distance # Distance
)

# Join DARs with the nearest DEG promoters
nearest_degs <- dars_gr %>%
  join_nearest(degs_promoters) %>%
  mutate(distance = mcols(dist_to_nearest)$distance)

# Combine DARs, distances, and DEG promoter metadata into a single data frame
dars_with_distances_and_metadata <- cbind(
  as.data.frame(dars_gr), 
  distance = mcols(dist_to_nearest)$distance, 
  as.data.frame(degs_promoters[subjectHits(dist_to_nearest)])  # Extract promoter metadata
)

colnames(dars_with_distances_and_metadata) <- c("seqnames_dar", "start_dar", "end_dar", "width_dar",
                                                "strand_dar", "cell_type_dar", "logFC_dar","distance",
                                                "seqnames_degs","start_deg", "end_deg", "width_deg","strand_deg",
                                                "cell_type_deg", "gene_name", "gene_biotype", "DAR_deg_overlap", "log2FC_deg")

write_csv(dars_with_distances_and_metadata,"working_data/dars_with_distances_and_metadata.csv")

metadata <- mcols(dist_to_nearest)
dars_gr_with_metadata <- dars_gr
mcols(dars_gr_with_metadata) <- cbind(mcols(dars_gr), metadata)

print(dars_gr_with_metadata)


