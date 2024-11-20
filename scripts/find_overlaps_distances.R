library(GenomicRanges)
library(GenomicDistributions)

# Make chromosome metadata consistent
seqlevels(dars_gr) <- gsub("^chr", "", seqlevels(dars_gr))
seqlevels(degs_gr) <- gsub("^chr", "", seqlevels(degs_gr))

# Calculate overlaps
overlaps <- findOverlaps(dars_gr, degs_gr)

# Count peaks that overlap any DEG regions
dars_near_degs <- length(unique(queryHits(overlaps)))
dars_near_degs

# calculate distance to nearest - dars & degs
dist_to_nearest <- distanceToNearest(dars_gr, degs_gr)

# get results
results <- data.frame(
  query = queryHits(dist_to_nearest),   # ATAC-seq peak index
  subject = subjectHits(dist_to_nearest), # Genomic region index
  distance = mcols(dist_to_nearest)$distance # Distance
)

glimpse(results)

# Merge metadata with original GRanges

# Extract the distances (distances between DARs and nearest DEGs)
distances <- mcols(dist_to_nearest)$distance

# Extract the indices of the nearest DEGs for each DAR
nearest_deg_indices <- subjectHits(dist_to_nearest)

# Extract the metadata of the nearest DEGs (metadata from degs_gr)
nearest_degs_metadata <- degs_gr[nearest_deg_indices]@elementMetadata

# Create a data frame to combine DARs with distances and DEGs metadata
dars_with_distances_and_metadata <- cbind(as.data.frame(dars_gr), distance = distances, nearest_degs_metadata)

colnames(dars_with_distances_and_metadata) <- c("seqnames", "start", "end", "width", "strand", "cell_type_dar", "distance", "cell_type_deg",
                                                "gene_name", "gene_biotype", "DAR_deg_overlap")

write_csv(dars_with_distances_and_metadata,"~/BIOS784-Functional_Genomics/working_data/dars_with_distances_and_metadata.csv")

metadata <- mcols(dist_to_nearest)
dars_gr_with_metadata <- dars_gr
mcols(dars_gr_with_metadata) <- cbind(mcols(dars_gr), metadata)

print(dars_gr_with_metadata)



