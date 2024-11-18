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
metadata <- mcols(dist_to_nearest)
dars_gr_with_metadata <- dars_gr
mcols(dars_gr_with_metadata) <- cbind(mcols(dars_gr), metadata)

print(dars_gr_with_metadata)



