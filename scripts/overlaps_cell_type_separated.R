cell_types <- intersect(unique(degs_gr$cell_type), unique(dars_gr$cell_type))

# Step 1: Initialize lists to store results by cell type
results_by_cell_type <- list()
dars_with_distances_and_metadata_by_cell_type <- list()

# Step 2: Loop over each cell type and compute necessary steps
for (ct in cell_types) {
  
  # Subset DEGs promoters and DARs for the current cell type
  degs_subset <- degs_promoters[degs_promoters$cell_type == ct]
  dars_subset <- dars_gr[dars_gr$cell_type == ct]
  
  # Check if both subsets are non-empty to avoid errors
  if (length(degs_subset) == 0 || length(dars_subset) == 0) {
    warning(paste("Skipping cell type", ct, "due to empty DEGs or DARs."))
    next
  }
  
  # Step 3: Calculate the distance to the nearest DEG promoter for each DAR
  dist_to_nearest <- distanceToNearest(dars_subset, degs_subset)
  
  # Create a vector of distances with NA for unmatched DARs
  distances <- rep(NA, length(dars_subset))
  distances[queryHits(dist_to_nearest)] <- mcols(dist_to_nearest)$distance
  
  # Retrieve indices of matched DARs and DEGs
  matched_dars_indices <- queryHits(dist_to_nearest)
  matched_degs_indices <- subjectHits(dist_to_nearest)
  
  # Initialize empty DataFrame for DEG metadata with NA values
  n_dars <- length(dars_subset)
  deg_metadata_full <- DataFrame(
    seqnames_deg = rep(NA_character_, n_dars),
    start_deg = rep(NA_integer_, n_dars),
    end_deg = rep(NA_integer_, n_dars),
    strand_deg = rep(NA_character_, n_dars),
    cell_type_deg = rep(NA_character_, n_dars),
    gene_name = rep(NA_character_, n_dars),
    gene_biotype = rep(NA_character_, n_dars),
    DAR_deg_overlap = rep(NA, n_dars),
    log2FC_deg = rep(NA_real_, n_dars)
  )
  
  # Only proceed if there are any matches
  if (length(matched_dars_indices) > 0) {
    # Extract DEG metadata for matched indices
    degs_seqnames <- as.character(seqnames(degs_subset))
    degs_start <- start(degs_subset)
    degs_end <- end(degs_subset)
    degs_strand <- as.character(strand(degs_subset))
    degs_cell_type <- as.character(mcols(degs_subset)$cell_type)
    degs_gene_name <- as.character(mcols(degs_subset)$gene_name)
    degs_gene_biotype <- as.character(mcols(degs_subset)$gene_biotype)
    DAR_deg_overlap <- as.character(mcols(degs_subset)$DAR_deg_overlap)
    degs_log2FC <- as.numeric(mcols(degs_subset)$log2FC_deg)
    
    # Ensure that the lengths of these vectors match the length of degs_subset
    n_degs <- length(degs_subset)
    
    # Now extract the matched metadata
    deg_metadata_matched <- DataFrame(
      seqnames_deg = degs_seqnames[matched_degs_indices],
      start_deg = degs_start[matched_degs_indices],
      end_deg = degs_end[matched_degs_indices],
      strand_deg = degs_strand[matched_degs_indices],
      cell_type_deg = degs_cell_type[matched_degs_indices],
      gene_name = degs_gene_name[matched_degs_indices],
      gene_biotype = degs_gene_biotype[matched_degs_indices],
      DAR_deg_overlap = DAR_deg_overlap[matched_degs_indices],
      log2FC_deg = degs_log2FC[matched_degs_indices]
    )
    
    # Assign the matched DEG metadata to the corresponding DAR indices
    deg_metadata_full[matched_dars_indices, ] <- deg_metadata_matched
  }
  
  # Combine DARs, distances, and DEG promoter metadata into a GRanges object
  mcols(dars_subset)$distance <- distances
  mcols(dars_subset) <- cbind(mcols(dars_subset), deg_metadata_full)
  
  # Store the GRanges object by cell type
  dars_with_distances_and_metadata_by_cell_type[[ct]] <- dars_subset
  
  # Optional: Print a summary for debugging
  message(paste("Processed cell type:", ct))
}

# Final: To view the GRanges results by cell type
dars_with_distances_and_metadata_by_cell_type  # GRanges with distances and metadata by cell type
# write_csv(dars_with_distances_and_metadata,"working_data/dars_with_distances_and_metadata_by_cell_type.csv")


