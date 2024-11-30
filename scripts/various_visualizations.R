#visualizing the overlaps 



#### look at how many genes overlap each DAR
# tells us if there are few dars with many overlaps, or many dars with few overlaps

overlaps_per <- count_overlaps(dars_gr, degs_promoters) %>% as.data.frame()

overlaps_df <- data.frame(
  DAR = seq_along(overlaps_per$.),
  Overlaps = overlaps_per$.
)

overlaps_per_hist <- ggplot(overlaps_df, aes(x = Overlaps)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
  scale_y_continuous(transform = "log10") +
  labs(
    title = "",
    x = "Number of Overlaps",
    y = "Frequency"
  ) +
  theme_classic()

ggsave('figures/hist_freq_num_overlaps_per_DAR_DEG.png', overlaps_per_hist, units = 'in', width = 5, height = 3)


# look at the MHC regions to see if this is where the high level of overlaps occurs 

mhc_region <- GRanges(seqnames = "chr6", ranges = IRanges(start = 28477797, end = 33448354))
overlaps_in_mhc <- subsetByOverlaps(dars_gr, mhc_region)

genes_in_mhc <- subset(degs_promoters, seqnames == "chr6" & start >= 28477797 & end <= 33448354)
print(as.data.frame(genes_in_mhc))


# check if the the nearest gene is a DEG

dars_with_distances_and_metadata %>%
  mutate(is_nearest_deg = gene_name %in% degs_all$feature) %>%
  group_by(is_nearest_deg) %>%
  summarise(count_true = sum(is_nearest_deg),
            count_false = sum(!is_nearest_deg))
  #NOTE: dar is harder to determine cell type


dars_with_distances_and_metadata %>%
  ggplot(aes(x = seqnames_dar, y = distance, color = cell_type_dar))+
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  theme_classic()

# log distance by the cell_type (identified by the dar)
# notes that a difference in the log at this scale is pretty significant
dist_by_celltype_plot <- dars_with_distances_and_metadata %>%
  ggplot(aes(x = cell_type_dar, y = log(distance), color = cell_type_dar))+
  geom_boxplot() +
  ggtitle('Distance to Nearest DEG by Cell Type') +
  labs(x = "Cell Type", y = "Log Distance to Nearest DEG")+
  scale_color_manual(values = natparks.pals('Yellowstone', 5)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("figures/dist_by_celltype_plot.png", plot = dist_by_celltype_plot)

# log distance by the gene_biotype 

dist_by_biotype_plot <- dars_with_distances_and_metadata %>%
  ggplot(aes(x = gene_biotype, y = log(distance), color = gene_biotype))+
  geom_boxplot() + 
  ggtitle('Distance to Nearest DEG by Gene Biotype') +
  labs(x = "Gene Biotype", y = "Log Distance to Nearest DEG")+
  scale_color_manual(values = natparks.pals('KingsCanyon', 11)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("figures/dist_by_biotype_plot.png", plot = dist_by_biotype_plot)


dars_with_distances_and_metadata %>%
  group_by(cell_type_dar) %>%
  summarise(mean = mean(distance, na.rm = T))

pca <- prcomp(dars_with_distances_and_metadata[, c("distance","logFC_dar")], scale. = TRUE)
pca_df <- data.frame(pca$x)
pca_df$cell_type_dar <- dars_with_distances_and_metadata$cell_type_dar
#pca_df$peak_type <- dars_all$peakType


#so much data
ggplot(pca_df, aes(x = PC1, y = PC2, color = cell_type_dar)) +
  geom_jitter() +
  labs(title = "PCA Plot by Cell Type") 

#randomly subset to 1000 
pca_df %>%
  sample_n(1000) %>%
  ggplot(aes(x = PC1, y = PC2, color = cell_type_dar))+
  geom_jitter() +
  labs(title = "PCA Plot by Cell Type") 
#still no real clustering

dars_with_distances_and_metadata %>%
  ggplot(aes(x = cell_type_dar, y = log(distance), color = cell_type_dar))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dars_with_distances_and_metadata %>%
  group_by(cell_type_dar) %>%
  summarise(mean = mean(distance))


#checking if DAR_deg_overlap matches with our distance

look <- dars_with_distances_and_metadata %>%
  filter(DAR_deg_overlap == T) %>%
  select(seqnames_dar, start_dar, end_dar, distance, DAR_deg_overlap)

dars_with_distances_and_metadata %>%
  filter(DAR_deg_overlap == T) %>%
  count()
#there are 169 with overlap

dars_with_distances_and_metadata %>%
  filter(distance == 0) %>%
  select(distance, DAR_deg_overlap) %>%
  table() %>%
  as.data.frame() %>%
  write_csv("figures/distance_overlap_check.csv")
  
dars_with_distances_and_metadata %>%
  mutate(less_than100 = distance < 100) %>%
  select(less_than100, DAR_deg_overlap) %>%
  table() %>%
  as.data.frame() %>%
  write_csv("figures/distance_less100_overlap_check.csv")

dars_with_distances_and_metadata_by_cell_type %>%
  filter(distance == 0) %>%
  select(distance, DAR_deg_overlap) %>%
  table()

#CORRELATION BETWEEN LOG2FC of DEG and DAR (IN REPORT)
#check correlation between between log2FC and log2FC

cor(dars_with_distances_and_metadata$logFC_dar,dars_with_distances_and_metadata$log2FC_deg)

corr_by_celltype_plot <- dars_with_distances_and_metadata %>%
  ggplot(aes(logFC_dar, log2FC_deg, color = cell_type_dar)) +
  geom_point() + 
  geom_smooth(color = "black") +
  scale_color_manual(values = natparks.pals('DeathValley',5)) +
  facet_wrap(~cell_type_dar)+ 
  theme(legend.position="none") + 
  labs(x = "Log2FC DAR", y = "Log2FC DEG", title = "Correlation between LOG 2FC by Cell Type")

ggsave("figures/corr_by_celltype_plot.png", plot = corr_by_celltype_plot)

corr_by_gene_biotype_plot <- dars_with_distances_and_metadata %>%
  ggplot(aes(logFC_dar, log2FC_deg, color = gene_biotype)) +
  geom_point() + 
  geom_smooth(color = "black") + 
  scale_color_manual(values = natparks.pals('DeathValley',11)) +
  facet_wrap(~gene_biotype) + 
  theme(legend.position="none")+ 
  labs(x = "Log2FC DAR", y = "Log2FC DEG", title = "Correlation between LOG 2FC by Gene Biotype")

ggsave("figures/corr_by_gene_biotype_plot.png", plot = corr_by_gene_biotype_plot)

# BARPLOT OF MATCHING SIGNS (NOT IN REPORT)
dars_with_distances_and_metadata %>%
  mutate(sign_match = ifelse(sign(logFC_dar) == sign(log2FC_deg), 1, 0)) %>%
  ggplot(aes(as.factor(sign_match)))+
  geom_bar() + 
  scale_color_manual(values = natparks.pals('SmokyMtns',2)) +
  theme_classic() +
  theme(axis.text.x = element_text(hjust = 1))

#SIGN MATCHING TABLE (IN REPORT)
# Create contingency table of signs
sign_table <- dars_with_distances_and_metadata %>%
  mutate(signs_dar = sign(logFC_dar),
         signs_deg = sign(log2FC_deg)) %>%
  select(signs_dar, signs_deg) %>%
  table() %>%
  as.data.frame()

write_csv(sign_table,"figures/sign_table.csv")

#SIGN MATCHING PLOT (IN REPORT)
sign_match_plot <- ggplot(as.data.frame(sign_table), aes(x = signs_dar, y = signs_deg, fill = Freq)) + 
  geom_tile()+
  scale_fill_gradient(low = "lightgrey", high = "blue") + 
  labs(x = "Sign of the DAR", y = "Sign of the DEG", title = "Sign Matching between DAR and DEG Log2FC")

ggsave("figures/sign_match_plot.png", plot = sign_match_plot)
