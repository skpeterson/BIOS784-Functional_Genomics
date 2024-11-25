#visualizing the overlaps 

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
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# log distance by the cell_type (identified by the dar)
# notes that a difference in the log at this scale is pretty significant
dars_with_distances_and_metadata %>%
  ggplot(aes(x = cell_type_dar, y = log(distance), color = cell_type_dar))+
  geom_boxplot()

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
  table()
  
dars_with_distances_and_metadata %>%
  mutate(less_than100 = distance < 100) %>%
  select(less_than100, DAR_deg_overlap) %>%
  table()

dars_with_distances_and_metadata_by_cell_type %>%
  filter(distance == 0) %>%
  select(distance, DAR_deg_overlap) %>%
  table()

#check correlation between between log2FC and log2FC

cor(dars_with_distances_and_metadata$logFC_dar,dars_with_distances_and_metadata$log2FC_deg)

dars_with_distances_and_metadata %>%
  ggplot(aes(logFC_dar, log2FC_deg, color = DAR_deg_overlap)) +
  geom_point() + geom_smooth(color = "black") + facet_wrap(~cell_type_dar) #+ 
  #labs(x = "Log2FC from the DAR", y = "Log2FC from the DEG")

# make a plot for if they agree in sign
  




