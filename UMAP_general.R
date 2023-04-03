#Loading libraries
library(tidyverse)
library(umap)
library(RColorBrewer)

#Loading in data
metafp <- "SpinyT3_fish_metadata.txt"
meta <- read_delim(metafp, delim="\t")

otufp <- "feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

#### FILTERING & DATA MANIPULATION #### 

#Processing the otu data frame and filtering it for UMAP
otu_fish <- otu %>% 
  select(where(is.numeric)) %>% 
  t() %>% 
  as.data.frame()

otu_fish$rowsum <- rowSums(otu_fish)

otu_fish <- otu_fish %>% 
  filter(rowsum >= 4395) %>% 
  select(-rowsum) %>% 
  t() %>% 
  as.data.frame()

otu_fish$colsum <- rowSums(otu_fish)

otu_fishy <- otu_fish %>% 
  filter(colsum != 0) %>% 
  select(-colsum) %>% 
  t() %>% 
  as.data.frame()

survivors <- rownames(otu_fishy)

rownames(otu_fishy) <- NULL

#Pruning the metadata file to match our truncated samples

meta_fish <- meta %>% 
  select(where(is.character)) %>%
  filter(X.SampleID %in% survivors) %>% 
  mutate(ID=row_number())
  

#### UMAP GENERATION ####

set.seed(0)

umap_fit_fish <- otu_fishy %>%
  mutate(ID=row_number()) %>% 
  column_to_rownames("ID") %>%
  scale() %>% 
  umap()

# Generating UMAP dataframe by joining it with the metadata
umap_df_fish <- umap_fit_fish$layout %>%
  as.data.frame() %>%
  rename(UMAP1="V1",
         UMAP2="V2") %>%
  mutate(ID=row_number()) %>%
  inner_join(meta_fish, by="ID")

#Plotting UMAP Object
umap_df_fish %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2,
             color = sample_type,
             shape = substrata_collection)) +
  geom_point(size = 3) +theme(legend.position = c(0.87, 0.87)) +
  theme_classic() +
  scale_color_discrete(name = "Body Site",
                       breaks = unique(umap_df_fish$sample_type),
                       labels = c("Gill", "Midgut", "Skin", "Hindgut")) +
  scale_shape_discrete(name = "Substrata",
                       breaks = unique(umap_df_fish$substrata_collection),
                       labels = c("Rocky Reef", "Sandy Bottom", "Kelp Forest", "Sandy Mud Bottom")) +
  labs(x = "UMAP1",
       y = "UMAP2",
       title = "UMAP Projection of 2216 ASVs") +
  theme(plot.title = element_text(hjust = 0.5)) #+
  #scale_color_brewer(palette = "Set2")

