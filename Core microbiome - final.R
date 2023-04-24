# Core microbiome on substrate (eliminated sandy bottom)

#!/usr/bin/env Rscript
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)

#### Load data ####
load("/Users/shutong/Desktop/fish475_project/rare_fish.RData")

fish_rare <- fish_rare %>%
  subset_samples(sample_type=="midgut")

# Filter dataset by substrata
substrata_sandy_bottom <- subset_samples(fish_rare, substrata_collection=="sandy bottom")
substrata_sandy_mud_bottom <- subset_samples(fish_rare, substrata_collection=="sandy mud bottom")
substrata_rocky_reef <- subset_samples(fish_rare, substrata_collection=="rocky reef")
substrata_kelp_forrest <- subset_samples(fish_rare, substrata_collection=="kelp forrest")

# What ASVs are found in more than 50% of samples in each vegetation category?
sandy_bottom_ASVs <- core_members(substrata_sandy_mud_bottom, detection=0.001, prevalence = 0.5)
rocky_reef_ASVs <- core_members(substrata_rocky_reef, detection=0.001, prevalence = 0.5)
kelp_forrest_ASVs <- core_members(substrata_kelp_forrest, detection=0.001, prevalence = 0.5)

# What are these ASVs?
prune_taxa(sandy_mud_bottom_ASVs,fish_rare) %>%
  tax_table()

prune_taxa(sandy_mud_bottom_ASVs,fish_rare) %>%
  plot_bar(, fill="Genus")+
  facet_wrap(.~substrata_collection, scales="free")

prune_taxa(kelp_forrest_ASVs,fish_rare) %>%
  tax_table()

prune_taxa(kelp_forrest_ASVs,fish_rare) %>%
  plot_bar(, fill="Genus")+
  facet_wrap(.~substrata_collection, scales="free")

#subset by midgut 
fish_rare <- fish_rare %>%
  subset_samples(sample_type=="midgut")

# Filter dataset by substrata
substrata_sandy_bottom <- subset_samples(fish_rare, substrata_collection=="sandy bottom")
substrata_rocky_reef <- subset_samples(fish_rare, substrata_collection=="rocky reef")
substrata_kelp_forrest <- subset_samples(fish_rare, substrata_collection=="kelp forrest")

### Venn diagram of all the ASVs in three substrata: rocky reef, sandy bottom and kelp forest 
sandy_bottom_list <- core_members(substrata_sandy_bottom, detection=0.001, prevalence = 0.5)
rocky_reef_list <- core_members(substrata_rocky_reef, detection=0.001, prevalence = 0.5)
kelp_forest_list <- core_members(substrata_kelp_forrest, detection=0.001, prevalence = 0.5)

list(SandyBottom = sandy_bottom_list, RockyReef = rocky_reef_list, KelpForest=kelp_forest_list)

taxonomy <- as.data.frame(tax_table(fish_rare))

sandy_bottom_core <- taxonomy %>% 
  rownames_to_column("ASVs") %>%
  filter(ASVs %in% sandy_bottom_list)
write.table(sandy_bottom_core, file = "sandy_bottom_core.tsv", sep = "\t", col.names = NA, row.names = TRUE)

rocky_reef_core <- taxonomy %>% 
  rownames_to_column("ASVs") %>%
  filter(ASVs %in% rocky_reef_list)
write.table(rocky_reef_core, file = "rocky_reef_core.tsv", sep = "\t", col.names = NA, row.names = TRUE)

kelp_forest_core <- taxonomy %>% 
  rownames_to_column("ASVs") %>%
  filter(ASVs %in% kelp_forest_list)
write.table(kelp_forest_core, file = "kelp_forest_core.tsv", sep = "\t", col.names = NA, row.names = TRUE)

ggVennDiagram(x=list(SandyBottom = sandy_bottom_list, RockyReef = rocky_reef_list, KelpForest=kelp_forest_list)
              , filename = "/Users/shutong/Desktop/fish475_project/venndiagram_substrata.png", output=TRUE) +
  labs(fill="Count") 

### As rocky reef has no ASVs that are present over 50% of samples, so eliminate rocky reef
sandy_mud_bottom_list <- core_members(substrata_sandy_mud_bottom, detection=0.001, prevalence = 0.5)
kelp_forrest_list <- core_members(substrata_kelp_forrest, detection=0.001, prevalence = 0.5)

list(SandyMudBottom = sandy_mud_bottom_list, KelpForrest=kelp_forrest_list)

ggVennDiagram(x=list(SandyMudBottom = sandy_mud_bottom_list, KelpForrest=kelp_forrest_list)
              , filename = "/Users/shutong/Desktop/fish475_project/venndiagram_substrata.png", output=TRUE) +
  labs(fill="Count") 


citation("phyloseq")
