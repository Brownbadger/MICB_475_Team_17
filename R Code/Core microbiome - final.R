# Core microbiome on substrata (eliminated sandy bottom)

#!/usr/bin/env Rscript
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)

#### Load data ####
load("fish_midgut.RData")

# Filter dataset by substrata
substrata_sandy_bottom <- subset_samples(fishmidgut, substrata_collection=="sandy bottom")
substrata_rocky_reef <- subset_samples(fishmidgut, substrata_collection=="rocky reef")
substrata_kelp_forest <- subset_samples(fishmidgut, substrata_collection=="kelp forest")

#### Venn diagram of all the ASVs in three substrata: rocky reef, sandy bottom and kelp forest  ####
sandy_bottom_list <- core_members(substrata_sandy_bottom, detection=0.001, prevalence = 0.5)
rocky_reef_list <- core_members(substrata_rocky_reef, detection=0.001, prevalence = 0.5)
kelp_forest_list <- core_members(substrata_kelp_forest, detection=0.001, prevalence = 0.5)

list(SandyBottom = sandy_bottom_list, RockyReef = rocky_reef_list, KelpForest=kelp_forest_list)

taxonomy <- as.data.frame(tax_table(fishmidgut))

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
              , filename = "venndiagram_substrata.png", output=TRUE) +
  labs(fill="Count") 

### As rocky reef has no ASVs that are present over 50% of samples, so eliminate rocky reef

# sandy_mud_bottom_list <- core_members(substrata_sandy_mud_bottom, detection=0.001, prevalence = 0.5)
# kelp_forrest_list <- core_members(substrata_kelp_forrest, detection=0.001, prevalence = 0.5)
# 
# list(SandyMudBottom = sandy_mud_bottom_list, KelpForrest=kelp_forrest_list)
# 
# ggVennDiagram(x=list(SandyMudBottom = sandy_mud_bottom_list, KelpForrest=kelp_forrest_list)
#               , filename = "/Users/shutong/Desktop/fish475_project/venndiagram_substrata.png", output=TRUE) +
#   labs(fill="Count") 


citation("phyloseq")
