#### Preliminary Steps ####

# Load libraries
library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)
library(indicspecies)

# Load phyloseq object
load("rare_fish.RData")

#### Indicator Species Analysis ####

# Set seed for ISA
set.seed(0)

# Run ISA for substrata
isa_fish <- multipatt(t(otu_table(fish_rare)), cluster = sample_data(fish_rare)$substrata_collection)

# Look at results
summary(isa_fish)

# Extract taxonomy table
taxtable <- tax_table(fish_rare) %>% as.data.frame() %>% rownames_to_column(var="ASV")

# Merge taxonomy table with phyloseq object and filter by significant p-value
isa_res <- isa_fish$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>% 
  filter(p.value <= 0.05) 

# View results
View(isa_res)

# Save results
save(isa_res, file = "isa_res_substrata.RData")
