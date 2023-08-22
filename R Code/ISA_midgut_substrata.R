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

midgut <- fish_rare %>% 
  subset_samples(sample_type=="midgut") %>% 
  subset_samples(substrata_collection!="sandy mud bottom")
  

# Run ISA for substrata
isa_fish_midgut <- multipatt(t(otu_table(midgut)), cluster = sample_data(midgut)$substrata_collection)

# Look at results
summary(isa_fish_midgut)

# Extract taxonomy table
taxtable <- tax_table(midgut) %>% as.data.frame() %>% rownames_to_column(var="ASV")

# Merge taxonomy table with phyloseq object and filter by significant p-value
isa_res_midgut <- isa_fish_midgut$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>% 
  filter(p.value <= 0.05) 

# View results
View(isa_res_midgut)

# Save results
save(isa_res_midgut, file = "isa_res_midgut_substrata.RData")

#### Skin Analysis ####

skin <- fish_rare %>% 
  subset_samples(sample_type=="skin") %>% 
  subset_samples(substrata_collection!="sandy mud bottom")


# Run ISA for substrata
isa_fish_skin <- multipatt(t(otu_table(skin)), cluster = sample_data(skin)$substrata_collection)

# Look at results
summary(isa_fish_skin)

# Extract taxonomy table
taxtable <- tax_table(skin) %>% as.data.frame() %>% rownames_to_column(var="ASV")

# Merge taxonomy table with phyloseq object and filter by significant p-value
isa_res_skin <- isa_fish_skin$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>% 
  filter(p.value <= 0.05) 

# View results
View(isa_res_skin)

# Save results
save(isa_res_skin, file = "isa_res_skin_substrata.RData")

#### Gill Analysis ####

gill <- fish_rare %>% 
  subset_samples(sample_type=="gill") %>% 
  subset_samples(substrata_collection!="sandy mud bottom")


# Run ISA for substrata
isa_fish_gill <- multipatt(t(otu_table(gill)), cluster = sample_data(gill)$substrata_collection)

# Look at results
summary(isa_fish_gill)

# Extract taxonomy table
taxtable <- tax_table(gill) %>% as.data.frame() %>% rownames_to_column(var="ASV")

# Merge taxonomy table with phyloseq object and filter by significant p-value
isa_res_gill <- isa_fish_gill$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>% 
  filter(p.value <= 0.05) 

# View results
View(isa_res_gill)

# Save results
save(isa_res_gill, file = "isa_res_gill_substrata.RData")

#### Hindgut Analysis

hindgut <- fish_rare %>% 
  subset_samples(sample_type=="hindgut") %>% 
  subset_samples(substrata_collection!="sandy mud bottom")


# Run ISA for substrata
isa_fish_hindgut <- multipatt(t(otu_table(hindgut)), cluster = sample_data(hindgut)$substrata_collection)

# Look at results
summary(isa_fish_hindgut)

# Extract taxonomy table
taxtable <- tax_table(hindgut) %>% as.data.frame() %>% rownames_to_column(var="ASV")

# Merge taxonomy table with phyloseq object and filter by significant p-value
isa_res_hindgut <- isa_fish_hindgut$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>% 
  filter(p.value <= 0.05) 

# View results
View(isa_res_hindgut)

# Save results
save(isa_res_hindgut, file = "isa_res_hindgut_substrata.RData")
