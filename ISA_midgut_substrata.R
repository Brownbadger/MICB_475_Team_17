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
