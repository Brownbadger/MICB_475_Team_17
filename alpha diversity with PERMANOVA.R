###Alpha diversity PERMANOVA analyses of sampling locations and substrata - Ellie Gang
library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)
library(phyloseq)

#### Load in RData ####
load("rare_fish.RData")
load("phyloseq_fish.RData")

#### Alpha diversity - sampling site ######
plot_richness(fish_rare) 

plot_richness(fish_rare, measures = c("Shannon","Chao1")) 

gg_richness_sampling <- plot_richness(fish_rare, x = "SAMPLE_TYPE", measures = c("Shannon","Chao1")) +
  xlab("Body sampling locations") +
  geom_boxplot() 
gg_richness_sampling +theme_classic()

ggsave(filename = "gg_richness_sampling.png"
       , gg_richness_sampling
       , height=4, width=6)

estimate_richness(fish_rare)

### Evaluating significance of alpha diversity - sampling site ###
#### Statistical analysis using PERMANOVA
# gill vs midgut 
samp_dat_wdiv <- data.frame(sample_data(fish_rare), estimate_richness(fish_rare))
gill_midgut_sub <- subset_samples(fish_rare, sample_type %in% c("gill", "midgut"))
dm_chao <- vegdist(t(otu_table(gill_midgut_sub)), method="chao") 
adonis2(dm_chao ~ sample_type, data=filter(samp_dat_wdiv, sample_type %in% c("gill", "midgut")))

# mindgut vs midgut 
samp_dat_wdiv <- data.frame(sample_data(fish_rare), estimate_richness(fish_rare))
hindgut_midgut_sub <- subset_samples(fish_rare, sample_type %in% c("hindgut", "midgut"))
dm_chao <- vegdist(t(otu_table(hindgut_midgut_sub)), method="chao") 
adonis2(dm_chao ~ sample_type, data=filter(samp_dat_wdiv, sample_type %in% c("hindgut", "midgut")))

# skin vs midgut 
samp_dat_wdiv <- data.frame(sample_data(fish_rare), estimate_richness(fish_rare))
skin_midgut_sub <- subset_samples(fish_rare, sample_type %in% c("skin", "midgut"))
dm_chao <- vegdist(t(otu_table(skin_midgut_sub)), method="chao") 
adonis2(dm_chao ~ sample_type, data=filter(samp_dat_wdiv, sample_type %in% c("skin", "midgut")))


# gill vs hindgut 
samp_dat_wdiv <- data.frame(sample_data(fish_rare), estimate_richness(fish_rare))
gill_hindgut_sub <- subset_samples(fish_rare, sample_type %in% c("gill", "hindgut"))
dm_chao <- vegdist(t(otu_table(gill_hindgut_sub)), method="chao") 
adonis2(dm_chao ~ sample_type, data=filter(samp_dat_wdiv, sample_type %in% c("gill", "hindgut")))


# gill vs skin 
samp_dat_wdiv <- data.frame(sample_data(fish_rare), estimate_richness(fish_rare))
gill_skin_sub <- subset_samples(fish_rare, sample_type %in% c("gill", "skin"))
dm_chao <- vegdist(t(otu_table(gill_skin_sub)), method="chao") 
adonis2(dm_chao ~ sample_type, data=filter(samp_dat_wdiv, sample_type %in% c("gill", "skin")))


# hindgut vs skin 
samp_dat_wdiv <- data.frame(sample_data(fish_rare), estimate_richness(fish_rare))
hindgut_skin_sub <- subset_samples(fish_rare, sample_type %in% c("hindgut", "skin"))
dm_chao <- vegdist(t(otu_table(hindgut_skin_sub)), method="chao") 
adonis2(dm_chao ~ sample_type, data=filter(samp_dat_wdiv, sample_type %in% c("hindgut", "skin")))


dm_unifrac <- UniFrac(fish_rare, weighted=TRUE) # Weighted UniFrac
adonis2(dm_unifrac ~ SAMPLE_TYPE, data=samp_dat_wdiv)

# gill vs midgut - chao










### Alpha diversity - substrata ######
plot_richness(fish_rare) 

plot_richness(fish_rare, measures = c("Shannon","Chao1")) 

gg_richness_substrata <- plot_richness(fish_rare, x = "substrata_collection", measures = c("Shannon","Chao1")) +
  xlab("Substrata") +
  geom_boxplot() +theme_classic()
gg_richness_substrata

ggsave(filename = "substrata_alpha.png"
       , gg_richness_substrata
       , height=4, width=6)

estimate_richness(fish_rare)

### Evaluating significance of alpha diversity - substrata ###
#### Statistical analysis using PERMANOVA
samp_dat_wdiv <- data.frame(sample_data(fish_rare), estimate_richness(fish_rare))
dm_chao <- vegdist(t(otu_table(fish_rare)), method="chao") 
adonis2(dm_chao ~ substrata_collection, data=samp_dat_wdiv)

# TEMPLATE hindgut vs skin 
samp_dat_wdiv <- data.frame(sample_data(fish_rare), estimate_richness(fish_rare))
hindgut_skin_sub <- subset_samples(fish_rare, sample_type %in% c("hindgut", "skin"))
dm_chao <- vegdist(t(otu_table(hindgut_skin_sub)), method="chao") 
adonis2(dm_chao ~ sample_type, data=filter(samp_dat_wdiv, sample_type %in% c("hindgut", "skin")))

# Kelp forest vs rocky reef 
samp_dat_wdiv <- data.frame(sample_data(fish_rare), estimate_richness(fish_rare))
kelp_rocky_sub <- subset_samples(fish_rare, substrata_collection %in% c("kelp forrest", "rocky reef"))
dm_chao <- vegdist(t(otu_table(kelp_rocky_sub)), method="chao") 
adonis2(dm_chao ~ substrata_collection, data=filter(samp_dat_wdiv, substrata_collection %in% c("kelp forrest", "rocky reef")))

# Kelp forest vs sandy bottom
samp_dat_wdiv <- data.frame(sample_data(fish_rare), estimate_richness(fish_rare))
kelp_sandB_sub <- subset_samples(fish_rare, substrata_collection %in% c("kelp forrest", "sandy bottom"))
dm_chao <- vegdist(t(otu_table(kelp_sandB_sub)), method="chao") 
adonis2(dm_chao ~ substrata_collection, data=filter(samp_dat_wdiv, substrata_collection %in% c("kelp forrest", "sandy bottom")))


# Kelp forest vs sandy mud bottom
samp_dat_wdiv <- data.frame(sample_data(fish_rare), estimate_richness(fish_rare))
kelp_sandMB_sub <- subset_samples(fish_rare, substrata_collection %in% c("kelp forrest", "sandy mud bottom"))
dm_chao <- vegdist(t(otu_table(kelp_sandMB_sub)), method="chao") 
adonis2(dm_chao ~ substrata_collection, data=filter(samp_dat_wdiv, substrata_collection %in% c("kelp forrest", "sandy mud bottom")))


# rocky reef vs sandy bottom
samp_dat_wdiv <- data.frame(sample_data(fish_rare), estimate_richness(fish_rare))
rocky_sandB_sub <- subset_samples(fish_rare, substrata_collection %in% c("rocky reef", "sandy bottom"))
dm_chao <- vegdist(t(otu_table(rocky_sandB_sub)), method="chao") 
adonis2(dm_chao ~ substrata_collection, data=filter(samp_dat_wdiv, substrata_collection %in% c("rocky reef", "sandy bottom")))

# rocky reef vs sandy mud bottom
samp_dat_wdiv <- data.frame(sample_data(fish_rare), estimate_richness(fish_rare))
rocky_sandMB_sub <- subset_samples(fish_rare, substrata_collection %in% c("rocky reef", "sandy mud bottom"))
dm_chao <- vegdist(t(otu_table(rocky_sandMB_sub)), method="chao") 
adonis2(dm_chao ~ substrata_collection, data=filter(samp_dat_wdiv, substrata_collection %in% c("rocky reef", "sandy mud bottom")))

# sandy bottom vs sandy mud bottom
samp_dat_wdiv <- data.frame(sample_data(fish_rare), estimate_richness(fish_rare))
sandB_sandMB_sub <- subset_samples(fish_rare, substrata_collection %in% c("sandy bottom", "sandy mud bottom"))
dm_chao <- vegdist(t(otu_table(sandB_sandMB_sub)), method="chao") 
adonis2(dm_chao ~ substrata_collection, data=filter(samp_dat_wdiv, substrata_collection %in% c("sandy bottom", "sandy mud bottom")))


dm_unifrac <- UniFrac(fish_rare, weighted=TRUE) # Weighted UniFrac
adonis2(dm_unifrac ~ substrata_collection, data=samp_dat_wdiv)


