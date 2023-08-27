###Alpha diversity PERMANOVA analyses of sampling locations and substrata - Ellie Gang
library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)

#### Load in RData ####
load("rare_fish.RData")

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
samp_dat_wdiv <- data.frame(sample_data(fish_rare), estimate_richness(fish_rare))

# gill vs midgut 
gill_midgut_sub <- subset_samples(fish_rare, sample_type %in% c("gill", "midgut"))
dm_chao <- vegdist(t(otu_table(gill_midgut_sub)), method="chao")
set.seed(0)
adonis2(dm_chao ~ sample_type, data=filter(samp_dat_wdiv, sample_type %in% c("gill", "midgut")))

# hindgut vs midgut
hindgut_midgut_sub <- subset_samples(fish_rare, sample_type %in% c("hindgut", "midgut"))
dm_chao <- vegdist(t(otu_table(hindgut_midgut_sub)), method="chao") 
set.seed(0)
adonis2(dm_chao ~ sample_type, data=filter(samp_dat_wdiv, sample_type %in% c("hindgut", "midgut")))

# skin vs midgut
skin_midgut_sub <- subset_samples(fish_rare, sample_type %in% c("skin", "midgut"))
dm_chao <- vegdist(t(otu_table(skin_midgut_sub)), method="chao")
set.seed(0)
adonis2(dm_chao ~ sample_type, data=filter(samp_dat_wdiv, sample_type %in% c("skin", "midgut")))


# gill vs hindgut
gill_hindgut_sub <- subset_samples(fish_rare, sample_type %in% c("gill", "hindgut"))
dm_chao <- vegdist(t(otu_table(gill_hindgut_sub)), method="chao") 
set.seed(0)
adonis2(dm_chao ~ sample_type, data=filter(samp_dat_wdiv, sample_type %in% c("gill", "hindgut")))


# gill vs skin 
gill_skin_sub <- subset_samples(fish_rare, sample_type %in% c("gill", "skin"))
dm_chao <- vegdist(t(otu_table(gill_skin_sub)), method="chao") 
set.seed(0)
adonis2(dm_chao ~ sample_type, data=filter(samp_dat_wdiv, sample_type %in% c("gill", "skin")))


# hindgut vs skin 
hindgut_skin_sub <- subset_samples(fish_rare, sample_type %in% c("hindgut", "skin"))
dm_chao <- vegdist(t(otu_table(hindgut_skin_sub)), method="chao") 
set.seed(0)
adonis2(dm_chao ~ sample_type, data=filter(samp_dat_wdiv, sample_type %in% c("hindgut", "skin")))


dm_unifrac <- UniFrac(fish_rare, weighted=TRUE) # Weighted UniFrac
adonis2(dm_unifrac ~ SAMPLE_TYPE, data=samp_dat_wdiv)

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
dm_chao <- vegdist(t(otu_table(fish_rare)), method="chao") 
adonis2(dm_chao ~ substrata_collection, data=samp_dat_wdiv)

# Kelp forest vs rocky reef 
kelp_rocky_sub <- subset_samples(fish_rare, substrata_collection %in% c("kelp forrest", "rocky reef"))
dm_chao <- vegdist(t(otu_table(kelp_rocky_sub)), method="chao")
set.seed(0)
adonis2(dm_chao ~ substrata_collection, data=filter(samp_dat_wdiv, substrata_collection %in% c("kelp forest", "rocky reef")))

# Kelp forest vs sandy bottom
kelp_sandB_sub <- subset_samples(fish_rare, substrata_collection %in% c("kelp forrest", "sandy bottom"))
dm_chao <- vegdist(t(otu_table(kelp_sandB_sub)), method="chao") 
set.seed(0)
adonis2(dm_chao ~ substrata_collection, data=filter(samp_dat_wdiv, substrata_collection %in% c("kelp forest", "sandy bottom")))


# Kelp forest vs sandy mud bottom
kelp_sandMB_sub <- subset_samples(fish_rare, substrata_collection %in% c("kelp forrest", "sandy mud bottom"))
dm_chao <- vegdist(t(otu_table(kelp_sandMB_sub)), method="chao") 
set.seed(0)
adonis2(dm_chao ~ substrata_collection, data=filter(samp_dat_wdiv, substrata_collection %in% c("kelp forest", "sandy mud bottom")))


# Rocky reef vs sandy bottom
rocky_sandB_sub <- subset_samples(fish_rare, substrata_collection %in% c("rocky reef", "sandy bottom"))
dm_chao <- vegdist(t(otu_table(rocky_sandB_sub)), method="chao") 
set.seed(0)
adonis2(dm_chao ~ substrata_collection, data=filter(samp_dat_wdiv, substrata_collection %in% c("rocky reef", "sandy bottom")))


# Rocky reef vs sandy mud bottom
rocky_sandMB_sub <- subset_samples(fish_rare, substrata_collection %in% c("rocky reef", "sandy mud bottom"))
dm_chao <- vegdist(t(otu_table(rocky_sandMB_sub)), method="chao") 
set.seed(0)
adonis2(dm_chao ~ substrata_collection, data=filter(samp_dat_wdiv, substrata_collection %in% c("rocky reef", "sandy mud bottom")))


# Sandy bottom vs sandy mud bottom
sandB_sandMB_sub <- subset_samples(fish_rare, substrata_collection %in% c("sandy bottom", "sandy mud bottom"))
dm_chao <- vegdist(t(otu_table(sandB_sandMB_sub)), method="chao") 
set.seed(0)
adonis2(dm_chao ~ substrata_collection, data=filter(samp_dat_wdiv, substrata_collection %in% c("sandy bottom", "sandy mud bottom")))


dm_unifrac <- UniFrac(fish_rare, weighted=TRUE) # Weighted UniFrac
adonis2(dm_unifrac ~ substrata_collection, data=samp_dat_wdiv)


