####Alpha and Beta Analysis on Midgut Samples####
# https://mibwurrepo.github.io/R_for_Microbial_Ecology/Microbiome_tutorial_V2.html Referred to this tutorial for help with Faith's PD and unifrac in R.

library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)

## Amend phyloseq object for sample type and substrata
# Merge all into a phyloseq object
fishmidgut <- fish_rare %>% 
  subset_samples(sample_type == "midgut") %>% 
  subset_samples(substrata_collection != "sandy mud bottom")

## Looking at phyloseq object
# View components of phyloseq object with the following commands
otu_table(fishmidgut)
sample_data(fishmidgut)
tax_table(fishmidgut)
phy_tree(fishmidgut)

save(fishmidgut, file = "fish_midgut.RData")

## Faith's PD

load("fish_midgut.RData")

library(picante)

OTU_table <- otu_table(fishmidgut)
tree <- phy_tree(fishmidgut)
meta_pd <- sample_data(fishmidgut)

# Tree is rooted.

faith_pd <- pd(t(OTU_table), tree, include.root=TRUE)

meta_pd$Phyogenetic_diversity <- faith_pd$PD 

plot_pd <- ggplot(meta_pd, aes(substrata_collection, Phyogenetic_diversity)) + geom_boxplot() + geom_point(size = 2) + theme_classic() + xlab("Substrata") + ylab("Alpha Diversity Measure")
print(plot_pd)

ggsave(filename = "TRUE faith's pd midgut.png"
       , plot_pd
       , height=5, width=6)

## Weighted Unifrac

weighted_uni_ordi <- ordinate(fishmidgut, "PCoA", "unifrac", weighted=T)

weighted_uni <- plot_ordination(fishmidgut, weighted_uni_ordi, color="substrata_collection") +
  ggtitle("Weighted UniFrac") + geom_point(size = 2) + theme_classic() +
  scale_color_brewer("Substrata", palette = "Set2")
print(weighted_uni)

ggsave(filename = "TRUE weighted unifrac midgut.png"
       , weighted_uni
       , height=5, width=6)

## Unweighted Unifrac

unweighted_uni_ordi <- ordinate(fishmidgut, "PCoA", "unifrac", weighted=F)

unweighted_uni <- plot_ordination(fishmidgut, unweighted_uni_ordi, color="substrata_collection") +
  ggtitle("Unweighted UniFrac") + geom_point(size = 2) + theme_classic() + 
  scale_color_brewer("Substrata", palette = "Set2")
print(unweighted_uni)

ggsave(filename = "TRUE unweighted unifrac midgut.png"
       , unweighted_uni
       , height=5, width=6)

## Shannon diversity

shannon <- plot_richness(fishmidgut, measures = c("Shannon")) 

shannon_plot <- plot_richness(fishmidgut, x = "substrata_collection", measures = c("Shannon")) +
  xlab("Substrata") + theme_classic() +
  geom_boxplot()
shannon_plot

ggsave(filename = "TRUE shannon midgut.png"
       , shannon_plot
       , height=5, width=6)

## Chao diversity

chao <- plot_richness(fishmidgut, measures = c("Chao1")) 

chao_plot <- plot_richness(fishmidgut, x = "substrata_collection", measures = c("Chao1")) +
  xlab("Substrata") + theme_classic() +
  geom_boxplot()
chao_plot

ggsave(filename = "TRUE chao midgut.png"
       , chao_plot
       , height=5, width=6)

chosen_subs <- c("kelp forest", "sandy bottom", "rocky reef")
meta_chosen_subs <- filter(filter(data.frame(sample_data(fishmidgut)), substrata_collection %in% chosen_subs), sample_type == "midgut")

dm_chao <- vegdist(t(otu_table(fishmidgut)), method="chao")
set.seed(0)
chao_stats <- adonis2(dm_chao ~ substrata_collection, data=meta_chosen_subs)
chao_stats

# rocky reef & kelp forest chao

rr_kf_phylo <- fishmidgut %>% 
  subset_samples(substrata_collection != "sandy bottom")

rr_kf_dm <- vegdist(t(otu_table(rr_kf_phylo)), method="chao")
set.seed(0)
chao_stats_rr_kf <- adonis2(rr_kf_dm ~ substrata_collection, data=filter(meta_chosen_subs, substrata_collection != "sandy bottom"))
chao_stats_rr_kf

# sandy bottom & kelp forest chao

sb_kf_phylo <- fishmidgut %>% 
  subset_samples(substrata_collection != "rocky reef")

sb_kf_dm <- vegdist(t(otu_table(sb_kf_phylo)), method="chao")
set.seed(0)
chao_stats_sb_kf <- adonis2(sb_kf_dm ~ substrata_collection, data=filter(meta_chosen_subs, substrata_collection != "rocky reef"))
chao_stats_sb_kf

# rocky reef & sandy bottom chao

rr_sb_phylo <- fishmidgut %>% 
  subset_samples(substrata_collection != "kelp forest")


rr_sb_dm <- vegdist(t(otu_table(rr_sb_phylo)), method="chao")
set.seed(0)
chao_stats_rr_sb <- adonis2(rr_sb_dm ~ substrata_collection, data=filter(meta_chosen_subs, substrata_collection != "kelp forest"))
chao_stats_rr_sb

# post hoc testing due to significant Pr (>F) value: not implemented, cannot use phyloseq object
#install.packages('devtools')
#library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis") 
#library(pairwiseAdonis)
#pair.mod <- pairwise.adonis2(fishmidgut ~ "substrata_collection", data= OTU_table)
#pair.mod

