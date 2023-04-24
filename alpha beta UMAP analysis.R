###Alpha, beta, UMAP analyses of sampling locations and substrata - Ellie Gang 03.14.23
library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)
library(phyloseq)


#### Load in RData ####
load("/Users/shutong/Desktop/fish475_project/rare_fish.RData")
load("/Users/shutong/Desktop/fish475_project/phyloseq_fish.RData")

#### Alpha diversity - sampling site ######
plot_richness(fish_rare) 

plot_richness(fish_rare, measures = c("Shannon","Chao1")) 

gg_richness_sampling <- plot_richness(fish_rare, x = "SAMPLE_TYPE", measures = c("Shannon","Chao1")) +
  xlab("Body sampling locations") +
  geom_boxplot()
gg_richness_sampling


ggsave(filename = "sampling_site_alpha.png"
       , gg_richness_sampling
       , height=4, width=6)

estimate_richness(fish_rare)

### Evaluating significance of alpha diversity - sampling site ###
#### Statistical analysis using PERMANOVA
samp_dat_wdiv <- data.frame(sample_data(fish_rare), estimate_richness(fish_rare))
dm_chao <- vegdist(t(otu_table(fish_rare)), method="chao") 
adonis2(dm_chao ~ substrata_collection*SAMPLE_TYPE, data=samp_dat_wdiv)
?vegdist

#### Comparison of >2 means with ANOVA (Parametric) by sample site####
alphadiv <- estimate_richness(fish_rare)
samp_dat <- sample_data(fish_rare)
samp_dat_wdiv <- data.frame(samp_dat, alphadiv)

# Set up our linear model 
lm_shannon_vs_SAMPLE <- lm(Shannon ~ SAMPLE_TYPE, dat=samp_dat_wdiv)
# Calculate AOV
anova_shannon_vs_SAMPLE <- aov(lm_shannon_vs_SAMPLE)
# Summarize
summary(anova_shannon_vs_SAMPLE)

TukeyHSD(anova_shannon_vs_SAMPLE)

#### Beta diversity - sampling site#####
bc_dm <- distance(fish_rare, method="Weighted Unifrac")

pcoa_bc <- ordinate(fish_rare, method="PCoA", distance=bc_dm)

plot_ordination(fish_rare, pcoa_bc, color="SAMPLE_TYPE")

gg_pcoa <- plot_ordination(fish_rare,pcoa_bc, color = "SAMPLE_TYPE") +
  labs( col = "body sampling site")
gg_pcoa

ggsave("sampling_site_beta_plot_pcoa.png"
       , gg_pcoa
       , height=4, width=5)

#### Beta diversity - substrata#####
bc_dm_substrata <- distance(fish_rare, method="Weighted Unifrac")

pcoa_bc_substrata <- ordinate(fish_rare, method="PCoA", distance=bc_dm_substrata)

plot_ordination(fish_rare, pcoa_bc_substrata, shape ="substrata_collection")

gg_pcoa_sub <- plot_ordination(fish_rare,pcoa_bc, shape="substrata_collection") +
  labs(pch="substrata")
gg_pcoa_sub

ggsave("substrata_beta_plot_pcoa.png"
       , gg_pcoa_sub
       , height=4, width=5)

#### Beta diversity - sampling site and substrata COMBINED#####
bc_dm_com <- distance(fish_rare, method="Weighted Unifrac")

pcoa_bc_com <- ordinate(fish_rare, method="PCoA", distance=bc_dm_com)

plot_ordination(fish_rare, pcoa_bc_com, color="SAMPLE_TYPE", shape ="substrata_collection")

gg_pcoa_com <- plot_ordination(fish_rare,pcoa_bc, color = "SAMPLE_TYPE", shape="substrata_collection") +
  labs(pch="substrata", col = "body sampling site")
gg_pcoa_com

ggsave("combined_beta_plot_pcoa.png"
       , gg_pcoa_com
       , height=4, width=5)

#### UMAP - sampling site##### -- ISSUES with codes 
# load your omic data here as mydata
library(M3C)
r <- M3C(fish_rare,method=2)
umap(fish_rare,labels=as.factor(r$SAMPLE_TYPE),printres = TRUE,printwidth = 24)

umap(pollen$data,colvec=c('skyblue'))

#### Alpha diversity - substrata ######
plot_richness(fish_rare) 

plot_richness(fish_rare, measures = c("Shannon","Chao1")) 

gg_richness_substrata <- plot_richness(fish_rare, x = "substrata_collection", measures = c("Shannon","Chao1")) +
  xlab("Substrata") +
  geom_boxplot()
gg_richness_substrata

ggsave(filename = "substrata_alpha.png"
       , gg_richness_substrata
       , height=4, width=6)

estimate_richness(fish_rare)

### Evaluating significance of alpha diversity - substrata ###
#### Comparison of >2 means with ANOVA (Parametric) ####
alphadiv <- estimate_richness(fish_rare)
samp_dat <- sample_data(fish_rare)
samp_dat_wdiv <- data.frame(samp_dat, alphadiv)

# Set up our linear model 
lm_shannon_vs_substrata <- lm(log(Shannon) ~ substrata_collection, dat=samp_dat_wdiv)
# Calculate AOV
anova_shannon_vs_substrata <- aov(lm_shannon_vs_substrata)
# Summarize
summary(anova_shannon_vs_substrata)

TukeyHSD(anova_shannon_vs_SAMPLE)

