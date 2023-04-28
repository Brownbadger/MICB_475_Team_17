####Alpha and Beta Analysis on Midgut Samples####
# https://mibwurrepo.github.io/R_for_Microbial_Ecology/Microbiome_tutorial_V2.html Referred to this tutorial for help with Faith's PD and unifrac in R.

library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)

## Load data 

metafp <- "SpinyT3_fish_metadata.txt"
meta <- read_delim(metafp, delim="\t")

otufp <- "midgut-feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

taxfp <- "taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "tree.nwk"
phylotree <- read.tree(phylotreefp)

## Format OTU table
# OTU tables should be a matrix with rownames and colnames as OTUs and sampleIDs, respectively
# Note: tibbles do not allow rownames so if you imported with read_delim, change back

# save everything except first column (OTU ID) into a matrix
otu_mat <- as.matrix(otu[,-1])
# Make first column (#OTU ID) the rownames of the new matrix
rownames(otu_mat) <- otu$`#OTU ID`
# Use the "otu_table" function to make an OTU table
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 
class(OTU)

## Filter sample metadata
meta <- filter(meta, sample_type == "midgut")
chosen_subs <- c("kelp forrest", "sandy bottom", "rocky reef")
meta_chosen_subs <- filter(meta, substrata_collection %in% chosen_subs)

## Format sample metadata 
# Save everything except sampleid as new data frame
samp_df <- as.data.frame(meta_chosen_subs[,-1])
# Make sampleids the rownames
rownames(samp_df)<- meta_chosen_subs$"#SampleID"
# Make phyloseq sample data with sample_data() function
SAMP <- sample_data(samp_df)
class(SAMP)

## Formatting taxonomy 
# Convert taxon strings to a table with separate taxa rank columns
tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() # Saving as a matrix
# Save everything except feature IDs 
tax_mat <- tax_mat[,-1]
# Make sampleids the rownames
rownames(tax_mat) <- tax$`Feature ID`
# Make taxa table
TAX <- tax_table(tax_mat)
class(TAX)

## Create phyloseq object
# Merge all into a phyloseq object
# fishmidgut <- phyloseq(OTU, SAMP, TAX, phylotree)

## Looking at phyloseq object
# View components of phyloseq object with the following commands
otu_table(fishmidgut)
sample_data(fishmidgut)
tax_table(fishmidgut)
phy_tree(fishmidgut)

save(fishmidgut, file = "fish_midgut.RData")

## Faith's PD  

load("C:/Users/Claudia/OneDrive - UBC/MICB 475 Bioinformatics/Project/fish/fish_midgut.RData")

library(picante)

OTU_table <- otu_table(fishmidgut)
tree <- phy_tree(fishmidgut)
meta_pd <- sample_data(fishmidgut)

# Tree is rooted.

faith_pd <- pd(t(OTU_table), tree, include.root=TRUE)

meta_pd$Phyogenetic_diversity <- faith_pd$PD 

plot_pd <- ggplot(meta_pd, aes(substrata_collection, Phyogenetic_diversity)) + geom_boxplot() + geom_point(size = 2) + theme(axis.text.x = element_text(size=14, angle = 90)) + theme_bw() + xlab("Substrata") + ylab("Phylogenetic Diversity")
print(plot_pd)

ggsave(filename = "faith's pd midgut.png" 
       , plot_pd
       , height=5, width=4)

## Weighted Unifrac 

weighted_uni_ordi <- ordinate(fishmidgut, "PCoA", "unifrac", weighted=T)

weighted_uni <- plot_ordination(fishmidgut, weighted_uni_ordi, color="substrata_collection") +
  ggtitle("Weighted UniFrac") + geom_point(size = 2) +
  scale_color_brewer("Substrata", palette = "Set2")
print(weighted_uni)

ggsave(filename = "weighted unifrac midgut.png"
       , weighted_uni
       , height=4, width=5)

## Unweighted Unifrac 

unweighted_uni_ordi <- ordinate(fishmidgut, "PCoA", "unifrac", weighted=F)

unweighted_uni <- plot_ordination(fishmidgut, unweighted_uni_ordi, color="substrata_collection") +
  ggtitle("Unweighted UniFrac") + geom_point(size = 2) +
  scale_color_brewer("Substrata", palette = "Set2")
print(unweighted_uni)

ggsave(filename = "unweighted unifrac midgut.png"
       , unweighted_uni
       , height=4, width=5)

## Shannon Diversity


shannon <- plot_richness(fishmidgut, measures = c("Shannon")) 

shannon_plot <- plot_richness(fishmidgut, x = "substrata_collection", measures = c("Shannon")) +
  xlab("Substrata") +
  geom_boxplot()
shannon_plot

ggsave(filename = "shannon midgut.png"
       , weighted_uni
       , height=4, width=5)

## Chao Diversity

chao <- plot_richness(fishmidgut, measures = c("Chao1")) 

chao_plot <- plot_richness(fishmidgut, x = "substrata_collection", measures = c("Chao1")) +
  xlab("Substrata") +
  geom_boxplot()
chao_plot

ggsave(filename = "chao midgut.png"
       , weighted_uni
       , height=4, width=5)

dm_chao <- vegdist(t(otu_table(fishmidgut)), method="chao")
chao_stats <- adonis2(dm_chao ~ substrata_collection, data=samp_dat_wdiv)
chao_stats

# rocky reef & kelp forrest chao

rr_kf <- c("rocky reef", "kelp forrest")
rr_kf_meta <- filter(meta_chosen_subs, substrata_collection %in% rr_kf) 

samp_df2 <- as.data.frame(rr_kf_meta[,-1])
rownames(samp_df2)<- rr_kf_meta$"#SampleID"
SAMP2 <- sample_data(samp_df2)
class(SAMP2)

rr_kf_phylo <- phyloseq(OTU, SAMP2, TAX, phylotree)

rr_kf_dm <- vegdist(t(otu_table(rr_kf_phylo)), method="chao")

chao_stats_rr_kf <- adonis2(rr_kf_dm ~ substrata_collection, data=rr_kf_meta)
chao_stats_rr_kf

# sandy bottom & kelp forrest chao

sb_kf <- c("sandy bottom", "kelp forrest")
sb_kf_meta <- filter(meta_chosen_subs, substrata_collection %in% sb_kf) 

samp_df3 <- as.data.frame(sb_kf_meta[,-1])
rownames(samp_df3)<- sb_kf_meta$"#SampleID"
SAMP3 <- sample_data(samp_df3)
class(SAMP3)

sb_kf_phylo <- phyloseq(OTU, SAMP3, TAX, phylotree)

sb_kf_dm <- vegdist(t(otu_table(sb_kf_phylo)), method="chao")
chao_stats_sb_kf <- adonis2(sb_kf_dm ~ substrata_collection, data=sb_kf_meta)
chao_stats_sb_kf

# rocky reef & sandy bottom chao

rr_sb <- c("rocky reef", "sandy bottom")
rr_sb_meta <- filter(meta_chosen_subs, substrata_collection %in% rr_sb) 

samp_df4 <- as.data.frame(rr_sb_meta[,-1])
rownames(samp_df4)<- rr_sb_meta$"#SampleID"
SAMP4 <- sample_data(samp_df4)
class(SAMP4)

rr_sb_phylo <- phyloseq(OTU, SAMP4, TAX, phylotree)

rr_sb_dm <- vegdist(t(otu_table(rr_sb_phylo)), method="chao")
chao_stats_rr_sb <- adonis2(rr_sb_dm ~ substrata_collection, data=rr_sb_meta)
chao_stats_rr_sb