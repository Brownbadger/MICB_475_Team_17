# Rarefaction - Brandon Connor - 13th March 2023 - Andreas Felber 24th April 2023


# Load libraries
library(phyloseq)
library(ape)
library(tidyverse)

# Load files

metafp <- "SpinyT3_fish_metadata.txt"
meta <- read_delim(metafp, delim="\t")

otufp <- "feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

taxfp <- "taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "tree.nwk"
phylotree <- read.tree(phylotreefp)

#### Format OTU table ####
# save everything except first column (OTU ID) into a matrix
otu_mat <- as.matrix(otu[,-1])
# Make first column (#OTU ID) the rownames of the new matrix
rownames(otu_mat) <- otu$`#OTU ID`
# Use the "otu_table" function to make an OTU table
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 
class(OTU)

#### Format sample metadata ####
# Save everything except sampleid as new data frame
samp_df <- as.data.frame(meta[,-1])
# Make sampleids the rownames
rownames(samp_df)<- meta$`X.SampleID`
# Make phyloseq sample data with sample_data() function
SAMP <- sample_data(samp_df)
class(SAMP)

#### Formatting taxonomy ####
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

#### Create phyloseq object ####
# Merge all into a phyloseq object
phyloseq_fish <- phyloseq(OTU, SAMP, TAX, phylotree)

phy_tree(phyloseq_fish)

# Rarefaction curve
library(vegan)

rarecurve(t(as.data.frame(otu_table(phyloseq_fish))), cex = 0.1, label = FALSE, xlab = "Sample Depth", ylab = "Species Count")
abline(v = 4395)

fish_rare <- rarefy_even_depth(phyloseq_fish, rngseed = 1, sample.size = 4395)

save(phyloseq_fish, file="phyloseq_fish.RData")
save(fish_rare, file="rare_fish.RData")

#write.table(meta, file='SpinyT3_final_metadata.tsv', quote=FALSE, sep='\t', col.names = NA)
