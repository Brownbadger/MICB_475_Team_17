library(tidyverse)

meta <- read.delim("SpinyT3_fish_metadata_NoHashtag.tsv")

skin_metadata <- meta %>%
  filter(sample_type == "skin")
skin_metadata

midgut_metadata <- meta %>%
  filter(sample_type == "midgut")

hindgut_metadata <- meta %>%
  filter(sample_type == "hindgut")

gill_metadata <- meta %>%
  filter(sample_type == "gill")

write.table(skin_metadata, file = './skin_metadata.tsv', col.names = TRUE,
            row.names = FALSE, sep = "\t")

write.table(midgut_metadata, file = './midgut_metadata.tsv', col.names = TRUE,
            row.names = FALSE, sep = "\t")

write.table(hindgut_metadata, file = './hindgut_metadata.tsv', col.names = TRUE,
            row.names = FALSE, sep = "\t")

write.table(gill_metadata, file = './gill_metadata.tsv', col.names = TRUE,
            row.names = FALSE, sep = "\t")