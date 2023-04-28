####AFEGCS523-XXX: Analysis of Data from Specific Sampling Locations####

#Remove duplicate column causing QIIME2 Error
metafp <- "SpinyT3_fish_metadata.txt"
meta <- read_delim(metafp, delim="\t")

meta2 <- meta %>%
  select(-SAMPLE_TYPE, -X)
write.table(meta2, file = './spinyT3_metadata_fixed.txt', col.names = TRUE,
            row.names = FALSE, quote = FALSE, sep = "\t")