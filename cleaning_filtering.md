---
title: "cleaning_filtering"
output: html_document
date: "2023-03-04"
---

library(tidyverse)

primer_1_txt <- "C:/Users/Claudia/OneDrive - UBC/MICB 475 Bioinformatics/Project/fish/sequencing_info/13414_prep_9748_20220714-121058.txt"
primer_1 <- read.delim(primer_1_txt) %>%
select(sample_name, barcode, primer)
# Primer for set 1 is: GTGTGYCAGCMGCCGCGGTAA

primer_2_txt <- "C:/Users/Claudia/OneDrive - UBC/MICB 475 Bioinformatics/Project/fish/sequencing_info/13414_prep_9770_20220714-121058.txt"
primer_2 <- read.delim(primer_2_txt) %>%
select(sample_name, barcode, primer)
# Primer for set 2 is: GTGTGYCAGCMGCCGCGGTAA

primer_3_txt <- "C:/Users/Claudia/OneDrive - UBC/MICB 475 Bioinformatics/Project/fish/sequencing_info/13414_prep_10454_20220714-121058.txt"
primer_3 <- read.delim(primer_3_txt) %>%
select(sample_name, barcode, primer)
# Primer for set 3 is: GTGTGYCAGCMGCCGCGGTAA

primer_4_txt <- "C:/Users/Claudia/OneDrive - UBC/MICB 475 Bioinformatics/Project/fish/sequencing_info/13414_prep_10466_20220714-115709.txt"
primer_4 <- read.delim(primer_4_txt) %>%
select(sample_name, barcode, primer)
# Primer for set 4 is: GTGTGYCAGCMGCCGCGGTAA

primer_5_txt <- "C:/Users/Claudia/OneDrive - UBC/MICB 475 Bioinformatics/Project/fish/sequencing_info/13414_prep_10522_20210721-151753.txt"
primer_5 <- read.delim(primer_5_txt) %>%
select(sample_name, barcode, primer)
# Primer for set 5 is: GTGTGYCAGCMGCCGCGGTAA

#They're all GTGTGYCAGCMGCCGCGGTAA
