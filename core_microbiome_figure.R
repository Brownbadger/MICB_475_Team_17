#Load Libraries
library(tidyverse)
library(RColorBrewer)

#Load Data
sandy_bottom_core <- read_tsv("sandy_bottom_core.tsv")
kelp_forest_core <- read_tsv("kelp_forest_core.tsv")
rocky_reef_core <- read_tsv("rocky_reef_core.tsv")

#### Data Processing ####
sandy_bottom_core_unique <- sandy_bottom_core %>% 
  slice(-9)
sandy_bottom_core_unique$Substrata <- "Sandy Bottom"

kelp_forest_core_unique <- kelp_forest_core %>% 
  slice(-10)
kelp_forest_core_unique$Substrata <- "Kelp Forest"

rocky_reef_core_unique <- rocky_reef_core %>% 
  slice(-14)
rocky_reef_core_unique$Substrata <- "Rocky Reef"

core_unique <- sandy_bottom_core %>% 
  slice(9)
core_unique$Substrata <- "Sandy Bottom & Kelp Forest & Rocky Reef"

core_microbes <- rbind (sandy_bottom_core_unique, kelp_forest_core_unique, rocky_reef_core_unique, core_unique)

#Making the Stacked Bar Plot

core_microbes %>% 
  ggplot(aes(fill=Substrata, y=Phylum)) +
  geom_col() +
  labs(x = "Phyla", y = "Phyla Count") +
  ggtitle("She calls me Steven") +
  scale_color_brewer(palette = "Set2")
