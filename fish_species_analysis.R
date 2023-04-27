#### Load Libraries and dataset ####
library(tidyverse)

fishMD <- "fish_metadata.txt"
fish <- read.delim(file=fishMD)

fishMF <- "fish_manifest.txt"
fish_man <- read.delim(file=fishMF)

#Ray-finned fish (Class:	Actinopterygii)
FishSpecies <- c("kelp bass",
                 "sargo",
                 "longjaw mudsucker",
                 "diamond turbot",
                 "Black croaker",
                 "Pacific chub mackerel",
                 "Pacific 'california' barracuda",
                 "barred sand bass",
                 "queenfish",
                 "pacific staghorn sculpin",
                 "treefish",
                 "CA sheephead",
                 "jack smelt",
                 "ocean whitefish",
                 "Pacific jack mackerel",
                 "zebra perch",
                 "opaleye",
                 "giant kelpfish",
                 "senorita",
                 "Bay pipefish unknown",
                 "longfin sanddab",
                 "calico rockfish",
                 "vermilion rockfish",
                 "half banded rockfish",
                 "honeycomb rockfish",
                 "california kingcroaker 'corbina'",
                 "spotted kelpfish",
                 "rockpool blenny",
                 "Spotted bay bass",
                 "blacksmith",
                 "California flounder",
                 "squarespot rockfish",
                 "greenspotted rockfish",
                 "starry rockfish",
                 "Northern anchovy",
                 "California scorpionfish",
                 "brown rockfish",
                 "california tonguefish",
                 "yellowfin croaker",
                 "california salema",
                 "spotfin croaker",
                 "specklefin midshipman",
                 "shortfin corvina",
                 "gopher rockfish",
                 "bonefish",
                 "California lizardfish",
                 "rubberlip perch",
                 "kelp perch",
                 "halfmoon",
                 "bonito",
                 "Pacific bonito",
                 "olive rockfish",
                 "South American pilchard",
                 "Mexican lampfish",
                 "slender sole",
                 "splitnose rockfish",
                 "dover sole",
                 "hundre fathom mora",
                 "rex sole",
                 "California moray",
                 "Blacktip poacher",
                 "white seaperch",
                 "Earned blacksmelt",
                 "pacific sanddab",
                 "plainfin midshipmen",
                 "english sole",
                 "Northern pacific hake",
                 "bigfin eelpout",
                 "bearded eelpout",
                 "California smoothtongue",
                 "Snipe eel",
                 "yellowtail amberjack",
                 "rock wrasse",
                 "black eelpout",
                 "black belly dragonfish",
                 "black belly eelpout",
                 "Dogface witch eel",
                 "Twospine bigscale",
                 "Broadfin lampfish",
                 "Highlight hatchetfish",
                 "black perch",
                 "Pacific hatchet fish",
                 "Dogtooth lampfish",
                 "giant black seabass",
                 "california grenadier",
                 "yellowfin tuna",
                 "arrow goby",
                 "California killifish",
                 "Dolphinfish",
                 "dolphinfish",
                 "skipjack tuna",
                 "topsmelt silverside",
                 "summer flounder",
                 "atlantic needlefish",
                 "oyster toadfish",
                 "american eel ",
                 "rainwater killifish (fresh-brackish)",
                 "blue runner",
                 "northern searobin",
                 "striped bass",
                 "silver perch",
                 "mummichog",
                 "bluefish",
                 "striped blenny",
                 "atlantic menhaden")

#### Filter manifest file to only contain samples of interest ####

spiny_fish <- fish %>% 
  filter(common_name %in% FishSpecies)

SpinyIDs <- spiny_fish$X.SampleID

write.table(spiny_fish, "Spiny_fish_metadata.txt", sep = "\t", quote = FALSE, row.names = FALSE)

fish_man_ray <- fish_man %>% 
  filter(sample.id %in% SpinyIDs)

write.table(fish_man_ray, "fish_ray_manifest.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#For trophic level 3 fish
spiny_fishT3 <- fish %>% 
  filter(common_name %in% FishSpecies & trophic == 3)

SpinyIDsT3 <- spiny_fishT3$X.SampleID

write.table(spiny_fishT3, "SpinyT3_fish_metadata.txt", sep = "\t", quote = FALSE, row.names = FALSE)

fish_man_rayT3 <- fish_man %>% 
  filter(sample.id %in% SpinyIDsT3)

write.table(fish_man_rayT3, "fish_rayT3_manifest.txt", sep = "\t", quote = FALSE, row.names = FALSE)