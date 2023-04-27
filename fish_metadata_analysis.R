###Used in exploratory analysis of metatdata variables###

#Load library
library(tidyverse)

#The fish dataset blub blub
fishMD <- "fish_metadata.txt"
fish <- read.delim(file=fishMD)

# Use it to look at column names to pick variables to explore
colnames(fish)

#Independent variables, i.e. columns whose relationships with other variables you want to measure
ind <- c("Depth_m", "conservation_status", "swim_mode", "trophic", "common_name", "substrata_collection", "habitat_depth_level2")

#Dependent variables, MUST BE NUMERIC
dep <- c("mass_g", "fl_cm", "gape_cm", "gi_cm", "tl_cm", "dist_to_dorsal_cm", "distance_from_shore_m", "distance_from_shore_m_log",
         "host_body_mass_index", "ratio_dorsal_to_tl", "ratio_gape_to_tl", "ratio_gi_to_tl")

#The true function. Give it your target dataset, a vector of independent variables, and a vector of numeric dependent variables and you will receive
#the largest range of standard scores across each independent variable with the correlated dependent variable listed. Larger gap = more variability.
NewJeans <- function(fish, ind, dep){
  p <- 1
  danielle <- c()
  angie <- c()
  jungle <- c()
  while(p <= length(ind) + 1){
    if(p > length(ind)){
      for(z in c(1:length(danielle))){
        print(paste("Factor", ind[z], "has", "a sscore gap of", as.character(jungle[z]), "for", dep[danielle[z]]))  
      }
    } else {
      for(q in c(1:length(dep))){
        #angie <- c()
        #jungle <- c()
        sscore_col <- fish %>% 
          mutate(clouds = as.numeric(fish[,dep[q]])) %>% 
          group_by(fish[,ind[p]]) %>% 
          filter(!is.na(clouds)) %>%
          summarize(average_clouds = mean(clouds)) %>% 
          mutate(sscore = (average_clouds - mean(average_clouds))/sd(average_clouds)) %>% 
          select(sscore)
        sscore_gap <- max(sscore_col) - min(sscore_col)
        angie <- append(angie, sscore_gap)
        #jungle <- append(jungle, max(angie))
      }
    }
    danielle <- append(danielle, which.max(angie))
    jungle <- append(jungle, max(angie))
    angie <- c()
    p <- p + 1
  }
}

#Functionally very similar to NewJeans, although you receive a stat for every combination of independent and dependent variables. Great when you
#want to look across a few variables, but it can get very messy very quickly if your list of variables starts to get longer.
NewJeans_Ditto <- function(fish, ind, dep){
  for(p in c(1:length(ind))){
    for(q in c(1:length(dep))){
      sscore_col <- fish %>% 
        mutate(clouds = as.numeric(fish[,dep[q]])) %>% 
        group_by(fish[,ind[p]]) %>% 
        filter(!is.na(clouds)) %>%
        summarize(average_clouds = mean(clouds)) %>% 
        mutate(sscore = (average_clouds - mean(average_clouds))/sd(average_clouds)) %>% 
        select(sscore)
      sscore_gap <- max(sscore_col) - min(sscore_col)
      print(paste("Factor", ind[p], "has", "a SScore gap of", as.character(sscore_gap), "for", dep[q]))
    } 
  }  
}


NewJeans(fish, ind, dep)


NewJeans_Ditto(fish, ind, dep)





# Prototype to start testing the loop for independent variables
# for(q in c(1:length(dep))){
#   sscore_col <- fish %>% 
#     mutate(clouds = as.numeric(fish[,dep[q]])) %>% 
#     group_by(fish[,ind[1]]) %>% 
#     filter(!is.na(clouds)) %>% 
#     summarize(average_clouds = mean(clouds)) %>% 
#     mutate(sscore = (average_clouds - mean(average_clouds))/sd(average_clouds)) %>% 
#     select(sscore)
#   sscore_gap <- max(sscore_col) - min(sscore_col)
#   #danielle <- append(danielle, which.min())
#   print(paste("Factor", ind[1], "has", "a SScore gap of", as.character(sscore_gap), "for", dep[q]))
# } 


# Prototype to test the loop for dependent variables
# for(q in c(1:length(dep))){
#   sscore_col <- fish %>% 
#     mutate(clouds = as.numeric(fish[,dep[q]])) %>% 
#     #filter(!is.na(clouds)) %>% 
#     group_by(conservation_status) %>% 
#     filter(!is.na(clouds)) %>% 
#     summarize(average_clouds = mean(clouds)) %>% 
#     mutate(sscore = (average_clouds - mean(average_clouds))/sd(average_clouds)) %>% 
#     select(sscore)
#   sscore_gap <- max(sscore_col) - min(sscore_col)
#   #danielle <- append(danielle, which.min())
#   print(paste("Factor", "Conservation Status", "has", "a SScore gap of", as.character(sscore_gap), "for", dep[q]))
# } 


# Initial prototype to understand the structure of my dataset manipulations
# scol <- fish %>% 
#   mutate(mass_g = as.numeric(mass_g)) %>% 
#   filter(!is.na(mass_g)) %>% 
#   group_by(conservation_status) %>% 
#   #mutate(sscore = mass_g-mean(mass_g)/sd(mass_g)) %>% 
#   summarize(average_mass_g = mean(mass_g)) %>% 
#   mutate(sscore = (average_mass_g - mean(average_mass_g))/sd(average_mass_g)) %>% 
#   select(sscore)


# Me initially messing around with the table
# fish %>% 
#   mutate(mass_g = as.numeric(mass_g)) %>% 
#   filter(!is.na(mass_g)) %>% 
#   group_by(conservation_status) %>% 
#   #mutate(sscore = mass_g-mean(mass_g)/sd(mass_g)) %>% 
#   summarize(average_mass_g = mean(mass_g)) %>% 
#   mutate(sscore = (average_mass_g - mean(average_mass_g))/sd(average_mass_g)) %>% 
#   select(sscore) %>% 
#   max()
# #filter(!is.na(average_mass_g))
# #mutate(tsd = average_mass_g/sd_mass_g)
