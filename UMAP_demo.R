#Practice Palmer Penguins UMAP
library(tidyverse)
library(palmerpenguins)
library(umap)

#Formating for UMAP
penguins <- penguins %>% 
  drop_na() %>% 
  select(-year) %>% 
  mutate(ID=row_number())

penguins_meta <- penguins %>% 
  select(ID, species, island, sex)

### UMAP GENERATION ###

set.seed(0)
umap_fit <- penguins %>%
  select(where(is.numeric)) %>%
  column_to_rownames("ID") %>%
  scale() %>% 
  umap()

umap_df <- umap_fit$layout %>%
  as.data.frame()%>%
  rename(UMAP1="V1",
         UMAP2="V2") %>%
  mutate(ID=row_number())%>%
  inner_join(penguins_meta, by="ID")

umap_df %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2, 
             color = species,
             shape = sex))+
  geom_point()+
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle = "UMAP plot")
