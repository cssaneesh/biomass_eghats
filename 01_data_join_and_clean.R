# Packages----
library(tidyverse)
library(patchwork)
library(brms)
library(bayesplot)
library(cmdstanr)
library(cowplot)
library(gt)
library(webshot)

# Data----
raw_dat <- read.csv(
  "biomass_data.csv",
  header = T,
  fill = TRUE,
  sep = ",",
  na.strings = c("", " ", "NA", "NA ", "na", "NULL")
)

# Data wrangling----
Site_dat <- raw_dat %>%
  mutate(
    Treatment = as.factor(Treatment),
    Life_form = as.factor(Life_form),
    Functional_groups = as.factor(Functional_groups)
  )

sp.list <-
  Site_dat %>% select(Family,
                      Scientific_name,
                      Functional_type,
                      Functional_groups,
                      Palatability) %>%
  mutate(
    Functional_type= case_when(
      Functional_type == "Annual undershrub"  ~ "Annual",
      Functional_type == "Annual herb"   ~ "Annual",
      Functional_type == "Annual graminoid"  ~ "Annual",
      Functional_type == "Annual/perennial graminoid "  ~ "Annual",
      Functional_type == "Annual/perennial herb"  ~ "Annual",
      Functional_type == "Perennial herb"   ~ "Perennial",
      Functional_type == "Perennial graminoid"   ~ "Perennial",
      Functional_type == "Perennial undershrub"  ~ "Perennial")) %>% 
  mutate(
    Scientific_name = recode(Scientific_name, 'Cymbopogon sp.' = 'Cymbopogon caesius (Hook. & Arn.) Stapf')
  ) %>%
  unique()

cymbopogons <- data.frame(
  Family = c('Poaceae', 'Poaceae'),
  Scientific_name = c(
    'Cymbopogon jwarancusa (Jones) Schult',
    'Cymbopogon martinii (Roxb.) J. F. Watson'
  ),
  Functional_type = c('Perennial', 'Perennial'),
  Functional_groups = c('Graminoid', 'Graminoid'),
  Palatability = c('No', 'No')
)

unique.sp.list <- bind_rows(sp.list, cymbopogons) %>% 
  mutate(Palatability= recode(Palatability, 'Cymbopogon sp.'= 'No')) %>% 
  mutate(Functional_groups= recode(Functional_groups, 'Cymbopogon'= 'Graminoid'))

# write.csv(unique.sp.list, 'sp.list_supp-1.csv')

# number of graminoids (annual/perennial) and forbs (annual/perennial)
anu.peri <- unique.sp.list %>%
  distinct(Functional_groups,
           Functional_type,
           Scientific_name,
           Palatability) %>%
  mutate(anu.peri = as.factor(Functional_groups)) %>%
  dplyr::count(Functional_groups, Functional_type, Palatability, name = 'Species')
anu.peri

Site_prep <- Site_dat %>%
  arrange(Site, Treatment) %>%
  mutate(
    Treatment = case_when(
      Treatment == "ab" ~ "Control",
      # Cymbopogon present fire present
      Treatment == "bgpnf" ~ "CPFA",
      # Cymbopogon present fire absent
      Treatment == "bgrnf" ~ "CAFA" # Cymbopogon absent fire absent
    )
  )

# create Site prep data
Site_prep <- Site_dat %>%
  arrange(Site, Treatment) %>%
  mutate(
    Treatment = case_when(
      Treatment == "ab" ~ "Control",
      # Cymbopogon present fire present
      Treatment == "bgpnf" ~ "CPFA",
      # Cymbopogon present fire absent
      Treatment == "bgrnf" ~ "CAFA" # Cymbopogon absent fire absent
    )
  )

# what is the summed biomass per Site with cymbopogon?
Site_sum <-
  Site_prep %>% group_by(Site, Treatment) %>%
  summarise(Site_biomass = sum(Weight)) %>%
  ungroup()

# Site_calc
Site_calc <- Site_prep %>% left_join(Site_sum) %>%
  mutate(
    relative_biomass = (Weight / Site_biomass) ,
    relative_biomass_p =  round(((
      Weight / Site_biomass
    ) * 100) , 2)
  )

# Absolute_biomass
absolute_weight <- Site_calc %>%
  select(Site, Treatment, Weight, Palatability) %>%
  group_by(Palatability, Treatment, Site) %>%
  summarise(Weight = sum(Weight)) %>%
  mutate(Treatment = factor(Treatment)) %>% # to order treatments in the plot
  mutate(Palatability = factor(Palatability)) %>%
  mutate(Palatability = fct_relevel(Palatability, c('Yes', 'No', 'Cymbopogon sp.'))) %>%
  mutate(Treatment = fct_relevel(Treatment, c("Control", "CPFA", "CAFA")))

# number of graminoids (annual/perennial) and forbs (annual/perennial)
species_list <- Site_dat %>% 
  select(Scientific_name, Family, Treatment, Functional_type, Functional_groups, Palatability) %>%
  mutate(
    Functional_type= case_when(
      Functional_type == "Annual undershrub"  ~ "Annual",
      Functional_type == "Annual herb"   ~ "Annual",
      Functional_type == "Annual graminoid"  ~ "Annual",
      Functional_type == "Annual/perennial graminoid "  ~ "Annual",
      Functional_type == "Annual/perennial herb"  ~ "Annual",
      Functional_type == "Perennial herb"   ~ "Perennial",
      Functional_type == "Perennial graminoid"   ~ "Perennial",
      Functional_type == "Perennial undershrub"  ~ "Perennial")) %>%
  mutate(Palatability= recode(Palatability, 'Cymbopogon sp.'= 'No')) %>%  
  distinct(Functional_groups, Treatment, Palatability, Functional_type, Scientific_name, Family)

# write.csv(species_list, 'Species list-Fam_ supplementary- 1.csv') # add 4 Cymb. grass manually.
species_list %>%
  mutate(anu.peri = as.factor(Functional_groups)) %>%
  dplyr::count(Functional_groups, Functional_type, Palatability) 
# 56 palatable 13 non palatable- including 4 Cymbopogons as No palatable.

cafa <- species_list %>%
  filter(Treatment=='bgrnf') %>% 
  mutate(anu.peri = as.factor(Functional_groups)) %>%
  dplyr::count(Functional_type, Functional_groups,Palatability) %>% 
  mutate(Treatment= as.factor('CAFA'))

cpfa <- species_list %>%
  filter(Treatment=='bgpnf') %>% 
  mutate(anu.peri = as.factor(Functional_groups)) %>%
  dplyr::count(Functional_type, Functional_groups,Palatability ) %>% 
  mutate(Treatment= as.factor('CPFA'))

cpfp <- species_list %>%
  filter(Treatment=='ab') %>% 
  mutate(anu.peri = as.factor(Functional_groups)) %>%
  dplyr::count(Functional_type, Functional_groups,Palatability ) %>% 
  mutate(Treatment= as.factor('CPFP'))
