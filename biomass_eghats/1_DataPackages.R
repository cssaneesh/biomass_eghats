# Packages----
rm(list = ls())
library(bayesplot)
library(brms)
library(cmdstanr)
library(cowplot)
library(gt)
library(patchwork)
library(tidyverse)
library(vegan)
library(webshot2)
library(gridExtra)

# Data----
data <- read.csv(
  "biomass_data.csv",
  header = T,
  fill = TRUE,
  sep = ",",
  na.strings = c("", " ", "NA", "NA ", "na", "NULL")
)

# Data wrangling----
data <- data %>%
  mutate(
    Treatment = as.factor(Treatment),
    Life_form = as.factor(Life_form),
    Functional_groups = as.factor(Functional_groups)
  )

# create Site prep data
Site_prep <- data %>%
  arrange(Site, Treatment) %>%
  mutate(
    Treatment = case_when(
      Treatment == "ab" ~ "Control",
      # Cymbopogon present fire present
      Treatment == "bgpnf" ~ "CPFA",
      # Cymbopogon present fire absent
      Treatment == "bgrnf" ~ "CAFA"  # Cymbopogon absent fire absent
    )
  )

# what is the summed biomass per Site without cymbopogon?
Site_sum_no_cym <- Site_prep %>%
  filter(!Sci_name == "Cymbopogon sp.") %>%
  group_by(Site, Treatment) %>%
  summarise(Site_biomass_no_cym = sum(Weight)) %>%
  ungroup()

# what is the summed biomass per Site with cymbopogon?
Site_sum_w_cym <-
  Site_prep %>% group_by(Site, Treatment) %>%
  summarise(Site_biomass = sum(Weight)) %>%
  ungroup()

# biomass of non cymbopogons and cymbopogons
Site_sum <-
  Site_sum_w_cym %>% left_join(Site_sum_no_cym)

# Site_calc
Site_calc <- Site_prep %>% left_join(Site_sum) %>%
  mutate(
    relative_biomass = (Weight / Site_biomass) ,
    relative_biomass_p = ((Weight / Site_biomass) * 100),
    relative_biomass_nc = (Weight / Site_biomass_no_cym) ,
    relative_biomass_nc_p = ((Weight / Site_biomass) * 100)
  )

# Absolute_biomass----
absolute_weight <- Site_calc %>%
  select(Site, Treatment, Weight, Palatability) %>%
  group_by(Palatability, Treatment, Site) %>%
  summarise(Weight = sum(Weight)) %>%
  mutate(Treatment = factor(Treatment)) %>% # to order treatments in the plot
  mutate(Palatability = factor(Palatability)) %>%
  mutate(Palatability = fct_relevel(Palatability, c('Yes', 'No', 'Cymbopogon sp.'))) %>%
  mutate(Treatment = fct_relevel(Treatment, c("Control", "CPFA", "CAFA")))

# Relative_biomass----
relative_weight <-
  Site_calc %>% select(Village, Site, Treatment, Palatability, relative_biomass) %>%
  group_by(Village, Treatment, Site, Palatability) %>%
  summarise(relative_biomass = sum(relative_biomass) * 100) %>%
  mutate(Treatment = factor(Treatment)) %>%
  mutate(Palatability = factor(Palatability)) %>%
  mutate(Palatability = fct_relevel(Palatability, c('Yes', 'No', 'Cymbopogon sp.'))) %>%
  mutate(Treatment = fct_relevel(Treatment, c("Control", "CPFA", "CAFA")))

# Fire history of sites-----
# EHM (biomass_trees) sites of FES from 2017 to 2019 Anantapur and Chittoor without missing data 
# and continuously monitored

fire_his <-
  read.csv('fire_interval.csv')

# names(fire_his)

# average fire return interval = (total number of years)/(total number of fires)


afri <- fire_his %>% # afri= average fire return interval
  group_by(Habitation, Plot_no) %>% # now make a logical vector
  summarize(
    Burnt_once = sum(Fire == "Yes" &
                       Year %in% c(2017, 2018, 2019)) == 1,
    # in these years look for 1 'yes'
    Burnt_twice = sum(Fire == "Yes" &
                        Year %in% c(2017, 2018, 2019)) == 2,
    Burnt_thrice = sum(Fire == "Yes" &
                         Year %in% c(2017, 2018, 2019)) == 3
  ) %>%
  mutate(
    # convert logical response to numeric
    Burnt_once = ifelse(Burnt_once, 1, 0),
    Burnt_twice = ifelse(Burnt_twice, 1, 0),
    Burnt_thrice = ifelse(Burnt_thrice, 1, 0)
  ) %>%
  summarise(
    # summaries data
    once = sum(Burnt_once),
    twice = sum(Burnt_twice),
    thrice = sum(Burnt_thrice)
  ) %>%
  summarise(across(where(is.numeric), # remove any missing values before calculating the sum
                   ~ sum(.x, na.rm = TRUE)))
afri %>% mutate(
  years= ncol(afri),
  afri= years/rowSums(.))

# Family figure----
data <- data.frame(
  category=c("Poaceae",
             'Fabaceae', 
             "Asteraceae", 
             "Acanthaceae", 
             'Rubiaceae', 
             'Euphorbiaceae', 
             'Others'),
  count=c(18,15,8,4,4,3,17 )
)

data$fraction = data$count / sum(data$count)
data$ymax = cumsum(data$fraction)
data$ymin = c(0, head(data$ymax, n=-1))

# Compute label position
data$labelPosition <- (data$ymax + data$ymin) / 2

# Compute a good label
data$label <- paste0(data$category, "\n Species: ", data$count)

families.pie <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect(col='darkgrey', linetype= 'dotted')+
  scale_fill_brewer(palette=4) +
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(legend.position = "none")+
  geom_text( x=2, aes(y=labelPosition, label=label),
             #color='black',
             color='brown',
             #color= c('#2b8cbe', '#7bccc4', '#ccebc5','#f0f9e8','#08589e', '#a8ddb5', '#4eb3d3'), 
             size=4)
families.pie

ggsave('families.jpg',
       width = 10,
       height = 6,
       dpi = 300)



