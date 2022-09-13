# Packages----
library(tidyverse)
library(patchwork)
library(vegan)
library(MetBrewer)
library(forcats)

# Data----
# raw data
raw_dat <- read.csv(
  "biomass_data.csv",
  header = T,
  fill = TRUE,
  sep = ",",
  na.strings = c("", " ", "NA", "NA ", "na", "NULL")
)

# Data wrangling----
transect_dat <- raw_dat %>%
  mutate(
    Treatment = as.factor(Treatment),
    Life_form = as.factor(Life_form),
    Functional_groups = as.factor(Functional_groups)
  )

# create transect prep data
transect_prep <- transect_dat %>%
  arrange(Transect, Treatment) %>%
  mutate(
    Treatment = case_when(
      Treatment == "ab" ~ "Control",
      # Cymbopogon present fire present
      Treatment == "bgpnf" ~ "CPFA",
      # Cymbopogon present fire absent
      Treatment == "bgrnf" ~ "CAFA" # Cymbopogon absent fire absent
    )
  )

# what is the summed biomass per transect with cymbopogon?
transect_sum <-
  transect_prep %>% group_by(Transect, Treatment) %>%
  summarise(transect_biomass = sum(Weight)) %>%
  ungroup()

# Transect_calc
transect_calc <- transect_prep %>% left_join(transect_sum) %>%
  mutate(
    relative_biomass = (Weight / transect_biomass) ,
    relative_biomass_p =  round(((
      Weight / transect_biomass
    ) * 100) , 2)
  )

# Absolute_biomass
absolute_weight <- transect_calc %>%
  select(Transect, Treatment, Weight, Palatability) %>%
  group_by(Palatability, Treatment, Transect) %>%
  summarise(Weight = sum(Weight)) %>%
  mutate(Treatment = factor(Treatment)) %>% # to order treatments in the plot
  mutate(Palatability = factor(Palatability)) %>%
  mutate(Palatability = fct_relevel(Palatability, c('Yes', 'No', 'Cymbopogon sp.'))) %>%
  mutate(Treatment = fct_relevel(Treatment, c("Control", "CPFA", "CAFA")))

head(transect_calc)

# Relative_biomass
relative_weight <-
  transect_calc %>% select( Sci_name, Site, Transect, Treatment, Palatability, relative_biomass) %>%
  group_by(Sci_name, Site, Treatment, Transect, Palatability) %>%
  summarise(relative_biomass = (relative_biomass) * 100) %>%
  mutate(Treatment = factor(Treatment)) %>%
  mutate(Palatability = factor(Palatability)) %>%
  mutate(Palatability = fct_relevel(Palatability, c('Yes', 'No', 'Cymbopogon sp.'))) %>%
  mutate(Treatment = fct_relevel(Treatment, c("Control", "CPFA", "CAFA"))) %>% 
  arrange(Transect) %>% ungroup()

head(relative_weight)

eghats_details <- relative_weight %>% select(Site, Treatment, Transect) %>% distinct()

eghats_species_details <- relative_weight %>% select(Sci_name, Palatability) %>% distinct()

# create a transect by species matrix
relative_weight_matrix <- relative_weight %>%  select(Transect, Sci_name, relative_biomass) %>%
  mutate(Sci_name = str_replace(Sci_name, " ", "_")) %>%
  group_by(Transect) %>%
  spread(Sci_name, relative_biomass) %>% replace(is.na(.), 0) %>%
  `row.names<-`(., NULL) %>% 
  column_to_rownames(var = "Transect")

head(relative_weight_matrix)

# run NMDS with  bray curtis distance
# https://chrischizinski.github.io/rstats/vegan-ggplot2/

eghats.mds <- metaMDS(relative_weight_matrix, distance = "bray", autotransform = FALSE)

plot(eghats.mds, type = "t")

# site scores
data.scores <- as.data.frame(scores(eghats.mds))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$Transect <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores

eghats.scores <- data.scores %>% 
  left_join(eghats_details) %>% arrange(Treatment)

head(eghats.scores)

# species scores
species.scores <- as.data.frame(scores(eghats.mds, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
head(species.scores)  #look at the data

eghats.species.scores <- species.scores %>% 
  mutate(Sci_name = species) %>% select(-species) %>%
  mutate(Sci_name = str_replace(Sci_name, "_", " ")) %>%
  left_join(eghats_species_details)

head(eghats.species.scores)

head(eghats.scores)

eghats.ctl <- eghats.scores[eghats.scores$Treatment == "Control", ][chull(eghats.scores[eghats.scores$Treatment == 
                                                                   "Control", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
eghats.cpfa <- eghats.scores[eghats.scores$Treatment == "CPFA", ][chull(eghats.scores[eghats.scores$Treatment == 
                                                                   "CPFA", c("NMDS1", "NMDS2")]), ]  # hull values for grp B
eghats.cafa<- eghats.scores[eghats.scores$Treatment == "CAFA", ][chull(eghats.scores[eghats.scores$Treatment == 
                                                                               "CAFA", c("NMDS1", "NMDS2")]), ]  # hull values for grp B

hull.data <- rbind(eghats.ctl, eghats.cpfa, eghats.cafa)  #combine groups
hull.data

nmdsplot <- ggplot() + 
  geom_polygon(data=hull.data, aes(x=NMDS1, y=NMDS2, fill=Treatment, group=Treatment), alpha=0.30) + # add the convex hulls
  #geom_text(data=eghats.species.scores, aes(x=NMDS1,y=NMDS2,label=Sci_name),alpha=0.5) +  # add the species labels
  geom_point(data=hull.data, aes(x=NMDS1, y=NMDS2, shape=Treatment, colour=Treatment),size=2) + # add the point markers
  scale_color_manual(values = c("Control" = "#BB9689",
                                "CPFA" = "#836656",
                                "CAFA" = "#6C3859"))+
  scale_fill_manual(values = c("Control" = "#BB9689",
                                "CPFA" = "#836656",
                                "CAFA" = "#6C3859"))+
  coord_equal() +
  theme_bw() + 
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=12), # remove x-axis labels
        axis.title.y = element_text(size=12), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())+
  theme(legend.position = 'top')

nmdsplot


# could do something like Michael's plots below his NMDS in Figure 4
#  https://doi.org/10.1111/rec.13006
# you already have somehting like this but I would suggest doing it for relative cover at the treatment scale so you have 3 different histograms below your NMDS
# one for each treatment

# Histogram of top xx biomass species

for_hist <-
  relative_weight %>% group_by(Treatment) %>% top_n(25) %>%
  arrange(Treatment, relative_biomass)

sp_hist <-
  ggplot(for_hist,
         aes(relative_biomass, Sci_name, col = Treatment, fill = Treatment)) +
  geom_histogram(stat = 'identity') +
  facet_wrap( ~ Treatment) +
  scale_color_manual(values = c(
    "Control" = "#BB9689",
    "CPFA" = "#836656",
    "CAFA" = "#6C3859"
  ))  +
  scale_fill_manual(values = c(
    "Control" = "#BB9689",
    "CPFA" = "#836656",
    "CAFA" = "#6C3859"
  ))+
  theme_bw()+
  theme(axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=12), # remove x-axis labels
        axis.title.y = element_text(size=12), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())+
  theme(legend.position = 'none')
  


sp_histogram <- sp_hist + aes(y = reorder(Sci_name, relative_biomass))+
  labs(x='Relative biomass', y= 'Scientific name')


sp_histogram+nmdsplot

# Save image (Evenness)
ggsave('fig_5.jpg', width = 10, height = 6, dpi = 300) 
