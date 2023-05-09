# Packages----
library(tidyverse)
library(patchwork)
library(vegan)

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

# Relative_biomass
relative_weight <-
  Site_calc %>% select( Sci_name, Village, Site, Treatment, Palatability, relative_biomass) %>%
  group_by(Sci_name, Village, Treatment, Site, Palatability) %>%
  summarise(relative_biomass = (relative_biomass) * 100) %>%
  mutate(Treatment = factor(Treatment)) %>%
  mutate(Palatability = factor(Palatability)) %>%
  mutate(Palatability = fct_relevel(Palatability, c('Yes', 'No', 'Cymbopogon sp.'))) %>%
  mutate(Treatment = fct_relevel(Treatment, c("Control", "CPFA", "CAFA"))) %>% 
  arrange(Site) %>% ungroup()

# df of top 10 species----
species_biomass <- Site_calc %>%
  select(Treatment, Weight, Sci_name) %>%
  group_by(Sci_name, Treatment,) %>%
  summarise(Weight = sum(Weight)) %>%
  mutate(Treatment = fct_relevel(Treatment, c("Control", "CPFA", "CAFA")))

treatment_total <-
  species_biomass %>% 
  group_by(Treatment) %>% 
  summarise(group_biomass= sum(Weight))

top10_biomass <- left_join(x= species_biomass, 
                           y = treatment_total, 
                           "Treatment") %>% 
  mutate(Treatment= as.factor(Treatment)) %>% 
  mutate(relative_biomass= (Weight/group_biomass)*100) %>% 
  group_by(Treatment) %>% top_n(10) %>%
  arrange(Treatment, relative_biomass) %>%
  mutate(Treatment = fct_relevel(Treatment, c("Control", "CPFA", "CAFA")))

# Making nmds metrics
eghats_details <- relative_weight %>% select(Village, Treatment, Site) %>% distinct()

eghats_species_details <- relative_weight %>% select(Sci_name, Palatability) %>% distinct()

# create a Site by species matrix
relative_weight_matrix <- relative_weight %>%  select(Site, Sci_name, relative_biomass) %>%
  mutate(Sci_name = str_replace(Sci_name, " ", "_")) %>%
  group_by(Site) %>%
  spread(Sci_name, relative_biomass) %>% replace(is.na(.), 0) %>%
  `row.names<-`(., NULL) %>% 
  column_to_rownames(var = "Site")


# run NMDS with  bray curtis distance
# https://chrischizinski.github.io/rstats/vegan-ggplot2/

set.seed(1980)
eghats.mds <- metaMDS(relative_weight_matrix, distance = "bray", autotransform = FALSE)
# stress is 0.18 ideal is less than 0.2

plot(eghats.mds, type = "t")

# permutation test----
relative_weight_matrix_1 <- relative_weight %>%  select(Site, Sci_name, relative_biomass) %>%
  mutate(Sci_name = str_replace(Sci_name, " ", "_")) %>%
  group_by(Site) %>%
  spread(Sci_name, relative_biomass) %>% replace(is.na(.), 0) %>%
  `row.names<-`(., NULL)

data <- eghats_details %>% left_join(relative_weight_matrix_1) %>% select(-Site)

dataM <- data %>% 
  select(-c(Village, Treatment)) %>% 
  as.matrix() 

rownames(dataM) <- data$Village

# Bray-Curtis
distances <- vegdist(dataM, 
                     distance = 'bray')
# order of samples
distances_groups <- data[match(labels(distances),
                               data$Village),]$Treatment 
# beta dispersion
distances_betadispersion <- betadisper(distances, 
                                       distances_groups)

distance <- data.frame(dist=distances_betadispersion$distances, group=distances_betadispersion$group)

boxplot(dist~group, data= distance)

anova(distances_betadispersion) # not met, p < 0.05

# test homogenous dispersion/Permutation test for F
pmod <- permutest(distances_betadispersion, permutations = 99, pairwise = TRUE)

## Tukey's Honest Significant Differences
(mod.HSD <- TukeyHSD(distances_betadispersion))
plot(mod.HSD)

# Village scores for ggplot
data.scores <- as.data.frame(scores(eghats.mds)$sites)  #Using the scores function from vegan to extract the Village scores and convert to a data.frame
data.scores$Site <- rownames(data.scores)  # create a column of Village names, from the rownames of data.scores

eghats.scores <- data.scores %>% 
  left_join(eghats_details) %>% arrange(Treatment)

# species scores
species.scores <- as.data.frame(scores(eghats.mds, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores

eghats.species.scores <- species.scores %>% 
  mutate(Sci_name = species) %>% select(-species) %>%
  mutate(Sci_name = str_replace(Sci_name, "_", " ")) %>%
  left_join(eghats_species_details)

# convex hull values for treatments
eghats.ctl <-
  eghats.scores[eghats.scores$Treatment == "Control", ][chull
                                                        (eghats.scores[eghats.scores$Treatment ==
                                                                         "Control",
                                                                       c("NMDS1", "NMDS2")]), ]
eghats.cpfa <-
  eghats.scores[eghats.scores$Treatment == "CPFA", ][chull(eghats.scores[eghats.scores$Treatment ==
                                                                           "CPFA", c("NMDS1", "NMDS2")]), ]
eghats.cafa <-
  eghats.scores[eghats.scores$Treatment == "CAFA", ][chull(eghats.scores[eghats.scores$Treatment ==
                                                                           "CAFA", c("NMDS1", "NMDS2")]), ]

hull.data <-
  rbind(eghats.ctl, eghats.cpfa, eghats.cafa)  #combine groups
hull.data

# nmdsplot----
treatment_colors <- c("Control" = "#3b5d4d", "CPFA" = "#c5af99","CAFA" = "#ffd365")

nmdsplot <- ggplot() +
  geom_point(
    data = hull.data,
    aes(
      x = NMDS1,
      y = NMDS2,
      shape = Treatment,
      colour = Treatment
    ),
    size = 2
  ) + # add the point markers
  geom_polygon(
    data = hull.data,
    aes(
      x = NMDS1,
      y = NMDS2,
      fill = Treatment,
      group = Treatment
    ),
    alpha = 0.30
  ) + # add the convex hulls
  #geom_text(data=eghats.species.scores, aes(x=NMDS1,y=NMDS2,label=Sci_name),alpha=0.5) +  # add the species labels
  scale_color_manual(values = treatment_colors) +
  scale_fill_manual(values = treatment_colors ) +
  coord_equal() +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    # remove x-axis text
    axis.text.y = element_blank(),
    # remove y-axis text
    axis.ticks = element_blank(),
    # remove axis ticks
    axis.title.x = element_text(size = 12),
    # remove x-axis labels
    axis.title.y = element_text(size = 12),
    # remove y-axis labels
    panel.grid.major = element_blank(),
    #remove major-grid labels
    panel.grid.minor = element_blank(),
    #remove minor-grid labels
  ) +
  theme(legend.position = 'right')

nmdsplot+
  annotate(
    geom = "text",
    x = -1.03,
    y = 0.5,
    size = 3,
    label = "Cymbopogon spp.",
    color = "red",
    fontface = 'italic'
  )+ 
  annotate(
    geom = "text",
    x = -0.4,
    y = 1.1,
    size = 3,
    label = "A. nervosa",
    color = "black",
    fontface = 'italic'
  )+ 
  annotate(
    geom = "text",
    x = 0.16,
    y = 1.5,
    size = 3,
    label = "A. nervosa",
    color = "black",
    fontface = 'italic'
  )+ 
  annotate(
    geom = "text",
    x = 0.6,
    y = 1.65,
    size = 3,
    label = "P. wightiana",
    color = "black",
    fontface = 'italic'
  )+ 
  annotate(
    geom = "text",
    x = -0.3,
    y = -0.85,
    size = 3,
    label = "A. mutica",
    color = "black",
    fontface = 'italic'
  )+ 
  annotate(
    geom = "text",
    x = 0.4,
    y = -1.05,
    size = 3,
    label = "H. contortus",
    color = "black",
    fontface = 'italic'
  )+ 
  annotate(
    geom = "text",
    x = 0.9,
    y = -1.04,
    size = 3,
    label = "H. contortus",
    color = "black",
    fontface = 'italic'
  )+ 
  annotate(
    geom = "text",
    x = 1.01,
    y = 0.5,
    size = 3,
    label = "L. cristata",
    color = "black",
    fontface = 'italic'
  )

# Save image----
ggsave('fig_7_nmds.jpg', width = 10, height = 6, dpi = 300) 

top10_hist <- ggplot(top10_biomass, aes(relative_biomass, Sci_name, fill=Treatment))+
  geom_histogram(stat = 'identity')+
  facet_grid(~Treatment)+
  scale_color_manual(values = treatment_colors)  +
  scale_fill_manual(values = treatment_colors)+
  theme_bw()+
  theme(axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=12), # remove x-axis labels
        axis.title.y = element_text(size=12), # remove y-axis labels
        axis.text.y= element_text(size = 11 ,face = 'italic'),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank())+  #remove minor-grid labels
  theme(legend.position = 'none')+
  aes(y = reorder(Sci_name, relative_biomass))+
  labs(x=' Relative biomass', y= 'Scientific names')
top10_hist
# ggsave('fig_7_topten.jpg', width = 10, height = 6, dpi = 300)

# For my own understanding:
# In this ordination, the closer two points are, the more similar the corresponding samples are with respect to the variables that went into making the NMDS plot.
# The closer the points/samples are together in the ordination space, the more similar their plant communities.
# In the histogram we can see the species responsible for the result in NMDS
# We can see the names of top 10 species responsible for this result

# statistical analysis----
# ? is composition different across our treatments- may be not because there is 
# a lot of overlap.

