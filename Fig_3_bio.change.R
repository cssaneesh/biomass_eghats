# Packages----
library(tidyverse)
library(brms)
library(gt)

# Data----
raw_dat <- read.csv(
  "biomass_data.csv",
  header = T,
  fill = TRUE,
  sep = ",",
  na.strings = c("", " ", "NA", "NA ", "na", "NULL")
)

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

# Absolute weight----
bio.change_prep <- Site_prep %>% select(Village, Treatment, Sci_name, Palatability, Weight) %>%
  group_by(Village) %>% 
  nest(Treat.nest=c(Treatment, Sci_name, Palatability, Weight)) %>% 
  ungroup()

bio.change.nest <- bio.change_prep %>%
  mutate(Treat.nest = purrr::map(Treat.nest, ~ spread(., key = Treatment, value = Weight))) %>%
  mutate(Treat.nest = purrr::map(Treat.nest, ~ replace(., is.na(.), 0))) %>%
  unnest(Treat.nest)

head(bio.change.nest, 3)

bio.change <- bio.change.nest %>% 
  mutate(B.m.CAFA= (CAFA-Control), B.m.CPFA= (CPFA-Control))

head(bio.change, 3)

bio.change.treatment <- bio.change %>% select(-c(CAFA, CPFA, Control)) %>%
  filter(!is.na(B.m.CAFA)|!is.na(B.m.CPFA)) %>% 
  gather(Treatment, Biomass, B.m.CAFA:B.m.CPFA)

head(bio.change.treatment)  

bio.change.treatment %>% select(Village) %>% 
  distinct() # 14 villages

# Analysis----
ghats.biomass.change <-
  brm(
    Biomass ~   Treatment * Palatability+
      ( 1 | Village/Sci_name ) ,
    family = student(),
    data = bio.change.treatment,
    iter = 10000,
    warmup = 1000,
    cores = 4,
    chains = 4,
    control = list(adapt_delta = 0.99, max_treedepth = 12 )
  )

# save(ghats.biomass.change, file = "ghats.biomass.change.Rdata")

load("ghats.biomass.change.Rdata")

summary(ghats.biomass.change)

fig_s4 <- pp_check(ghats.biomass.change)+
  xlab("Relative biomass change") + ylab("Density") +
  labs(title = "Site-level",
       #subtitle = "a)"
  ) +
  theme_classic() + xlim(-40,150) +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        legend.position = "none")# predicted vs. observed values

fig_s4

ggsave('sup_Fig_4.jpg',
       width = 10,
       height = 6,
       dpi = 300)

plot(ghats.biomass.change)

ghats_change_biomass <-
  conditional_effects(
    ghats.biomass.change,
    effects = 'Treatment:Palatability',
    re_formula = NA,
    method = 'fitted'
  )  # conditional effects
ghats_change_biomass

ghats.biomass.change_df <-
  as.data.frame(ghats_change_biomass$`Treatment:Palatability`)

ghats.biomass.change_df
#View(ghats.biomass.change_df)

# ghats.biomass.change_df
# names(ghats.biomass.change_df)
# View(ghats.biomass.change_df)

# names(bio.change.treatment)
# Plot----
fig_change.biomass <- ggplot() +
  geom_point(
    data = bio.change.treatment,
    aes(
      x = Treatment,
      y = Biomass,
      group = Palatability,
      colour = 	Palatability
    ),
    size = 1,
    alpha = 0.8,
    position = position_jitterdodge(
      jitter.width = 0.05,
      jitter.height = 0.45,
      dodge.width = 0.75,
      seed = NA
    ) # jitter points are species
  ) +
  geom_point(
    data = ghats_change_biomass$`Treatment:Palatability`,
    aes(
      x = Treatment,
      y = estimate__, # mean value from of the model output
      group = Palatability),
    size = 2,
    alpha= 0.6,
    position = position_dodge(width = 0.75,)
  )+
  geom_errorbar(
    data = ghats_change_biomass$`Treatment:Palatability`,
    aes(
      x = Treatment,
      ymin = lower__, # lowest value from of the model output
      ymax = upper__, # highest value from of the model output
      group = Palatability),
    position = position_dodge(width = 0.75),
    linewidth = 0.4,
    width = 0.09,
    alpha= 0.6
  ) + labs(x = '', y = '') +
  scale_color_manual(values =  c("#1a6Ba8", '#b57d70', "#9a9a9a"))+
  theme_bw(base_size = 12) + theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(
      t = 0.2,
      r = 0.2,
      b = 0.2,
      l = 0.2,
      unit = "cm"
    ),
    plot.title = element_text(size = 14, hjust =
                                0.5),
    strip.background = element_blank(),
    #legend.position = "bottom"
  ) + labs(subtitle = '') + ylab("Biomass change (g)") +
  theme(
    panel.grid.major = element_line(colour = "gray86", size = 0.1),
    panel.background = element_rect(fill = "white")
  )+
  theme(axis.ticks = element_blank())+
  geom_hline(yintercept = 0, linetype= 'dashed', color='darkgrey')+
  ylim(-1000, 1000)

fig_change.biomass+ 
  annotate(geom="text", 
           x=1.02, 
           y=1000, 
           size= 3,
           label="Apluda mutica",
           color="black")+
  annotate(
    geom = 'segment',
    x=1.24, xend= 1.07,
    y= 860, yend= 960)+
  annotate(geom="text", 
           x=1.02, 
           y=-500, 
           size= 3,
           label="Arundinella nervosa",
           color="black")+
  annotate(
    geom = 'segment',
    x=1.24, xend= 1.07,
    y= -300, yend= -470)


# library(plotly)
# ggplotly(fig_change.biomass)

# Relative weight----
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
  Site_calc %>% select(Village, Treatment, Sci_name, Palatability, relative_biomass) %>%
  group_by(Village, Treatment, Sci_name, Palatability) %>%
  summarise(relative_biomass = sum(relative_biomass) * 100) %>%
  mutate(Treatment = factor(Treatment)) %>%
  mutate(Palatability = factor(Palatability)) %>%
  mutate(Palatability = fct_relevel(Palatability, c('Yes', 'No', 'Cymbopogon sp.'))) %>%
  mutate(Treatment = fct_relevel(Treatment, c("Control", "CPFA", "CAFA")))

rel.bio.change_prep <- relative_weight %>% select(Village, Treatment, Sci_name, Palatability, relative_biomass) %>%
  group_by(Village) %>% 
  nest(Treat.nest=c(Treatment, Sci_name, Palatability, relative_biomass)) %>% 
  ungroup()

rel.bio.change.nest <- rel.bio.change_prep %>%
  mutate(Treat.nest = purrr::map(Treat.nest, ~ spread(., key = Treatment, value = relative_biomass))) %>%
  mutate(Treat.nest = purrr::map(Treat.nest, ~ replace(., is.na(.), 0))) %>%
  unnest(Treat.nest)

head(rel.bio.change.nest, 3)

rel.bio.change <- rel.bio.change.nest %>% 
  mutate(B.m.CAFA= (CAFA-Control), B.m.CPFA= (CPFA-Control))
# View(rel.bio.change)
head(rel.bio.change, 3)

rel.bio.change.treatment <- rel.bio.change %>% select(-c(CAFA, Control, CPFA)) %>%
  filter(!is.na(B.m.CAFA)|!is.na(B.m.CPFA)) %>% 
  gather(Treatment, Rel.Biomass, B.m.CAFA:B.m.CPFA) %>% 
  mutate(Treatment=recode(Treatment, 'B.m.CAFA'= 'CAFA',
                          'B.m.CPFA'='CPFA')) %>%
  mutate(Treatment = fct_relevel(Treatment, c("CPFA", "CAFA")))

head(rel.bio.change.treatment)  

# Analysis----
# rel.ghats.biomass.change <-
#   brm(
#     Rel.Biomass ~   Treatment * Palatability+
#       ( 1 | Village/Sci_name ) ,
#     family = student(),
#     data = rel.bio.change.treatment,
#     iter = 3000,
#     warmup = 1000,
#     cores = 4,
#     chains = 4
#   )

# save(rel.ghats.biomass.change, file = "rel.ghats.biomass.change.Rdata")

load("rel.ghats.biomass.change.Rdata")

summary(rel.ghats.biomass.change)
pp_check(rel.ghats.biomass.change)
plot(rel.ghats.biomass.change)

rel.ghats_change_biomass <-
  conditional_effects(
    rel.ghats.biomass.change,
    effects = 'Treatment:Palatability',
    re_formula = NA,
    method = 'fitted'
  )  # conditional effects

rel.ghats.biomass.change_df <-
  as.data.frame(rel.ghats_change_biomass$`Treatment:Palatability`)

# View(rel.ghats.biomass.change_df)

rel.ghats.biomass.change_df
names(rel.ghats.biomass.change_df)
# View(ghats.biomass.change_df)

names(rel.bio.change.treatment)
# Plot----
set.seed(1980)

fig_rel.change.biomass <- ggplot() +
  geom_point(
    data = rel.bio.change.treatment,
    aes(
      x = Treatment,
      y = Rel.Biomass,
      group = Palatability,
      colour = 	Palatability
    ),
    size = 1,
    alpha = 0.7,
    position = position_jitterdodge(
      jitter.width = 0.09,
      jitter.height = 0.45,
      dodge.width = 0.75,
      seed = NA
    ) # jitter points are species
  ) +
  geom_point(
    data = rel.ghats_change_biomass$`Treatment:Palatability`,
    aes(
      x = Treatment,
      y = estimate__, # mean value from of the model output
      group = Palatability,
      colour = Palatability
    ),
    size = 3,
    alpha= 0.6,
    position = position_dodge(width = 0.75,)
  ) +
  geom_errorbar(
    data = rel.ghats_change_biomass$`Treatment:Palatability`,
    aes(
      x = Treatment,
      ymin = lower__, # lowest value from of the model output
      ymax = upper__, # highest value from of the model output
      group = Palatability,
      #colour = Palatability
    ),
    position = position_dodge(width = 0.75),
    linewidth = 0.8,
    width = 0.1,
    alpha= 0.8
  ) + labs(x = '', y = '') +
  scale_color_manual(values =  c("#1a6Ba8", '#b57d70', "#9a9a9a"))+
  theme_bw(base_size = 12) + theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(
      t = 0.2,
      r = 0.2,
      b = 0.2,
      l = 0.2,
      unit = "cm"
    ),
    plot.title = element_text(size = 14, hjust =
                                0.5),
    strip.background = element_blank(),
    #legend.position = "bottom"
  ) + labs(subtitle = '') + ylab("Relative biomass change (g)") +
  theme(
    panel.grid.major = element_line(colour = "gray86", size = 0.1),
    panel.background = element_rect(fill = "white")
  )+
  theme(axis.ticks = element_blank())+
  geom_hline(yintercept = 0, linetype= 'dashed', color='darkgrey')+
  ylim(-70, 70)

fig_rel.change.biomass 

fig_rel.change.biomass <-
  fig_rel.change.biomass +
  annotate(
    geom = "text",
    x = 0.8,
    y = 28,
    size = 3,
    label = "C. fulvus",
    color = "black",
    fontface = 'italic'
  ) +
  annotate(
    geom = 'segment',
    x = 0.76,
    xend = 0.78,
    y = 20,
    yend = 24.5,
    size= 0.3
  ) + # cpfa 1.81 yes
  annotate(
    geom = "text",
    x = 0.8,
    y = -37,
    size = 3,
    label = "C. fulvus",
    color = "black",
    fontface = 'italic'
  ) +
  annotate(
    geom = 'segment',
    x =0.77,
    xend = 0.81,
    y = -29,
    yend = -32,
    size= 0.3
  ) + # cpfa -2.71 yes
  annotate(
    geom = "text",
    x = 0.93,
    y = 10,
    size = 3,
    label = "L. zeylanica",
    color = "black",
    fontface = 'italic'
  ) +
  annotate(
    geom = 'segment',
    x = .98,
    xend = .95,
    y = 3,
    yend = 7,
    size= 0.3
  )+ # cpfa no 1.53
  annotate(
    geom = "text",
    x = 1,
    y = -18,
    size = 3,
    label = "L. zeylanica",
    color = "black",
    fontface = 'italic'
  ) +
  annotate(
    geom = 'segment',
    x = 1.01,
    xend = 1,
    y = -10,
    yend = -15,
    size= 0.3
  )+ # cpfa no -7.72
  annotate(
    geom = "text",
    x = 1.65,
    y = 58,
    size = 3,
    label = "H. contortus",
    color = "black",
    fontface = 'italic'
  ) +
  annotate(
    geom = 'segment',
    x = 1.74,
    xend = 1.65,
    y = 68,
    yend = 60,
    size= 0.3
  ) + # cafa 7.47 yes
  annotate(
    geom = "text",
    x = 1.65,
    y = -57,
    size = 3,
    label = "S. pusilla",
    color = "black",
    fontface = 'italic'
  ) +
  annotate(
    geom = 'segment',
    x = 1.73,
    xend = 1.67,
    y = -52,
    yend = -55,
    size= 0.3
  ) + # cafa -5.03 yes
  annotate(
    geom = "text",
    x = 2.03,
    y = 25.5,
    size = 3,
    label = "T. villosa",
    color = "black",
    fontface = 'italic'
  ) +
  annotate(
    geom = 'segment',
    x = 2.01,
    xend = 2.03,
    y = 19,
    yend = 22,
    size= 0.3
  ) + # cafa no 1.74
  annotate(
    geom = "text",
    x = 2.01,
    y = -14,
    size = 3,
    label = "L. zeylanica",
    color = "black",
    fontface = 'italic'
  ) +
  annotate(
    geom = 'segment',
    x = 2.01,
    xend = 2.03,
    y = -9,
    yend = -11,
    size= 0.3
  ) # cafa no -7.72
  

fig_rel.change.biomass

ggsave('Fig_3_biochng.jpg',
       width = 10,
       height = 6,
       dpi = 300)


# plotly::ggplotly(fig_rel.change.biomass)


rel.ghats_biochange_df <- as.data.frame(rel.ghats_change_biomass$`Treatment:Palatability`)
# View(rel.ghats_biochange_df)

table_re.change <- rel.ghats_biochange_df %>% 
  select(Treatment, Palatability, estimate__, lower__, upper__) %>% 
  rename(Estimate= estimate__, 
         Lower= lower__,
         Upper= upper__) %>% 
  mutate(Estimate = round(Estimate, 2)) %>% 
  mutate(Lower = round(Lower, 2)) %>%
  mutate(Upper = round(Upper, 2)) %>%
  gt() %>% tab_options(column_labels.font.size = 11,
                       table.font.size = 10,
                       column_labels.font.weight = "bold") %>% 
  opt_table_font(default_fonts()) %>%  # Fonts: Roboto Mono,IBM Plex Mono, Red Hat Mono
  opt_table_outline(style = "solid", width = px(2))

table_re.change
table_re.change %>% gtsave('table_3_biochange.png', expand = 5)

# list of species 
rel.bio.change.treatment %>% drop_na() %>% 
  filter(Palatability=='Yes') %>% filter(Treatment=='CAFA') %>% View()

rel.bio.change.treatment %>% drop_na() %>% 
  filter(Palatability=='No') %>% filter(Treatment=='CAFA') %>% View()

rel.bio.change.treatment %>% drop_na() %>% 
  filter(Palatability=='Yes') %>% filter(Treatment=='CPFA') %>% View()

rel.bio.change.treatment %>% drop_na() %>% 
  filter(Palatability=='No') %>% filter(Treatment=='CPFA') %>% View()

# Winners and losers----
# winners of palatable species in cafa

rel.bio.change.treatment %>% drop_na() %>% 
  filter(Palatability=='Yes') %>% filter(Treatment=='CAFA') %>% 
  filter(Rel.Biomass>0) %>% distinct(Sci_name) # 40 species 

# losers of palatable species in cafa
rel.bio.change.treatment %>% drop_na() %>% 
  filter(Palatability=='Yes') %>% filter(Treatment=='CAFA') %>% 
  filter(Rel.Biomass<0) %>% distinct(Sci_name) # 25 species 

# winners of unpalatable species in cafa
rel.bio.change.treatment %>% drop_na() %>% 
  filter(Palatability=='No') %>% filter(Treatment=='CAFA') %>% 
  filter(Rel.Biomass>0) %>% distinct(Sci_name) # 7 species 

# losers of unpalatable species in cafa
rel.bio.change.treatment %>% drop_na() %>% 
  filter(Palatability=='No') %>% filter(Treatment=='CAFA') %>% 
  filter(Rel.Biomass<0) %>% distinct(Sci_name) # 4 species

# winners of palatable species in cpfa
rel.bio.change.treatment %>% drop_na() %>% 
  filter(Palatability=='Yes') %>% filter(Treatment=='CPFA') %>% 
  filter(Rel.Biomass>0) %>% distinct(Sci_name) # 13 species 

# losers of palatable species in cpfa
rel.bio.change.treatment %>% drop_na() %>% 
  filter(Palatability=='Yes') %>% filter(Treatment=='CPFA') %>% 
  filter(Rel.Biomass<0) %>% distinct(Sci_name) # 13 species 

# winners of unpalatable species in cpfa
rel.bio.change.treatment %>% drop_na() %>% 
  filter(Palatability=='No') %>% filter(Treatment=='CPFA') %>% 
  filter(Rel.Biomass>0) %>% distinct(Sci_name) # 1 species 

# losers of unpalatable species in cpfa
rel.bio.change.treatment %>% drop_na() %>% 
  filter(Palatability=='No') %>% filter(Treatment=='CPFA') %>% 
  filter(Rel.Biomass<0) %>% distinct(Sci_name)# 3 species 
