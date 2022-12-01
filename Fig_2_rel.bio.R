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
# raw data
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


allsites <- bind_rows(cafa, cpfa, cpfp) %>% 
  mutate(Functional_groups=recode(Functional_groups,'Cymbopogon'='Graminoid'))

sup_fig_2 <- ggplot(allsites, aes(Functional_type,n, fill= Functional_groups))+
  geom_col()+
  facet_wrap(~Treatment)+
  labs(y= 'Number of species', x= 'Functional type') +
  scale_fill_manual (values =  c('#bfbfbf', "#5a9bd5"))+
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
  ) +
  theme(
    panel.grid.major = element_line(colour = "gray86", size = 0.1),
    panel.background = element_rect(fill = "white")
  )

ggsave('sup_fig_2 .jpg',
      width = 10,
      height = 6,
      dpi = 300)



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
  Site_calc %>% select(Village, Site, Treatment, Palatability, relative_biomass) %>%
  group_by(Village, Treatment, Site, Palatability) %>%
  summarise(relative_biomass = sum(relative_biomass) * 100) %>%
  mutate(Treatment = factor(Treatment)) %>%
  mutate(Palatability = factor(Palatability)) %>%
  mutate(Palatability = fct_relevel(Palatability, c('Yes', 'No', 'Cymbopogon sp.'))) %>%
  mutate(Treatment = fct_relevel(Treatment, c("Control", "CPFA", "CAFA")))

# Exploratory plots----
# Absolute weight
fig_e1 <-
  ggplot(data = absolute_weight) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(
    aes(
      x = Treatment,
      y = Weight,
      group = Palatability,
      fill = Palatability
    ),
    size = 0.75,
    position = position_jitterdodge(
      jitter.width = 0.05,
      jitter.height = 0,
      dodge.width = 0.75,
      seed = NA
    )
  ) +
  geom_boxplot(aes(x = Treatment,
                   y = Weight,
                   fill = Palatability),
               position = position_dodge(width = .80)) +
  theme_bw(base_size = 18) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(
      t = 0.2,
      r = 0.2,
      b = 0.2,
      l = 0.2,
      unit = "cm"
    ),
    plot.title = element_text(size = 18, hjust = 0.5),
    strip.background = element_blank(),
    legend.position = "none"
  ) +
  ylim(-10, 800) +
  ylab("Biomass weight (g)") + theme(legend.direction = "none") +
  xlab("") +
  labs(title = "Absolute biomass",
       subtitle = 'a)') +
  scale_fill_manual(values = c("#f8bb49", '#a45200', "#4A2300")) +
  theme(plot.caption = element_text(size = 8, face = "italic",
                                    hjust = 0))

# relative biomass
fig_e2 <-
  ggplot(data = relative_weight) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(
    aes(
      x = Treatment,
      y = relative_biomass,
      group = Palatability,
      fill = Palatability
    ),
    size = 0.75,
    position = position_jitterdodge(
      jitter.width = 0.05,
      jitter.height = 0,
      dodge.width = 0.75,
      seed = NA
    )
  ) +
  geom_boxplot(aes(x = Treatment,
                   y = relative_biomass,
                   fill = Palatability)) +
  theme_bw(base_size = 18) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  ) +
  labs(
    title = 'Relative biomass',
    subtitle = 'b)',
    x = 'Treatment',
    y = 'Relative biomass',
    caption = "Control= Cymbopogon present and fire present
CPFA= Cymbopogon present and fire absent
CAFA= Cymbopogon absent and fire absent"
  ) +
  scale_fill_manual(values = c("#f8bb49", '#a45200', "#4A2300")) +
  theme(plot.caption = element_text(size = 8, face = "italic",
                                    hjust = 0)) + theme(plot.title = element_text(hjust = 0.5,
                                                                                  vjust = 0))

fig_e <- (fig_e1 / fig_e2) # use patchwork to stick plots together

fig_e

head(relative_weight)

# Analysis----
# ghats.rel_biomass <-
#   brm(
#     relative_biomass ~   Treatment * Palatability  +
#       ( 1 | Site ),
#     family = gaussian(),
#      # family = student(),
#     # family = lognormal(),
#     data = relative_weight,
#     iter = 2000,
#     warmup = 1000,
#     cores = 4,
#     chains = 4
#   )
# 
# save(ghats.rel_biomass, file = "ghats.rel_biomass.Rdata")

load("ghats.rel_biomass.Rdata")

color_scheme_set("darkgray")

fig_s1a <- pp_check(ghats.rel_biomass) +
  xlab("Functional group Relative biomass") + ylab("Density") +
  labs(title = "Site-level",
    subtitle = "a)") +
  theme_classic() + xlim(-20,150) +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        legend.position = "none")# predicted vs. observed values

fig_s1a

summary(ghats.rel_biomass)

ghats_rel_biomass <-
  conditional_effects(
    ghats.rel_biomass,
    effects = 'Treatment:Palatability',
    re_formula = NA,
    method = 'fitted'
  )  # conditional effects

ghats_rel_biomass

ghats_rel_biomass_df <-
  as.data.frame(ghats_rel_biomass$`Treatment:Palatability`)# to make a df to report statistics results

# View(ghats_rel_biomass_df)

# Table_1 ----
# Biomass
table_2 <- ghats_rel_biomass_df %>%
  select(Treatment, Palatability, estimate__, lower__, upper__) %>%
  rename(Estimate = estimate__,
         Lower = lower__,
         Upper = upper__) %>%
  mutate_if(is.numeric, round, 2) %>% 
  gt() %>% 
  tab_options(column_labels.font.size = 11,
              table.font.size = 10,
              column_labels.font.weight = "bold") %>% 
  opt_table_font(default_fonts()) %>%  # Fonts: Roboto Mono,IBM Plex Mono, Red Hat Mono
  opt_table_outline(style = "solid", width = px(2))

table_2
table_2 %>% gtsave('Table_2_rel_bio.png', expand = 5) # expand to set white space

# Plot----
fig_rel_biomass <- ggplot() +
  geom_point(
    data = relative_weight,
    aes(
      x = Treatment,
      y = relative_biomass,
      group = Palatability,
      colour = 	Palatability
    ),
    size = 1,
    alpha = 0.7,
    position = position_jitterdodge(
      jitter.width = 0.05,
      jitter.height = 0.45,
      dodge.width = 0.75,
      seed = NA
    )
  ) +
  geom_point(
    data = ghats_rel_biomass$`Treatment:Palatability`,
    aes(
      x = Treatment,
      y = estimate__,
      group = Palatability,
      colour = Palatability
    ),
    size = 3,
    position = position_dodge(width = 0.75,)
  ) +
  geom_errorbar(
    data = ghats_rel_biomass$`Treatment:Palatability`,
    aes(
      x = Treatment,
      ymin = lower__,
      ymax = upper__,
      group = Palatability,
      colour = Palatability
    ),
    position = position_dodge(width = 0.75),
    linewidth = 1.3,
    width = 0.1
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
  ) + labs(subtitle = '') + ylab("Relative biomass (g)") +
  theme(
    panel.grid.major = element_line(colour = "gray86", size = 0.1),
    panel.background = element_rect(fill = "white")
  )+
  theme(axis.ticks = element_blank())

fig_rel_biomass

ggsave('Fig_2.jpg',
       width = 10,
       height = 6,
       dpi = 300)
