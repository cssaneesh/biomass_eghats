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
transect_dat <- raw_dat %>%
  mutate(
    Treatment = as.factor(Treatment),
    Life_form = as.factor(Life_form),
    Functional_groups = as.factor(Functional_groups)
  )

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

# Relative_biomass
relative_weight <-
  transect_calc %>% select(Transect, Treatment, Palatability, relative_biomass) %>%
  group_by(Treatment, Transect, Palatability) %>%
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

# Analysis----

# ghats.rel_biomass <-
#   brm(
#     relative_biomass ~   Treatment * Palatability  +
#       (Treatment * Palatability  | Transect) ,
#     family = student(),
#     data = relative_weight,
#     iter = 2000,
#     warmup = 1000,
#     cores = 4,
#     chains = 4,
#     backend = 'rstan'
#   )
#
# save(ghats.rel_biomass, file = "ghats.rel_biomass.Rdata")
load("ghats.rel_biomass.Rdata")

color_scheme_set("darkgray")

density_plot <- pp_check(ghats.rel_biomass) +
  xlab("Relative biomass") + ylab("Density") +
  labs(subtitle = "a)") +
  theme_classic() +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        legend.position = "none")# predicted vs. observed values

density_plot

summary(ghats.rel_biomass)

ghats_rel_biomass <-
  conditional_effects(
    ghats.rel_biomass,
    effects = 'Treatment:Palatability',
    re_formula = NA,
    method = 'fitted'
  )  # conditional effects

ghats_rel_biomass_df <-
  as.data.frame(ghats_rel_biomass$`Treatment:Palatability`)# to make a df to report statistics results

# View(ghats_rel_biomass_df)

# Table_1 ----
# Biomass
table_1 <- ghats_rel_biomass_df %>%
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

table_1
table_1 %>% gtsave('Table_1.png', expand = 5) # expand to set white space

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
    alpha = 0.5,
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
    position = position_dodge(width = 0.75)
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
    size = 1,
    width = 0
  ) + labs(x = '', y = '') +
  scale_color_manual(values =  c("#f8bb49", '#a45200', "#4A2300"))  +
  # scale_colour_grey()+
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
    legend.position = "bottom"
  ) + labs(subtitle = '') + ylab("Relative biomass") +
  theme(
    panel.grid.major = element_line(colour = "gray86", size = 0.1),
    panel.background = element_rect(fill = "white")
  )

fig_rel_biomass + plot_annotation(title = "Relative biomass",
                                  theme = theme(plot.title = element_text(size = 14, hjust = 0.5)))

# add treatment icons to x axis
treats <- axis_canvas(fig_rel_biomass, axis = 'x') +
  cowplot::draw_image('CPFP.png', x = 0.5, scale = 0.5) +
  cowplot::draw_image('CPFA.png', x = 1.5, scale = 0.5) +
  cowplot::draw_image('CAFA.png', x = 2.5, scale = 0.5)

Fig_2 <-
  ggdraw(insert_xaxis_grob(fig_rel_biomass, treats, position = "center"))

Fig_2 <- Fig_2 + plot_annotation(title = "Relative biomass",
                        theme = theme(plot.title = element_text(size = 14, hjust = 0.5)))
Fig_2

# Save image (Biomass/Fig_2)
ggsave('Fig_2.jpg',
       width = 10,
       height = 6,
       dpi = 300)


# table test----


test <- ghats_rel_biomass_df %>%
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

test

test %>% gtsave('test.png', expand = 5) # expand to set white spacegtsave("test.jpg")

  