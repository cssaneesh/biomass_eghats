# run the following line to load libraries and add data.
source('1_DataPackages.R')

# Exploratory plots----
# Absolute weight
fig_e1_ab <-
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
  labs(title = "Absolute biomass") +
  scale_fill_manual(values = c("#f8bb49", '#a45200', "#4A2300")) +
  theme(plot.caption = element_text(size = 8, face = "italic",
                                    hjust = 0))
fig_e1_ab

# Analysis----

# ghats.ab_biomass <-
#   brm(
#     Weight ~   Treatment * Palatability  +
#       ( 1 | Site ) ,
#     family = student(),
#     data = absolute_weight,
#     iter = 2000,
#     warmup = 1000,
#     cores = 4,
#     chains = 4
#   )
# save(ghats.ab_biomass, file = "ghats.ab_biomass.Rdata")

load("ghats.ab_biomass.Rdata")

color_scheme_set("darkgray")

fig_s1ab <- pp_check(ghats.ab_biomass) +
  xlab("Functional group absolute biomass") + ylab("Density") +
  labs(title = "Site-level",
       subtitle = "a)") +
  theme_classic() + xlim(-20,150) +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        legend.position = "none")# predicted vs. observed values

fig_s1ab

summary(ghats.ab_biomass)

ghats_ab_biomass <-
  conditional_effects(
    ghats.ab_biomass,
    effects = 'Treatment:Palatability',
    re_formula = NA,
    method = 'fitted'
  )  # conditional effects

ghats_ab_biomass_df <-
  as.data.frame(ghats_ab_biomass$`Treatment:Palatability`)# to make a df to report statistics results

# View(ghats_ab_biomass_df)
# Table_1 ----
# Biomass
TableS2 <- ghats_ab_biomass_df %>%
  select(Treatment, Palatability, estimate__, lower__, upper__) %>%
  rename(Estimate = estimate__,
         Lower = lower__,
         Upper = upper__) %>%
  mutate_if(is.numeric, round, 2) %>% 
  gt() %>% 
  tab_options(column_labels.font.size = 11,
              table.font.size = 10,
              column_labels.font.weight = "bold") %>% 
  tab_header(subtitle = '', 'Absolute biomass') %>% 
  opt_table_font(default_fonts()) %>%  # Fonts: Roboto Mono,IBM Plex Mono, Red Hat Mono
  opt_table_outline(style = "solid", width = px(2))

TableS2
TableS2 %>% gtsave('Table_S2_ab_biomass.png', expand = 5) # expand to set white space

# Plot----
fig_ab_biomass <- ggplot() +
  geom_point(
    data = absolute_weight,
    aes(
      x = Treatment,
      y = Weight,
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
    data = ghats_ab_biomass$`Treatment:Palatability`,
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
    data = ghats_ab_biomass$`Treatment:Palatability`,
    aes(
      x = Treatment,
      ymin = lower__,
      ymax = upper__,
      group = Palatability,
      colour = Palatability
    ),
    position = position_dodge(width = 0.75),
    linewidth = 1,
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
  ) + labs(subtitle = '') + ylab("Absolute biomass (g)") +
  theme(
    panel.grid.major = element_line(colour = "gray86", linewidth = 0.1),
    panel.background = element_rect(fill = "white")
  )

fig_ab_biomass + plot_annotation(title = "Absolute biomass",
                                 theme = theme(plot.title = element_text(size = 14, hjust = 0.5)))

# add treatment icons to x axis
treats <- axis_canvas(fig_ab_biomass, axis = 'x') +
  cowplot::draw_image('CPFP.png', x = 0.5, scale = 0.5) +
  cowplot::draw_image('CPFA.png', x = 1.5, scale = 0.5) +
  cowplot::draw_image('CAFA.png', x = 2.5, scale = 0.5)

Fig_2_ab <-
  ggdraw(insert_xaxis_grob(fig_ab_biomass, treats, position = "center"))

Fig_2_ab <- Fig_2_ab + plot_annotation(title = "Absolute biomass",
                                 theme = theme(plot.title = element_text(size = 14, hjust = 0.5)))
Fig_2_ab

# Save image (Biomass/Fig_2)
ggsave('Fig_2_ab.jpg',
       width = 10,
       height = 6,
       dpi = 300)