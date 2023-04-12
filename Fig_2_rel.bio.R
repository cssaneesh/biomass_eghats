source('1_DataPackages.R')

# Exploratory plots----

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

fig_e2

head(relative_weight)

# Analysis----
# ghats.rel_biomass <-
#   brm(
#     relative_biomass ~   Treatment * Palatability  +
#       ( 1 | Site ),
#     #family = gaussian(),
#       family = student(),
#     # family = lognormal(),
#     data = relative_weight,
#     iter = 10000,
#     warmup = 1000,
#     cores = 4,
#     chains = 4,
#     control = list(adapt_delta = 0.99 , max_treedepth = 12 )
#   )

# save(ghats.rel_biomass, file = "ghats.rel_biomass.Rdata")

load("ghats.rel_biomass.Rdata")

color_scheme_set("darkgray")

fig_s3 <- pp_check(ghats.rel_biomass) +
  xlab("Relative biomass") + ylab("Density") +
  labs(title = "Site-level",
    #subtitle = "a)"
    ) +
  theme_classic() + xlim(-20,150) +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        legend.position = "none")# predicted vs. observed values

fig_s3

# ggsave('sup_Fig_3.jpg',
#        width = 10,
#        height = 6,
#        dpi = 300)


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
# write.csv(ghats_rel_biomass_df, 'ghats_rel_biomass_df.csv')

# Table_1 ----
# Biomass
TableS3 <- ghats_rel_biomass_df %>%
  select(Treatment, Palatability, estimate__, lower__, upper__) %>%
  rename(Estimate = estimate__,
         Lower = lower__,
         Upper = upper__) %>%
  mutate_if(is.numeric, round, 2) %>% 
  gt() %>% 
  tab_options(column_labels.font.size = 11,
              table.font.size = 10,
              column_labels.font.weight = "bold") %>% 
  tab_header(subtitle = '', 'Relative biomass') %>% 
  opt_table_font(default_fonts()) %>%  # Fonts: Roboto Mono,IBM Plex Mono, Red Hat Mono
  opt_table_outline(style = "solid", width = px(2))

TableS3
TableS3 %>% gtsave('Table_S3_rel_bio.png', expand = 5) # expand to set white space

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
    alpha= 0.6,
    position = position_dodge(width = 0.75,)
  ) +
  geom_errorbar(
    data = ghats_rel_biomass$`Treatment:Palatability`,
    aes(
      x = Treatment,
      ymin = lower__,
      ymax = upper__,
      group = Palatability,
      # colour = Palatability
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
  ) + labs(subtitle = '') + ylab("Relative biomass (g)") +
  theme(
    panel.grid.major = element_line(colour = "gray86", linewidth = 0.1),
    panel.background = element_rect(fill = "white")
  )+
  theme(axis.ticks = element_blank())+
  theme(legend.position = ' ')

# save(fig_rel_biomass, file= 'fig_rel_biomass.Rdata')
# 
# legend_rel_biomass <-  fig_rel_biomass +
#   theme(legend.position = 'bottom')
# 
# save(legend_rel_biomass, file= 'legend_rel_biomass.Rdata')

# ggsave('Fig_2.jpg',
#        width = 10,
#        height = 6,
#        dpi = 300)
