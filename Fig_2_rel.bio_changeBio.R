source('1_DataPackages.R')

# relative biomass-----
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

fig_s2 <- pp_check(ghats.rel_biomass) +
  xlab("Relative biomass") + ylab("Density") +
  labs(title = "Site-level",
    #subtitle = "a)"
    ) +
  theme_classic() + xlim(-20,150) +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        legend.position = "none")# predicted vs. observed values

fig_s2

# ggsave('Figure S2.jpg',
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
  ) + labs(subtitle = '') + ylab("Relative biomass") +
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

# biomass change----
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

# # Analysis----
# ghats.biomass.change <-
#   brm(
#     Biomass ~   Treatment * Palatability+
#       ( 1 | Village/Sci_name ) ,
#     family = student(),
#     data = bio.change.treatment,
#     iter = 10000,
#     warmup = 1000,
#     cores = 4,
#     chains = 4,
#     control = list(adapt_delta = 0.99, max_treedepth = 12 )
#   )

# save(ghats.biomass.change, file = "ghats.biomass.change.Rdata")

load("ghats.biomass.change.Rdata")

summary(ghats.biomass.change)

# fig_sx <- pp_check(ghats.biomass.change)+
#   xlab("Relative biomass change") + ylab("Density") +
#   labs(title = "Site-level",
#        #subtitle = "a)"
#   ) +
#   theme_classic() + xlim(-40,150) +
#   theme(plot.title = element_text(size = 18, hjust = 0.5),
#         legend.position = "none")# predicted vs. observed values
# 
# fig_sx
# 
# ggsave('Figure Sx.jpg',
#        width = 10,
#        height = 6,
#        dpi = 300)

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
    panel.grid.major = element_line(colour = "gray86", linewidth = 0.1),
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

# ggplotly(fig_change.biomass)

# Relative weight----
# what is the summed biomass per Site with cymbopogon?


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
#     iter = 10000,
#     warmup = 1000,
#     cores = 4,
#     chains = 4,
#     control = list(adapt_delta = 0.99, max_treedepth = 12 )
#   )

# save(rel.ghats.biomass.change, file = "rel.ghats.biomass.change.Rdata")

load("rel.ghats.biomass.change.Rdata")

summary(rel.ghats.biomass.change)



# fig_s3 <- pp_check(rel.ghats.biomass.change)+
#   xlab("Relative biomass change") + ylab("Density") +
#   labs(title = "Site-level",
#        #subtitle = "a)"
#   ) +
#   theme_classic() + xlim(-10,10) +
#   theme(plot.title = element_text(size = 18, hjust = 0.5),
#         legend.position = "none")# predicted vs. observed values
# 
# fig_s3
# 
# ggsave('Figure S3.jpg',
#        width = 10,
#        height = 6,
#        dpi = 300)


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
  ) + labs(subtitle = '') + ylab("Relative biomass change") +
  theme(
    panel.grid.major = element_line(colour = "gray86", linewidth  = 0.1),
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
    linewidth = 0.3
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
    linewidth = 0.3
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
    linewidth = 0.3
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
    linewidth = 0.3
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
    linewidth = 0.3
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
    linewidth = 0.3
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
    linewidth = 0.3
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
    linewidth = 0.3
  )+ # cafa no -7.72
  theme(legend.position = '')

# fig_rel.change.biomass
# save(fig_rel.change.biomass, file= 'fig_rel.change.biomass.Rdata')

# fig 2& 3-----
load('fig_rel.change.biomass.Rdata')
load('fig_rel_biomass.Rdata')
load('legend_rel_biomass.Rdata')

# library("gridExtra")
# Create user-defined function, which extracts legends from ggplots
extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}

# Apply user-defined function to extract legend
shared_legend <- extract_legend(legend_rel_biomass)

# Draw plots with shared legend
fig_2 <- grid.arrange(arrangeGrob(fig_rel_biomass, fig_rel.change.biomass, ncol = 2),
                       shared_legend, nrow = 2, heights = c(10, 1))

# ggsave('Fig_x_biochng.jpg',
#        width = 10,
#        height = 6,
#        dpi = 300)

# ggsave('Fig_2_biomassAndchng.jpg', fig_2, 
#        width = 10,
#        height = 6,
#        dpi = 300)



# plotly::ggplotly(fig_rel.change.biomass)

rel.ghats_biochange_df <- as.data.frame(rel.ghats_change_biomass$`Treatment:Palatability`)
# View(rel.ghats_biochange_df)

TableS4 <- rel.ghats_biochange_df %>% 
  select(Treatment, Palatability, estimate__, lower__, upper__) %>% 
  rename(Estimate= estimate__, 
         Lower= lower__,
         Upper = upper__) %>% 
  mutate(Estimate = round(Estimate, 2)) %>% 
  mutate(Lpwer = round(Lower, 2)) %>%
  mutate(Upper = round(Upper, 2)) %>%
  gt() %>% tab_options(column_labels.font.size = 11,
                       table.font.size = 10,
                       column_labels.font.weight = "bold") %>% 
  tab_header(subtitle = '', 'Relative biomass change') %>% 
  opt_table_font(default_fonts()) %>%  # Fonts: Roboto Mono,IBM Plex Mono, Red Hat Mono
  opt_table_outline(style = "solid", width = px(2))

TableS4
TableS4 %>% gtsave('Table_S4_biochange.png', expand = 5)

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

