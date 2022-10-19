# Install packages----
# install.packages("remotes")
# remotes::install_github("KaiHsiangHu/iNEXT.3D")
# vignette: https://cran.r-project.org/web/packages/iNEXT/vignettes/Introduction.pdf

# library ----
library(tidyverse)
library(patchwork)
library(iNEXT.3D)
library(viridis)
library(gt)
# Data----
Site_dat <- read.csv(
  "biomass_data.csv",
  header = T,
  fill = TRUE,
  sep = ",",
  na.strings = c("", " ", "NA", "NA ", "na", "NULL")
)

# Data wrangling----
Site_dat <- Site_dat %>%
  mutate(
    Treatment = as.factor(Treatment),
    Life_form = as.factor(Life_form),
    Functional_groups = as.factor(Functional_groups)
  )

# Analysis----
Site_prep_iNext <- Site_dat %>%
  mutate(species = Sci_name) %>%
  arrange(Site, Treatment) %>%
  select(-c(Palatability, Sci_name)) %>%
  mutate(
    Treatment = case_when(Treatment == "ab" ~ "Control", # Cymbopogon present fire present
                          Treatment == "bgpnf" ~ "CPFA", # Cymbopogon present fire absent
                          Treatment == "bgrnf" ~ "CAFA" # Cymbopogon absent fire absent
    ),
  ) %>%
  mutate(Cymb_status = case_when(
    !species == "Cymbopogon sp." ~ "Other Sp.",
    TRUE ~ as.character(species)
  ),) %>%
  mutate(pres = as.numeric(1)) %>% 
  mutate(Treatment = factor(Treatment)) %>% # to order treatments in the plot
  mutate(Treatment = fct_relevel(Treatment, c("Control","CPFA","CAFA")))

Site.list <- Site_prep_iNext %>%
  split(.$Treatment)

Site.matrix.list <- purrr::map(
  Site.list,
  ~ .x %>%
    select(species, Site, pres) %>%
    distinct() %>%
    spread(key = Site, value = pres) %>%
    replace(is.na(.), 0) %>%
    column_to_rownames(var = "species")
)

#  Taxonomic diversity
TD_treat_out <-
  iNEXT3D(
    data = Site.matrix.list,
    diversity = 'TD',
    q = c(0, 1, 2),
    datatype = 'incidence_raw',
    size = c(1:60),
    nboot = 0
  )

TD_treat_out$DataInfo
# Assemblage = the treatment or groups
# T = Reference sample size
# U = Total number of incidents
# S.obs = Observed species richness
# SC = Sample coverage

TD_treat_out$AsyEst # to see the asymptotic diversity estimates

# Make df for ploting----
Site.TD.df <- TD_treat_out %>%
  purrr::pluck("iNextEst", "size_based")

Site_info <- Site_prep_iNext %>%
  distinct() %>% mutate(Assemblage = as.character(Treatment))

Site.hill.TD <- Site.TD.df %>% left_join(Site_info) %>%
  mutate(Order.q  = case_when(Order.q  == "0" ~ "q = 0", # q=0 species richness
                              Order.q == "1" ~ "q = 1", # q=1 Shannon diversity
                              Order.q == "2" ~ "q = 2")) %>% # q=2 Simpson diversity or evenness
  filter(!Order.q == "q = 1")

df.point <-
  Site.hill.TD[which(Site.hill.TD$Method == "Observed"), ]
df.line <-
  Site.hill.TD[which(Site.hill.TD$Method != "Observed"), ]
df.line$Method <- factor(df.line$Method,
                         c("Rarefaction", "Extrapolation"))
# Rarefaction or Interpolation?

# Plot----
treatment_colors <- c("Control" = "#3b5d4d", "CPFA" = "#c5af99","CAFA" = "#ffd365")

fig_5 <- ggplot(Site.hill.TD ,
                              aes(x = nt, y = qD,   color = Treatment)) +
  facet_wrap( ~ Order.q) +
  geom_point(aes(),
             shape = 1,
             size = 3,
             data = df.point) +
  geom_line(aes(linetype = Method), lwd = 0.75, data = df.line) +
  labs(
    x = "Number of sites",
    y = "Taxonomic diversity",) +
  scale_color_manual(values = treatment_colors) +
  theme_bw(base_size = 12) +
  theme(legend.text = element_text(size = 8)) +
  guides(col = guide_legend(ncol = 15)) +
  theme(plot.title = element_text(size = 14, hjust = 0.5)) +
  theme(
    panel.grid.major = element_line(colour = "gray86", size = 0.1),
    panel.background = element_rect(fill = "white")
  )

fig_5

(
  fig_5 +
    theme(plot.caption = element_text(
      size = 8, face = "italic",
      hjust = 0.0
    )) +theme(
      legend.position = c(1.3, .95),
     legend.justification = c('right', 'top')
    ) +
    theme(legend.background = element_rect(fill = NA))
)

ggsave('fig_5_acu_curve.jpg',
       width = 10,
       height = 6,
       dpi = 300)  

TD_treat_out$DataInfo # S.obs= Observed species richness will be the same, this is the q=0 in the plot
TD_treat_out$AsyEst
iNextEststimates <- as.data.frame(TD_treat_out$AsyEst)

# iNext table
table_4_iNext <- iNextEststimates %>% select(Assemblage, Diversity, Observed, Estimator) %>%
  filter(Diversity!="Shannon diversity") %>% 
  mutate_if(is.numeric, round, 2) %>% 
  gt()%>% 
  tab_options(column_labels.font.size = 11,
              table.font.size = 10,
              column_labels.font.weight = "bold") %>% 
  opt_table_font(default_fonts()) %>%  # Fonts: Roboto Mono,IBM Plex Mono, Red Hat Mono
  opt_table_outline(style = "solid", width = px(2))

table_4_iNext %>% gtsave('Table_4 (iNext_AsyEst).png', expand = 5) # expand to set white space

