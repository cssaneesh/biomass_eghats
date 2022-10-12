# Install packages----
# install.packages("remotes")
# remotes::install_github("KaiHsiangHu/iNEXT.3D")
# vignette: https://cran.r-project.org/web/packages/iNEXT/vignettes/Introduction.pdf

# library ----
library(tidyverse)
library(patchwork)
library(iNEXT.3D)
library(viridis)
# Data----
transect_dat <- read.csv(
  "biomass_data.csv",
  header = T,
  fill = TRUE,
  sep = ",",
  na.strings = c("", " ", "NA", "NA ", "na", "NULL")
)

# Data wrangling----
transect_dat <- transect_dat %>%
  mutate(
    Treatment = as.factor(Treatment),
    Life_form = as.factor(Life_form),
    Functional_groups = as.factor(Functional_groups)
  )

# Analysis----
transect_prep_iNext <- transect_dat %>%
  mutate(species = Sci_name) %>%
  arrange(Transect, Treatment) %>%
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

transect.list <- transect_prep_iNext %>%
  split(.$Treatment)

transect.matrix.list <- purrr::map(
  transect.list,
  ~ .x %>%
    select(species, Transect, pres) %>%
    distinct() %>%
    spread(key = Transect, value = pres) %>%
    replace(is.na(.), 0) %>%
    column_to_rownames(var = "species")
)

#  Taxonomic diversity
TD_treat_out <-
  iNEXT3D(
    data = transect.matrix.list,
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
transect.TD.df <- TD_treat_out %>%
  purrr::pluck("iNextEst", "size_based")

transect_info <- transect_prep_iNext %>%
  distinct() %>% mutate(Assemblage = as.character(Treatment))

transect.hill.TD <- transect.TD.df %>% left_join(transect_info) %>%
  mutate(Order.q  = case_when(Order.q  == "0" ~ "q = 0", # q=0 species richness
                              Order.q == "1" ~ "q = 1", # q=1 Shannon diversity
                              Order.q == "2" ~ "q = 2")) %>% # q=2 Simpson diversity or evenness
  filter(!Order.q == "q = 1")

df.point <-
  transect.hill.TD[which(transect.hill.TD$Method == "Observed"), ]
df.line <-
  transect.hill.TD[which(transect.hill.TD$Method != "Observed"), ]
df.line$Method <- factor(df.line$Method,
                         c("Rarefaction", "Extrapolation"))
# Rarefaction or Interpolation?

# Plot----
transect.TD.fig_op1 <- ggplot(transect.hill.TD ,
                              aes(x = nt, y = qD,   color = Treatment)) +
  facet_wrap( ~ Order.q) +
  geom_point(aes(),
             shape = 1,
             size = 2,
             data = df.point) +
  geom_line(aes(linetype = Method), lwd = 0.75, data = df.line) +
  labs(
    x = "Number of sampling units",
    y = "Taxonomic diversity",
    title =
      "Sample - based diversity accumulation ",
    caption = "Control= Cymbopogon present and fire present
CPFA= Cymbopogon present and fire absent
CAFA= Cymbopogon absent and fire absent"
  ) +
  scale_color_manual(values = c(
    "Control" = "#BB9689",
    "CPFA" = "#836656",
    "CAFA" = "#6C3859"
  )) +
  # scale_color_manual(values =  c("#A6BAd7", '#836656',"#6C3859"))  +
  # scale_colour_grey()+
  theme_bw(base_size = 12) +
  theme(legend.text = element_text(size = 8)) +
  guides(col = guide_legend(ncol = 15)) +
  theme(plot.title = element_text(size = 14, hjust = 0.5)) +
  theme(
    panel.grid.major = element_line(colour = "gray86", size = 0.1),
    panel.background = element_rect(fill = "white")
  )



(
  transect.TD.fig_op1 +
    #   labs(caption = "Control= Cymbopogon present and fire present
    # CPFA= Cymbopogon present and fire absent
    # CAFA= Cymbopogon absent and fire absent")+
    theme(plot.caption = element_text(
      size = 8, face = "italic",
      hjust = 0.0
    )) +
    theme(legend.position = 'bottom',) +
    theme(
      legend.position = c(0.95, .75),
      legend.direction = "vertical"
    ) +
    theme(legend.background = element_rect(fill = NA))
)

ggsave('fig_annex_acu_curve.jpg',
       width = 10,
       height = 6,
       dpi = 300)  

transect_prep_iNext  %>%
  group_by(Treatment) %>%
  distinct(species) %>%
  summarise(species = n()) # note the number of Family and run the following code

TD_treat_out$DataInfo # S.obs= Observed species richness will be the same, this is the q=0 in the plot

TD_treat_out$iNextEst

