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

alpha_dat <-  transect_dat %>% 
  group_by(Transect, Site, Treatment, Sci_name) %>%
  summarise(weight = sum(Weight)) %>% arrange(Transect, Sci_name)

plot_dat <- alpha_dat %>% 
  group_by(Transect, Site, Treatment) %>%
  summarise(plot_weight = sum(weight))

alpha_dat_prep <- alpha_dat %>% 
  left_join(plot_dat) %>%
  mutate(rel_weight = (weight/plot_weight))

transect_prep <- transect_dat %>%
  arrange(Transect, Treatment) %>%
  mutate(
    Treatment = case_when(
      Treatment == "ab" ~ "Control", # Cymbopogon present fire present
      Treatment == "bgpnf" ~ "CPFA", # Cymbopogon present fire absent
      Treatment == "bgrnf" ~ "CAFA" # Cymbopogon absent fire absent
    )
  )

# what is the summed biomass per transect without cymbopogon?
transect_sum_no_cym <- transect_prep %>%
  filter(!Sci_name == "Cymbopogon sp.") %>%
  group_by(Transect, Treatment) %>%
  summarise(transect_biomass_no_cym = sum(Weight)) %>%
  ungroup()

# what is the summed biomass per transect with cymbopogon?
transect_sum_w_cym <-
  transect_prep %>% group_by(Transect, Treatment) %>%
  summarise(transect_biomass = sum(Weight)) %>%
  ungroup()

transect_sum <-
  transect_sum_w_cym %>% left_join(transect_sum_no_cym)

# Transect_calc
transect_calc <- transect_prep %>% left_join(transect_sum) %>%
  mutate(
    relative_biomass = (Weight / transect_biomass) ,
    relative_biomass_p = ((Weight / transect_biomass) * 100),
    relative_biomass_nc = (Weight / transect_biomass_no_cym) ,
    relative_biomass_nc_p = ((Weight / transect_biomass) * 100)
  )

View(transect_calc)

alpha_div <-
  transect_calc %>% group_by(Transect, Treatment) %>%
  dplyr::summarise(
    alpha_rich = n_distinct(Sci_name),
    alpha_ENSPIE = vegan::diversity(relative_biomass,
                                    index = 'invsimpson')) %>% 
  mutate(Treatment = factor(Treatment)) %>% 
  mutate(Treatment = fct_relevel(Treatment, c("Control","CPFA","CAFA")))  %>%
  ungroup()


# Analysis-----

# ghats.alpha_ENSPIE----

ghats.alpha_ENSPIE <-
  brm(
    alpha_ENSPIE ~   Treatment + ( Treatment  | Transect) ,
    family = 'lognormal',
    data = alpha_div,
    iter = 6000,
    warmup = 1000,
    cores = 4,
    chains = 4,
    control = list(adapt_delta = 0.85) )
 
 save(ghats.alpha_ENSPIE, file = 'ghats.alpha_ENSPIE.Rdata')

load('ghats.alpha_ENSPIE.Rdata')


summary(ghats.alpha_ENSPIE) # summary of alpha richness model

color_scheme_set("darkgray")
# caterpillars/chains
plot(ghats.alpha_ENSPIE)
# you want these 'caterpillars to be 'hairy' (very evenly squiggly)

# check model residuals
head(alpha_div)
ma <- residuals(ghats.alpha_ENSPIE)
ma <- as.data.frame(ma)
ar.plot <- cbind(alpha_div, ma$Estimate)

#make sure they are factors
ar.plot$Treatment <- as.factor(ar.plot$Treatment )
ar.plot$Transect <- as.factor(ar.plot$Transect )

#plot residuals
par(mfrow=c(1,2))
with(ar.plot, plot(Treatment, ma$Estimate))
with(ar.plot, plot(Transect, ma$Estimate))
# you want these to be centrered on zero

fig_s4a <- pp_check(ghats.alpha_ENSPIE) +
    xlab( expression(paste(ENS[PIE])) ) + ylab("Density") + 
  ggtitle((expression(paste(italic(alpha), '-scale', sep = ''))))+
  labs(subtitle = "a)") +
  theme_classic() + xlim(-5,20) +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        legend.position = "none")# predicted vs. observed values

fig_s4a

ghats_alpha_ENSPIE <-
  conditional_effects(
    ghats.alpha_ENSPIE,
    effects = 'Treatment',
    re_formula = NA,
    method = 'fitted'
  )

# beta

# the solution, bootstrap resampling: prepare the data for bootstrap resampling of sites
gamma_dat <- alpha_dat_prep %>% 
  # collate relative weight of each species at each location (these are alpha-scale samples)
  group_by(Treatment, Transect) %>% 
  nest(c(Sci_name, rel_weight, weight, plot_weight)) %>% 
  ungroup()

# for n_samps, get 10 sites (alpha samples)

# n_sites = 10
# n_samps <- 200
# gamma_metrics <- tibble()
# for(i in 1:n_samps){
#   print(i)
#   # get these n_sites rows and calculate alpha S
#   alpha_sub_samp <- gamma_dat %>%
#     # from each group
#     group_by(Treatment) %>%
#     # get 10 rows
#     sample_n(n_sites, replace = F) %>%
#     # unnest
#     unnest() %>%
#     # calculate PIE, S for each site
#     group_by(Treatment, Transect) %>%
#     mutate(alphaS = n_distinct(Sci_name),
#            alpha_Spie = vegan::diversity(rel_weight, index = 'invsimpson')) %>%
#     ungroup() %>%
#     # get the minimum N and mean S for each treatment
#     group_by(Treatment) %>%
#     mutate(mean_alpha_S = mean(alphaS),
#            mean_alpha_Spie = mean(alpha_Spie)) %>%
#     ungroup()
# 
#   # aggregate same sub sample for gamma calculations
#   sub_samp <- alpha_sub_samp %>%
#     # aggregate data to gamma scale
#     group_by(Treatment, Sci_name) %>%
#     summarise(sp_trt_weight = sum(weight)) %>%
#     ungroup() %>%
#     # get minimum N for Sn
#     group_by(Treatment) %>%
#     mutate(trt_weight = sum(sp_trt_weight),
#            gamma_rel_weight = (sp_trt_weight/trt_weight)) %>%
#     ungroup() %>%
#     mutate(minrel = min(gamma_rel_weight))
# 
# 
#   # calculate the metrics we want
#   gamma_metrics <- gamma_metrics %>%
#     bind_rows(sub_samp %>%
#                 group_by(Treatment) %>%
#                 summarise(S = n_distinct(Sci_name),
#                           ENSPIE = vegan::diversity(gamma_rel_weight, index = 'invsimpson'))  %>%
#                 # add counter for sample based rarefaction
#                 left_join( alpha_sub_samp %>%
#                              select(Treatment, mean_alpha_S, mean_alpha_Spie) %>%
#                              distinct() %>%
#                              group_by(Treatment) %>%
#                              mutate(
#                                alpha_S = mean_alpha_S,
#                                alpha_Spie = mean_alpha_Spie,
#                                resample = i) ) )
# }
# 
# save(gamma_metrics, file= 'gamma_metrics.Rdata')

load('gamma_metrics.Rdata')
gamma_boot_results <- gamma_metrics %>% # calculate beta-diversities (beta = gamma/alpha) 
  mutate(beta_S = S/alpha_S,
         beta_S_PIE = ENSPIE/alpha_Spie ) %>% 
  group_by(Treatment) %>% 
  summarise(S_mean = mean(S),
            S_median = median(S),
            S_Q95 = quantile(S, probs = 0.95, names = F),
            S_Q5 = quantile(S, probs = 0.05, names = F),
            ENSPIE_mean = mean(ENSPIE),
            ENSPIE_median = median(ENSPIE),
            ENSPIE_Q95 = quantile(ENSPIE, probs = 0.95, names = F),
            ENSPIE_Q5 = quantile(ENSPIE, probs = 0.05, names = F),
            beta_S_mean = mean(beta_S),
            beta_S_median = median(beta_S),
            beta_S_Q95 = quantile(beta_S, probs = 0.95, names = F),
            beta_S_Q5 = quantile(beta_S, probs = 0.05, names = F),
            beta_S_PIE_mean = mean(beta_S_PIE),
            beta_S_PIE_median = median(beta_S_PIE),
            beta_S_PIE_Q95 = quantile(beta_S_PIE, probs = 0.95, names = F),
            beta_S_PIE_Q5 = quantile(beta_S_PIE, probs = 0.05, names = F))  %>%
  mutate( Treatment = case_when( 
    Treatment == "ab" ~ "Control", # Cymbopogon present fire present
    Treatment == "bgpnf" ~ "CPFA", # Cymbopogon present fire absent
    Treatment == "bgrnf" ~ "CAFA" # Cymbopogon absent fire absent
  )) %>% 
  mutate(Treatment = factor(Treatment)) %>% # to order treatments in the plot
  mutate(Treatment = fct_relevel(Treatment, c("Control","CPFA","CAFA"))) 

# Tables----
# table alpha_ENSPIE
ghats_alpha_ENSPIE_df <- as.data.frame(ghats_alpha_ENSPIE$Treatment)
View(ghats_alpha_ENSPIE_df)

table_3_ens_alpha <- ghats_alpha_ENSPIE_df %>% select(Treatment, estimate__, lower__, upper__) %>% 
  rename(Estimate=estimate__, Lower=lower__, Upper= upper__) %>% 
  mutate_if(is.numeric, round, 2) %>% 
  gt()%>% 
  tab_options(column_labels.font.size = 11,
              table.font.size = 10,
              column_labels.font.weight = "bold") %>% 
  opt_table_font(default_fonts()) %>%  # Fonts: Roboto Mono,IBM Plex Mono, Red Hat Mono
  opt_table_outline(style = "solid", width = px(2))

table_3_ens_alpha %>% gtsave('Table_3 (ens_alpha).png', expand = 5) # expand to set white space

# Table beta_ENSPIE
table_3_ens_beta <- gamma_boot_results %>% 
  select(Treatment, beta_S_PIE_mean, beta_S_PIE_Q5, beta_S_PIE_Q95) %>% 
  rename(Estimate= beta_S_PIE_mean, Lower= beta_S_PIE_Q5, Upper= beta_S_PIE_Q95) %>% 
  mutate_if(is.numeric, round, 2) %>% 
  gt()%>% 
  tab_options(column_labels.font.size = 11,
              table.font.size = 10,
              column_labels.font.weight = "bold") %>% 
  opt_table_font(default_fonts()) %>%  # Fonts: Roboto Mono,IBM Plex Mono, Red Hat Mono
  opt_table_outline(style = "solid", width = px(2))

table_3_ens_beta %>% gtsave('Table_3 (ens_beta).png', expand = 5) # expand to set white space

# Table gamma_ENSPIE
table_3_ens_gamma <- gamma_boot_results %>% 
  select(Treatment, ENSPIE_mean, ENSPIE_Q5, ENSPIE_Q95) %>% 
  rename(Estimate= ENSPIE_mean, Lower= ENSPIE_Q5, Upper= ENSPIE_Q95) %>% 
  mutate_if(is.numeric, round, 2) %>% 
  gt()%>% 
  tab_options(column_labels.font.size = 11,
              table.font.size = 10,
              column_labels.font.weight = "bold") %>% 
  opt_table_font(default_fonts()) %>%  # Fonts: Roboto Mono,IBM Plex Mono, Red Hat Mono
  opt_table_outline(style = "solid", width = px(2))

table_3_ens_gamma %>% gtsave('Table_3 (ens_gamma).png', expand = 5) # expand to set white space


# Plot----
fig_alpha_ENSPIE <- ggplot() +
  geom_point(
    data = alpha_div,
    aes(x = Treatment, y = alpha_ENSPIE, colour = 	'#A6BAd7'),
    size = 1,
    alpha = 0.7,
    position = position_jitter(width = 0.05, height = 0.45)
  ) +
  geom_point(
    data = ghats_alpha_ENSPIE$Treatment,
    aes(x = Treatment, y = estimate__, colour = Treatment),
    size = 3
  ) +
  geom_errorbar(
    data = ghats_alpha_ENSPIE$Treatment,
    aes(
      x = Treatment,
      ymin = lower__,
      ymax = upper__,
      colour = Treatment
    ),
    size = 1,
    width = 0
  ) + labs(x = '', y = '') +
  scale_color_manual(values = c("#A6BAd7",
                                "Control" = "#BB9689",
                                "CPFA" = "#836656",
                                "CAFA" = "#6C3859"))+
  labs(subtitle = 'a)')+
  ylab(expression(paste(italic(alpha), '-',
                        ENS[PIE])))+
  theme_bw(base_size = 12) +
  theme(legend.position = 'none', 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.tag.position = c(0.3, 0.8))+
  theme(panel.grid.major = element_line(colour = "gray86", size = 0.1),
        panel.background = element_rect(fill = "white"))

fig_4a <- fig_alpha_ENSPIE


beta_S_PIE_all <- ggplot() +
  geom_point(data = gamma_boot_results,
             aes(x = Treatment, y = beta_S_PIE_mean, colour = Treatment),
             size = 4) +
  geom_errorbar(data = gamma_boot_results,
                aes(x = Treatment, ymin = beta_S_PIE_Q5, ymax = beta_S_PIE_Q95, 
                    colour = Treatment),
                size = 1.3,
                width = 0.1) +
  scale_color_manual(values = c("Control" = "#BB9689",
                                "CPFA" = "#836656",
                                "CAFA" = "#6C3859"))+
labs(subtitle = 'b)', x = '',
       y = expression(paste(italic(beta), '-', ENS[PIE]))) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none', 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.tag.position = c(0.3, 0.8))+
  theme(panel.grid.major = element_line(colour = "gray86", size = 0.1),
        panel.background = element_rect(fill = "white"))

fig_4b <- beta_S_PIE_all

gamma_S_PIE_all <- ggplot() +
  geom_point(data = gamma_boot_results,
             aes(x = Treatment, y = ENSPIE_mean, colour = Treatment),
             size = 4) +
  geom_errorbar(data = gamma_boot_results,
                aes(x = Treatment, ymin = ENSPIE_Q5, ymax = ENSPIE_Q95, 
                    colour = Treatment),
                size = 1.3,
                width = 0.1)+
  scale_color_manual(values =  c("#BB9689", '#836656',"#6C3859"))+
 labs(subtitle = 'c)', x = '',
       y = expression(paste(italic(gamma),'-', ENS[PIE]))) +
  theme_bw() +
  theme(legend.position = 'none', 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.tag.position = c(0.3, 0.8))

fig_4c <- gamma_S_PIE_all

Evenness <- fig_4a+fig_4b+fig_4c

# To add images to x axis
treats <- axis_canvas(Evenness,axis = 'x')+
  cowplot::draw_image('CPFP.png', x = 0.5, scale = 0.5)+
  cowplot::draw_image('CPFA.png', x = 1.5, scale = 0.5)+
  cowplot::draw_image('CAFA.png', x = 2.5, scale = 0.5)

Fig_4a <- ggdraw(insert_xaxis_grob(fig_4a, treats, position = "bottom"))
Fig_4b <- ggdraw(insert_xaxis_grob(fig_4b, treats, position = "bottom"))
Fig_4c <- ggdraw(insert_xaxis_grob(fig_4c, treats, position = "bottom"))


Evenness <- Fig_4a+Fig_4b+Fig_4c

Evenness+plot_annotation(title ="Species evenness",
                         theme = theme(plot.title = element_text(size = 14, hjust = 0.5)))



# Save image (Evenness)
ggsave('fig_4.jpg', width = 10, height = 6, dpi = 300) 
