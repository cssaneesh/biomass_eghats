library(tidyverse)
library(betapart)
library(brms)
library(bayesplot)
library(patchwork)
library(cmdstanr)

transect_dat <-
  read.csv(
    "biomass_data.csv",
    header = T,
    fill = TRUE,
    sep = ",",
    na.strings = c("", " ", "NA", "NA ", "na", "NULL")
  )

head(transect_dat)

transect_weight <- transect_dat %>%
  group_by(Site, Transect, Treatment) %>%
  summarise(Plot_Weight = sum(Weight))


transect_calc <- transect_dat %>%
  left_join(transect_weight) %>%
  mutate(relative_weight = (Weight / Plot_Weight))


head(transect_calc)

transect_details <- transect_calc %>%
  select(Transect, Treatment, Site) %>% distinct()  %>%
  mutate(
    Treatment = case_when(
      Treatment == "ab" ~ "Control",
      # Cymbopogon present fire present
      Treatment == "bgpnf" ~ "CPFA",
      # Cymbopogon present fire absent
      Treatment == "bgrnf" ~ "CAFA" # Cymbopogon absent fire absent
    )
  ) %>%
  mutate(Treatment = factor(Treatment)) %>%
  mutate(Treatment = fct_relevel(Treatment, c("Control", "CPFA", "CAFA")))  %>%
  ungroup() %>%
  arrange(Site, Treatment)

head(transect_details, 3)

control_prep <- transect_calc %>% as_tibble() %>%
  mutate(pres = 1) %>%
  select(Transect, Treatment, Sci_name, pres) %>%
  mutate(
    Treatment = case_when(
      Treatment == "ab" ~ "Control",
      # Cymbopogon present fire present
      Treatment == "bgpnf" ~ "CPFA",
      # Cymbopogon present fire absent
      Treatment == "bgrnf" ~ "CAFA" # Cymbopogon absent fire absent
    )
  ) %>%
  mutate(Treatment = factor(Treatment)) %>%
  filter(Treatment == "Control") %>%
  select(-Transect) %>%
  distinct() %>%
  mutate(Sci_name_sp = paste0('sp_', Sci_name)) %>%
  select(-Sci_name) %>%
  group_by(Treatment, Sci_name_sp) %>%
  spread(Sci_name_sp, pres, fill = 0)

head(control_prep)

ctls <- transect_details %>% select(Transect) %>% distinct() %>%
  bind_cols(control_prep)

head(ctls, 3)

species_prep <- transect_calc %>% as_tibble() %>%
  mutate(pres = 1) %>%
  select(Transect, Treatment, Sci_name, pres) %>%
  mutate(
    Treatment = case_when(
      Treatment == "ab" ~ "Control",
      # Cymbopogon present fire present
      Treatment == "bgpnf" ~ "CPFA",
      # Cymbopogon present fire absent
      Treatment == "bgrnf" ~ "CAFA" # Cymbopogon absent fire absent
    )
  ) %>%
  mutate(Treatment = factor(Treatment)) %>%
  mutate(Treatment = fct_relevel(Treatment, c("Control", "CPFA", "CAFA")))  %>%
  ungroup()

View(species_prep)

species_wide <- species_prep %>%
  mutate(Sci_name_sp = paste0('sp_', Sci_name)) %>%
  select(-Sci_name) %>%
  group_by(Transect, Sci_name_sp) %>%
  spread(Sci_name_sp, pres, fill = 0) %>%
  bind_rows(ctls) %>% arrange(Transect) %>%
  replace(is.na(.), 0) # %>% View()

head(species_wide)

beta_pairs <- function(x) {
  # function to return the dissimilarities (turnover and nestedness component)
  # for each control treatment comparison in x
  # return dataframe with dissimilarities and the treatment magnitude (seed.rich)
  # separate out the control and treatment plots
  contr.plots = x %>%
    filter(Treatment == "Control")
  # fix for treatment labels
  trt.plots = x %>%
    filter(Treatment == 'CPFA' |
             Treatment == 'CAFA')
  out <- tibble()
  if (nrow(contr.plots) > 0) {
    for (i in 1:nrow(contr.plots)) {
      beta = beta.pair(
        bind_rows(
          contr.plots %>%
            slice(i) %>%
            select(-Treatment),
          trt.plots %>%
            select(-Treatment)
        ),
        index.family = 'jaccard'
      )
      # build the data we want for analysis
      out <- bind_rows(
        out,
        tibble(
          Treatment = trt.plots$Treatment,
          jtu = as.matrix(beta$beta.jtu)[-1, 1],
          jne = as.matrix(beta$beta.jne)[-1, 1],
          group = i
        )
      )
    }
  }
  # escape for no controls
  else{
    out = tibble(
      Treatment = NA,
      jtu = NA,
      jne = NA,
      group = NA
    )
  }
  return(out)
}

head(transect_details)

wide.df <- species_wide %>%
  left_join(transect_details) %>%
  group_by(Transect) %>%
  nest_legacy(starts_with('sp_'), Treatment)

head(wide.df)

wide.df <- wide.df %>%
  mutate(beta = purrr::map(data, ~ beta_pairs(.x))) # css


head(wide.df)

beta.df <- wide.df %>%
  unnest_legacy(beta) %>%
  unite(col = pw_beta_group,
        c(Treatment),
        sep = '_',
        remove = F) %>%
  select(-group) %>% left_join(transect_details)


# write.csv(beta.df, "beta.csv")

# box plots
ggplot() +
  geom_boxplot(data = beta.df,
               aes(x = Treatment, y = jtu)) +
  labs(y = "Turnover")


ggplot() +
  geom_boxplot(data = beta.df,
               aes(x = Treatment, y = jne)) +
  labs(y = "Nestedness")


head(beta.df)


# models
 ghats.turnover <-
  brm(
    jtu ~   Treatment + (1 | Site) ,
    family = zero_one_inflated_beta(),
    data = beta.df,
    iter = 3000,
    warmup = 1000,
    cores = 4,
    chains = 4,
    control = list(adapt_delta = 0.99)
  )

  save(ghats.turnover, file = 'ghats.turnover.Rdata')
load('ghats.turnover.Rdata')

summary(ghats.turnover)

ghats.nest <- brm(
  jne ~  Treatment + (1 | Site) ,
  family = zero_inflated_beta(),
  data = beta.df,
  iter = 3000,
  warmup = 1000,
  cores = 4,
  chains = 4,
  control = list(adapt_delta = 0.99)
)

save(ghats.nest , file = 'ghats.nest .Rdata')
load('ghats.nest .Rdata')

summary(ghats.turnover) # model summary


color_scheme_set("darkgray")
fig_s4a <- pp_check(ghats.turnover) +
  xlab("Turnover") + ylab("Density") +
  labs(subtitle = "a)") +
  theme_classic() +  theme(plot.title = element_text(size = 18, hjust =
                                                       0.5),
                           legend.position = "none")# predicted vs. observed values

fig_s4a


ghats_t <-
  conditional_effects(
    ghats.turnover,
    effects = 'Treatment',
    re_formula = NA,
    method = 'fitted'
  )  # conditional effects

head(ghats_t)



#head(alpha_c)

fig_4a <- ggplot() +
  geom_point(
    data = beta.df,
    aes(x = Treatment, y = jtu, colour = 	"#C0C0C0"),
    size = 1,
    alpha = 0.2,
    position = position_jitter(width = 0.05, height = 0.45)
  ) +
  geom_point(
    data = ghats_t$Treatment,
    aes(x = Treatment, y = estimate__, colour = Treatment),
    size = 3
  ) +
  geom_errorbar(
    data = ghats_t$Treatment,
    aes(
      x = Treatment,
      ymin = lower__,
      ymax = upper__,
      colour = Treatment
    ),
    size = 1,
    width = 0
  ) +
  labs(x = '',
       y = '') +
  scale_color_manual(values =  c("#C0C0C0", "#228B22", 	"#6B8E23"))  +
  theme_bw(base_size = 18) + theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(
      t = 0.2,
      r = 0.2,
      b = -0.2,
      l = 0.2,
      unit = "cm"
    ),
    plot.title = element_text(size = 18, hjust =
                                0.5),
    strip.background = element_blank(),
    legend.position = "none"
  ) +
  labs(subtitle = 'a)') + ylab("Turnover")


fig_4a




summary(ghats.nest) # model summary



fig_s4b <- pp_check(ghats.nest) +
  xlab("Nestedness") + ylab("Density") +
  labs(subtitle = "a)") +
  theme_classic() +  theme(plot.title = element_text(size = 18, hjust =
                                                       0.5),
                           legend.position = "none")# predicted vs. observed values

fig_s4b


ghats_n <-
  conditional_effects(
    ghats.nest,
    effects = 'Treatment',
    re_formula = NA,
    method = 'fitted'
  )  # conditional effects

head(ghats_n)


#head(alpha_c)

fig_4b <- ggplot() +
  geom_point(
    data = beta.df,
    aes(x = Treatment, y = jne, colour = 	"#C0C0C0"),
    size = 1,
    alpha = 0.2,
    position = position_jitter(width = 0.05, height = 0.45)
  ) +
  geom_point(
    data = ghats_n$Treatment,
    aes(x = Treatment, y = estimate__, colour = Treatment),
    size = 3
  ) +
  geom_errorbar(
    data = ghats_n$Treatment,
    aes(
      x = Treatment,
      ymin = lower__,
      ymax = upper__,
      colour = Treatment
    ),
    size = 1,
    width = 0
  ) +
  labs(x = '',
       y = '') +
  scale_color_manual(values =  c("#C0C0C0", "#228B22", 	"#6B8E23"))  +
  theme_bw(base_size = 18) + theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(
      t = 0.2,
      r = 0.2,
      b = -0.2,
      l = 0.2,
      unit = "cm"
    ),
    plot.title = element_text(size = 18, hjust =
                                0.5),
    strip.background = element_blank(),
    legend.position = "none"
  ) +
  labs(subtitle = 'a)') + ylab("Nestedness")


fig_4b


(fig_4a |fig_4b)
