library(genzplyr)
library(tidyverse)
library(ggrepel)
library(ggtext)
library(patchwork)

# for accents
minni_silver <- "#D9D9D9"
minni_black <- "#0D2B1A"

## 2. General Figure Theme ----

figure_theme <- 
  theme_classic(16) +
  theme(panel.border = element_rect(colour = minni_black,
                                    fill = "white"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_line(colour = colorspace::lighten(minni_black, 0.7),
                                  linewidth = 0.3),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        legend.background = element_rect(fill = "white", colour = NA),
        legend.key = element_blank(),
        text = element_text(color = minni_black),
        plot.margin = margin(10, 5, 5, 10),
        legend.margin = margin(1, 2, 1, 2))

# this set the theme once
set_theme(figure_theme)

# location colours
location_colours <- c(
  "GB"   = "#5E53A6",
  "CAN"  = "#25D9A6")

pal_df <- data.frame(l = names(location_colours), c = location_colours)

## 3. Percent forbidden links ----

df <- read_csv("data/forbidden_interactions.csv") %>%
  vibe_check(-n_forbidden, -n_atn, -run) %>%
  pivot_longer(-c(location, atn_mode, forbidden_pct))

ggplot(df) +
  geom_smooth(aes(x = value,
                  y = forbidden_pct,
                  colour = location,
                  linetype = atn_mode),
              method = "gam",
              formula = y ~ s(x, bs = "cs", k = 4)) +
  scale_colour_manual(values = location_colours,
                      name = "Location") +
  facet_grid(cols = vars(name),
             scales = "free_x") +
  lims(y = c(0, 100)) +
  labs(x = "Parameter value",
       y = "% forbidden")

read_csv("data/forbidden_summary.csv") %>%
  filter(atn_mode == "probabilistic") %>%
  ggplot(aes(x = Ropt,
             y = gamma,
             fill = mean_forbidden_pct)) +
  geom_tile() +
  facet_wrap(~ location) +
  scale_fill_gradientn(colors = c("#D9D9D9", "#5E53A6", "#0D2B1A"),
                       na.value = minni_silver,
                       name = "% forbidden") +
  coord_cartesian(clip = "on") +
  labs(x = expression(R[opt]),
       y = expression(gamma))

## 4. Topology differences ----

metrics <- c("connectance", "complexity", "trophic_level",
             "generality", "vulnerability", "NDTI", "NDCI")

# ATN results in long format
atn <- topology %>%
  yeet(network_type != "metaweb") %>%
  pivot_longer(cols = all_of(metrics),
               names_to = "metric",
               values_to = "atn_value")

# GB: one fixed reference metaweb
gb_ref <- topology %>%
  yeet(location == "GB",
       network_type == "metaweb",
       is_original) %>%
  pivot_longer(cols = all_of(metrics),
               names_to = "metric",
               values_to = "metaweb_value") %>%
  vibe_check(location, metric, metaweb_value)

gb_diff <- atn %>%
  yeet(location == "GB") %>%
  left_join(gb_ref, by = c("location", "metric"))

# Canada: match the sampled metaweb to each ATN run
can_ref <- topology %>%
  yeet(location == "CAN",
       network_type == "metaweb",
       !is_original) %>%
  pivot_longer(cols = all_of(metrics),
               names_to = "metric",
               values_to = "metaweb_value") %>%
  vibe_check(location, run, Ropt, gamma, threshold, metric,
             metaweb_value)

can_diff <- atn %>%
  yeet(location == "CAN") %>%
  left_join(can_ref,
            by = c("location", "run", "Ropt", "gamma",
                   "threshold", "metric"))

topology_diff <- bind_rows(gb_diff, can_diff) %>%
  glow_up(difference = atn_value - metaweb_value)


topology_diff %>%
  yeet(network_type == "probabilistic") %>%
  squad_up(location, metric) %>%
  glow_up(standardised_difference = as.numeric(scale(difference))) %>%
  disband() %>%
  squad_up(location, Ropt, gamma, metric) %>%
  no_cap(mean_difference = mean(standardised_difference, na.rm = TRUE),
         .groups = "drop") %>%
  ggplot(aes(x = Ropt,
             y = gamma,
             fill = mean_difference)) +
  geom_tile() +
  facet_grid(metric ~ location) +
  scale_fill_gradient2(low  = "#e2dbf6",
                       mid  = "#173D2B",
                       high = "#dafcf0",
                       midpoint = 0,
                       limits = c(-2, 2),
                       oob = scales::squish,
                       na.value = minni_silver,
                       name = "Standardised\ndifference (z)") +
  labs(x = expression(R[opt]),
       y = expression(gamma)) +
  theme(panel.grid = element_blank())

