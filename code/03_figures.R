library(genzplyr)
library(tidyverse)
library(ggrepel)
library(ggtext)
library(patchwork)

# for accents
minni_silver <- "#D9D9D9"
minni_black <- "#0D2B1A"

## 2. General Figure Theme ----

figure_theme = 
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

# 1. DISCRETE MODEL COLORS
location_colours <- c(
  "GB"       = "#5E53A6",
  "Canada"   = "#25D9A6")

pal_df <- data.frame(l = names(location_colours), c = location_colours)

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
  labs(x = expression(R[opt]),
       y = expression(gamma))
