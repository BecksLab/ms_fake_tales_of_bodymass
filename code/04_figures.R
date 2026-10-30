library(genzplyr)
library(tidyverse)
library(ggrepel)
library(ggtext)
library(patchwork)

# for accents
minni_silver <- "#D9D9D9"
minni_black <- "#0D2B1A"
minni_white <- "#D1D1D1"

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
  pivot_longer(-c(location, atn_mode, forbidden_pct)) %>%
  glow_up(name = case_when(name == "Ropt"      ~ "R[opt]",
                           name == "gamma"     ~ "gamma",
                           name == "threshold" ~ "Threshold",))

# Panel A: sensitivity to individual ATN parameters
p_forbidden_trends <- 
  ggplot(df) +
  geom_smooth(aes(x = value, 
                  y = forbidden_pct,
                  colour = location, linetype = atn_mode,
                  group = interaction(location, atn_mode)),
              method = "gam",
              formula = y ~ s(x, bs = "cs", k = 4)) +
  scale_colour_manual(values = location_colours,
                      name = "Location") +
  facet_grid(cols = vars(name),
             scales = "free_x",
             labeller = labeller(name = label_parsed)) +
  scale_y_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 25)) +
  labs(x = "Parameter value", 
       y = "% forbidden links")

# Panel B: joint sensitivity of probabilistic ATNs
p_forbidden_heatmap <- 
  read_csv("data/forbidden_summary.csv") %>%
  yeet(atn_mode == "probabilistic") %>%
  ggplot(aes(x = Ropt, 
             y = gamma, 
             fill = mean_forbidden_pct)) +
  geom_tile() +
  facet_wrap(~location) +
  scale_fill_gradientn(colours = c(minni_silver, "#5E53A6", minni_black),
                       #limits = c(0, 100),
                       name = "% forbidden") +
  labs(x = expression(R[opt]), 
       y = expression(gamma))

# Combine and export
fig1 <- p_forbidden_trends / p_forbidden_heatmap +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold"))

ggsave("../figures/forbidden_links.png",
       fig1, 
       width = 10, 
       height = 8)

## 4. Topology differences ----

topology <- read_csv("data/topology.csv")

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


topology_diff_std <- topology_diff %>%
  squad_up(location,
           network_type,
           Ropt,
           gamma,
           threshold,
           metric) %>%
  glow_up(mean_difference = mean(difference, na.rm = TRUE),
          sd_difference = sd(difference, na.rm = TRUE),
          n = sum(is.finite(difference)), 
          standardised_difference = if_else(is.finite(sd_difference) & sd_difference > 0,
                                            mean_difference / sd_difference,
                                            NA_real_)) %>%
  disband() %>%
  yeet(network_type == "probabilistic") %>%
  squad_up(location, Ropt, gamma, metric) %>%
  no_cap(mean_difference = mean(standardised_difference, na.rm = TRUE),
         .groups = "drop")

p_structure <- topology_diff_std %>%
  ggplot(aes(x = Ropt, y = gamma, fill = mean_difference)) +
  geom_tile() +
  facet_grid(metric ~ location) +
  scale_fill_gradient2(low = "#381E40",
                       mid = minni_silver,
                       high = "#173D2B",
                       midpoint = 0,
                       limits = c(min(topology_diff_std$mean_difference, na.rm = TRUE),
                                  max(topology_diff_std$mean_difference, na.rm = TRUE)),
                       name = "Standardised\ndifference",
                       na.value = "white") +
  labs(x = expression(R[opt]), y = expression(gamma)) +
  theme(panel.grid = element_blank())

ggsave("../figures/strucutre.png",
       p_structure, 
       width = 10, 
       height = 13)

# Lower values indicate closer structural matches
structure_score <- topology_diff_std %>%
  filter(network_type == "probabilistic") %>%
  group_by(location, Ropt, gamma) %>%
  summarise(
    structural_discrepancy = mean(
      abs(standardised_difference), na.rm = TRUE
    ),
    .groups = "drop"
  )

agreement <- read_csv("data/interaction_agreement.csv") %>%
  filter(atn_mode == "probabilistic") %>%
  group_by(location, Ropt, gamma) %>%
  summarise(
    jaccard = mean(jaccard, na.rm = TRUE),
    .groups = "drop"
  )

agreement %>%
  ggplot(aes(x = Ropt, 
             y = gamma, 
             fill = jaccard)) +
  geom_tile() +
  facet_wrap(~location) +
  scale_fill_gradientn(colours = c(minni_silver, "#5E53A6", minni_black),
                       #limits = c(0, 100),
                       name = "Jaccard similarity") +
  labs(x = expression(R[opt]), 
       y = expression(gamma))

p_agreement <- topology_diff_std %>%
  left_join(agreement, by = c("location", "Ropt", "gamma")) %>%
  ggplot(aes(x = standardised_difference,
             y = jaccard,
             colour = location)) +
  geom_point(alpha = 0.5, 
             size = 2) +
  scale_colour_manual(values = location_colours, name = "Location") +
  geom_vline(xintercept = 0,
             colour = minni_silver) +
  facet_grid(metric ~ location) +
  labs(x = "Standardised structural difference",
       y = "Jaccard interaction similarity") +
  theme(legend.position = "none")
